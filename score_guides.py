"""
Score a SPECIFIC list of guides against a target locus.

WHY THIS EXISTS
---------------
run_pipeline.py discovers its own guides and ranks them. That is the right
behaviour when you are starting from scratch, but it is the wrong shape for
a very common workflow:

    1. Design + screen genome-wide in a tool that has a real genome index
       (CHOPCHOP, CRISPOR, Cas-OFFinder). You get back a shortlist of guides
       that are known-clean across the whole genome.
    2. Bring that shortlist here to get the newer Rule Set 3 on-target score,
       a real CFD specificity number, and a ready-to-order HDR donor.

Without this script you would have to eyeball two tables and hope the
sequences line up. This script closes that gap: give it your target FASTA
and the guides someone else recommended, and it scores exactly those.

INPUT FORMATS FOR --guides-file
-------------------------------
- A plain text file, one protospacer per line.
- A CHOPCHOP results TSV/CSV: any column named "Target sequence" (case
  insensitive) is used automatically. CHOPCHOP reports 23 nt including the
  PAM; the trailing PAM is stripped for you.
Lines starting with '#' are ignored in both cases.

WHAT IT DOES NOT DO
-------------------
It does not re-derive genome-wide off-target counts — that is precisely the
thing you used the other tool for. Off-target screening here is still
reference-scoped (see offtarget.py). Carry the genome-wide numbers across
yourself; there is a column for them.
"""

from __future__ import annotations
import argparse
import csv
import sys

from sgrna_design import find_guides, score_guide, revcomp
from knockin_donor import design_knockin_donor
from offtarget import (
    scan_reference, summarize as summarize_offtargets,
    build_index, scan_with_index, read_masked_fasta_records,
)
import rule_set_3
from run_pipeline import read_fasta, read_fasta_records


def _clean(seq: str) -> str:
    return "".join(c for c in seq.strip().upper() if c in "ACGTN")


def _looks_like_protospacer(s: str, guide_len: int = 20) -> bool:
    """A 20 nt (or 23 nt with PAM) pure-ACGT string."""
    s = _clean(s)
    return len(s) in (guide_len, guide_len + 3) and set(s) <= set("ACGT")


def read_guides_file(path: str, verbose: bool = True) -> list[str]:
    """
    Accepts:
      - a plain list, one protospacer per line;
      - any delimited table (CHOPCHOP, CRISPOR, Benchling, a spreadsheet export).

    For tables the guide column is found by NAME when the header is one we
    recognise, and otherwise BY CONTENT: the column whose values most often
    look like a protospacer wins. Detecting by content matters because
    different tools label this column differently and not all of them
    document their exact headers, so name matching alone is brittle.
    """
    with open(path) as fh:
        text = fh.read()

    lines = [ln for ln in text.splitlines() if ln.strip() and not ln.lstrip().startswith("#")]
    if not lines:
        return []

    for delim in ("\t", ","):
        if delim not in lines[0]:
            continue
        rows = list(csv.reader(lines, delimiter=delim))
        if len(rows) < 2:
            continue
        header = [h.strip().lower() for h in rows[0]]
        body = rows[1:]

        # 1. Recognised header name.
        col = None
        for wanted in ("target sequence", "target_sequence", "protospacer",
                        "guide sequence", "guide_sequence", "seq", "sequence", "guide"):
            if wanted in header:
                col = header.index(wanted)
                break

        # 2. Otherwise pick the column that most looks like guide sequences.
        if col is None:
            best, best_hits = None, 0
            ncols = max(len(r) for r in body)
            for c in range(ncols):
                hits = sum(1 for r in body
                            if len(r) > c and _looks_like_protospacer(r[c]))
                if hits > best_hits:
                    best, best_hits = c, hits
            # require it to work for at least half the rows
            if best is not None and best_hits >= max(1, len(body) // 2):
                col = best
                if verbose:
                    label = header[col] if col < len(header) else f"column {col + 1}"
                    print(f"[guide column auto-detected by content: {label!r}]")

        if col is None:
            continue

        out = []
        for r in body:
            if len(r) > col:
                s = _clean(r[col])
                if s:
                    out.append(s)
        if out:
            return out

    # Otherwise: one sequence per line.
    return [s for s in (_clean(ln) for ln in lines) if s]


def normalize_guide(seq: str, guide_len: int = 20) -> str:
    """CHOPCHOP and friends often report protospacer+PAM (23 nt). Trim to 20."""
    seq = _clean(seq)
    if len(seq) == guide_len + 3:
        return seq[:guide_len]
    return seq


def locate(guides, protospacer: str):
    """Find the Guide object matching this protospacer, on either strand."""
    for g in guides:
        if g.protospacer == protospacer:
            return g
    return None


def main():
    p = argparse.ArgumentParser(
        description="Score a specific list of guides (e.g. a CHOPCHOP shortlist) "
                     "against a target locus: real Rule Set 3 on-target + real CFD specificity.")
    p.add_argument("fasta", help="Target locus FASTA (single record).")
    p.add_argument("--guides", nargs="+", default=[], metavar="SEQ",
                    help="Protospacer sequences, 20 nt (or 23 nt with PAM — the PAM is trimmed).")
    p.add_argument("--guides-file", metavar="FILE",
                    help="File of guides: one per line, or a CHOPCHOP TSV/CSV with a "
                         "'Target sequence' column.")
    p.add_argument("--reference", "-r", nargs="+", default=[], metavar="FASTA",
                    help="Extra FASTA files to off-target-screen against (reference-scoped).")
    p.add_argument("--max-mismatches", type=int, default=6)
    p.add_argument("--seed-len", type=int, default=12)
    p.add_argument("--fast-index-threshold", type=int, default=2_000_000)
    p.add_argument("--build-donor", metavar="SEQ",
                    help="After scoring, build the HDR donor for this guide.")
    p.add_argument("--cargo-fasta", metavar="FASTA", help="Cargo to insert (for --build-donor).")
    p.add_argument("--cargo-name", default="cargo")
    p.add_argument("--arm-len", type=int, default=40)
    args = p.parse_args()

    wanted = [normalize_guide(s) for s in args.guides]
    if args.guides_file:
        wanted += [normalize_guide(s) for s in read_guides_file(args.guides_file)]
    # de-duplicate, preserve order
    seen, ordered = set(), []
    for w in wanted:
        if w and w not in seen:
            seen.add(w)
            ordered.append(w)
    if not ordered:
        p.error("No guides supplied. Use --guides and/or --guides-file.")

    header, seq, locus_mask = read_masked_fasta_records(args.fasta)[0]
    print("=" * 74)
    print(f"SCORING {len(ordered)} SUPPLIED GUIDE(S)  |  target: {header}")
    print("=" * 74)
    print(f"Locus length: {len(seq)} bp")
    print()

    # No CDS declared here: this script scores a supplied shortlist, and the
    # knockout position rule needs a real CDS window to mean anything.
    all_guides = [score_guide(g) for g in find_guides(seq)]

    found, missing = [], []
    for w in ordered:
        g = locate(all_guides, w)
        if g is None:
            missing.append(w)
        else:
            found.append(g)

    if missing:
        print(f"NOT FOUND in this locus ({len(missing)}):")
        for m in missing:
            rc = revcomp(m)
            hint = ""
            if m in seq or rc in seq:
                hint = "  (present in the sequence, but not adjacent to an NGG PAM)"
            print(f"  {m}{hint}")
        print("  -> Check you supplied the right locus, and that the guide is "
              "written 5'->3' in its own orientation.")
        print()
    if not found:
        sys.exit("None of the supplied guides were found in this locus.")

    # --- real Rule Set 3, batched ---
    using_rs3 = False
    for g in found:
        g.rs3_score = None
    if rule_set_3.available():
        scorable = [g for g in found if rule_set_3.context_30mer(seq, g) is not None]
        ctxs = [rule_set_3.context_30mer(seq, g) for g in scorable]
        if ctxs:
            for g, s in zip(scorable, rule_set_3.predict(ctxs)):
                g.rs3_score = s
            using_rs3 = True
            print(f"On-target: REAL Rule Set 3 (Doench lab 2022, {rule_set_3.mode()})")
        else:
            print("On-target: heuristic (no guide had enough flanking for a 30 nt context)")
    else:
        print("On-target: heuristic (real Rule Set 3 not available — see rule_set_3.py)")

    # --- off-target screen ---
    refs = [(header, seq, locus_mask)]
    for path in args.reference:
        refs.extend((f"{path}:{h}" if h else path, s, m)
                     for h, s, m in read_masked_fasta_records(path))
    print(f"Off-target: screening against {len(refs)} reference sequence(s)"
          f"{' — locus only' if len(refs) == 1 else ''}")
    print()

    hits_by = [[] for _ in found]
    for ri, (rname, rseq, rmask) in enumerate(refs):
        is_locus = ri == 0
        fast = len(rseq) > args.fast_index_threshold
        idx = build_index(rseq, ref_name=rname, seed_len=args.seed_len, mask=rmask) if fast else None
        for gi, g in enumerate(found):
            on_pos = g.guide_start_fwd if is_locus else None
            on_str = g.strand if is_locus else None
            if fast:
                hits_by[gi].extend(scan_with_index(
                    g.protospacer, idx, max_mismatches=args.max_mismatches,
                    on_target_position=on_pos, on_target_strand=on_str))
            else:
                hits_by[gi].extend(scan_reference(
                    g.protospacer, rseq, ref_name=rname,
                    max_mismatches=args.max_mismatches, seed_len=args.seed_len,
                    on_target_position=on_pos, on_target_strand=on_str, mask=rmask))

    summaries = [summarize_offtargets(h) for h in hits_by]
    rows = list(zip(found, summaries))
    rank_key = (lambda t: (t[0].rs3_score is not None, t[0].rs3_score)) if using_rs3 \
        else (lambda t: t[0].score)
    rows.sort(key=rank_key, reverse=True)

    rank_hdr = "RS3" if using_rs3 else "score"
    print(f"{'#':>2}  {'protospacer':<22} {'PAM':<4} {'str':<4} {'cut':>6} {'GC%':>5} "
          f"{rank_hdr:>7} {'CFD-spec':>8} {'exact-OT':>9}")
    for i, (g, sm) in enumerate(rows, 1):
        rv = f"{g.rs3_score:>7.2f}" if using_rs3 and g.rs3_score is not None else f"{g.score:>7.0f}"
        print(f"{i:>2}  {g.protospacer:<22} {g.pam:<4} {g.strand:<4} {g.cut_site_fwd:>6} "
              f"{g.gc:>5.0f} {rv} {sm['cfd_specificity_score']:>8.1f} {sm['exact_offtargets']:>9}")

    print()
    print("Reminder: exact-OT here counts only what is in the sequences you supplied. "
          "Genome-wide counts (e.g. CHOPCHOP MM0-MM3) are not reproduced — carry them across yourself.")

    # --- optional donor ---
    if args.build_donor:
        target = normalize_guide(args.build_donor)
        g = locate(found, target)
        if g is None:
            sys.exit(f"\n--build-donor {target} is not among the scored guides.")
        if args.cargo_fasta:
            _, cargo = read_fasta(args.cargo_fasta)
        else:
            sys.exit("\n--build-donor requires --cargo-fasta.")
        print()
        print("-" * 74)
        print(f"HDR DONOR for {g.protospacer} (cut {g.cut_site_fwd}, {args.arm_len} bp arms)")
        print("-" * 74)
        donor = design_knockin_donor(
            seq, cut_site=g.cut_site_fwd, insert_seq=cargo,
            insert_name=args.cargo_name, protospacer=g.protospacer,
            pam_pos=g.guide_end_fwd if g.strand == "+" else None,
            arm_len=args.arm_len)
        print(donor.annotated_map())
        print()
        print(donor.ascii_track())
        print()
        print("Donor sequence:")
        print(donor.full_donor)


if __name__ == "__main__":
    main()
