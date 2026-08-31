"""
End-to-end CRISPR knock-in pipeline runner.

  FASTA/sequence -> find guides -> off-target screen -> pick best -> build donor -> annotated map

Usage:
    python3 run_pipeline.py                          # runs the built-in demo locus
    python3 run_pipeline.py my_locus.fa               # runs on your FASTA
    python3 run_pipeline.py my_locus.fa \\
        --reference paralogs.fa pseudogenes.fa        # + off-target screen vs a reference panel

Run `python3 run_pipeline.py --help` for all options.

Drop in a real AAVS1 (or any) sequence and the same code produces the real map.
"""

from __future__ import annotations
import argparse
from sgrna_design import find_guides, score_guide, revcomp, cds_fraction
from knockin_donor import design_knockin_donor
from offtarget import (
    scan_reference, summarize as summarize_offtargets,
    build_index, scan_with_index, read_masked_fasta_records,
)
import rule_set_3


# ---------- FASTA reading ----------
def read_fasta_records(path: str) -> list[tuple[str, str]]:
    """Return every (header, sequence) record in a FASTA file."""
    records: list[tuple[str, str]] = []
    header, parts = None, []
    with open(path) as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    records.append((header, "".join(parts)))
                header, parts = line[1:], []
            else:
                parts.append(line.upper())
    if header is not None:
        records.append((header, "".join(parts)))
    for h, seq in records:
        bad = set(seq) - set("ACGTN")
        if bad:
            raise ValueError(f"Non-DNA characters in FASTA record {h!r}: {sorted(bad)}")
    return records


def read_fasta(path: str) -> tuple[str, str]:
    """Return the single record from a FASTA holding one target locus."""
    records = read_fasta_records(path)
    if not records:
        raise ValueError(f"No sequence records found in {path}")
    if len(records) > 1:
        raise ValueError(
            f"{path} has {len(records)} records but run_pipeline expects a single "
            f"target locus. Use --reference for multi-sequence off-target panels."
        )
    return records[0]


# ---------- The pipeline ----------
def run(seq: str, cargo: str, cargo_name: str = "CAR",
        arm_len: int = 40, top_n: int = 5, label: str = "locus",
        reference_seqs: list[tuple[str, str, list | None]] | None = None,
        locus_mask: list[bool] | None = None,
        screen_n: int = 15, max_mismatches: int = 6, seed_len: int = 12,
        allow_exact_repeats: bool = False, fast_index_threshold: int = 2_000_000,
        mode: str = "knockin", cds_start: int | None = None,
        cds_end: int | None = None, cds_upstream_bp: int = 0,
        cds_total_len: int | None = None):

    knockout = (mode == "knockout")
    have_cds = cds_start is not None and cds_end is not None

    def cds_frac(g):
        return cds_fraction(g, cds_start, cds_end, cds_upstream_bp, cds_total_len)

    print("=" * 70)
    print(f"CRISPR {'KNOCK-OUT' if knockout else 'KNOCK-IN'} PIPELINE  |  target: {label}")
    print("=" * 70)
    print(f"Input sequence: {len(seq)} bp")
    if knockout:
        if have_cds:
            _total = cds_total_len if cds_total_len else (cds_end - cds_start)
            print(f"Coding sequence: offsets {cds_start}-{cds_end} in this window "
                  f"({cds_end - cds_start} bp)")
            if cds_total_len or cds_upstream_bp:
                print(f"                 window covers coding bases "
                      f"{cds_upstream_bp}-{cds_upstream_bp + (cds_end - cds_start)} "
                      f"of a {_total} bp full CDS")
            else:
                print(f"                 NOTE: treated as the COMPLETE CDS. If this "
                      f"window is one exon of a multi-exon gene, pass "
                      f"--cds-total-len (and --cds-upstream-bp) or the "
                      f"position-in-CDS figures will be wrong.")
        else:
            print("Coding sequence: NOT DECLARED — see the warning in STEP 1")
        print("Repair pathway: NHEJ (no donor is designed for a knock-out)")
    else:
        print(f"Cargo to insert: {cargo_name} ({len(cargo)} bp)")
    print()

    # --- STEP 1: find guides on both strands ---
    print("-" * 70)
    print("STEP 1  Find candidate guides (SpCas9, NGG PAM, both strands)")
    print("-" * 70)
    guides = [score_guide(g, cds_start, cds_end) for g in find_guides(seq)]
    if not guides:
        print("No guides found. Is the sequence long enough / does it contain NGG?")
        return
    guides.sort(key=lambda g: g.score, reverse=True)
    print(f"Found {len(guides)} candidate guides.")

    if knockout:
        # A knock-out needs no homology arms, so the arm-length filter does not
        # apply. What matters instead is landing inside the coding sequence.
        if have_cds:
            usable = [g for g in guides
                      if cds_frac(g) is not None]
            dropped = len(guides) - len(usable)
            if dropped:
                print(f"Dropped {dropped} guide(s) cutting outside the declared "
                      f"CDS ({cds_start}-{cds_end}).")
            if not usable:
                print("No guide cuts inside the declared CDS. Check your "
                      "--cds-start/--cds-end offsets.")
                return
            guides = usable
        else:
            print("WARNING: no --cds-start/--cds-end given, so guides cannot be "
                  "ranked by position within the coding sequence, and none can be "
                  "filtered out for landing in an intron or UTR. Supply the CDS "
                  "offsets for a meaningful knock-out design.")
    else:
        # Guides too close to either end can't support full-length homology arms.
        usable = [g for g in guides
                  if g.cut_site_fwd >= arm_len and g.cut_site_fwd <= len(seq) - arm_len]
        dropped = len(guides) - len(usable)
        if dropped:
            print(f"Dropped {dropped} guide(s) too close to the sequence ends "
                  f"(need {arm_len} bp flanking for full homology arms).")
        if not usable:
            print(f"No guide leaves room for {arm_len} bp arms. "
                  f"Provide more flanking sequence or reduce arm_len.")
            return
        guides = usable

    # Prefer the REAL, published Rule Set 3 on-target score (Doench lab,
    # 2022) over the from-scratch heuristic when the optional rs3_venv is
    # set up (see rule_set_3.py) — one batched subprocess call scores
    # every usable guide, then re-ranks by it.
    for g in guides:
        g.rs3_score = None
    using_rs3 = False
    if rule_set_3.available():
        scorable = [g for g in guides if rule_set_3.context_30mer(seq, g) is not None]
        contexts = [rule_set_3.context_30mer(seq, g) for g in scorable]
        if contexts:
            for g, s in zip(scorable, rule_set_3.predict(contexts)):
                g.rs3_score = s
            guides.sort(key=lambda g: (g.rs3_score is not None, g.rs3_score), reverse=True)
            using_rs3 = True
            print(f"{len(guides)} usable — ranked by the REAL Rule Set 3 score "
                  f"(Doench lab, 2022, {rule_set_3.mode()}; {len(scorable)}/{len(guides)} "
                  f"had enough flanking for a 30nt context).")
        else:
            print(f"{len(guides)} usable — real Rule Set 3 available ({rule_set_3.mode()}) "
                  f"but no guide had enough flanking for a 30nt context; "
                  f"using the on-target heuristic.")
    else:
        print(f"{len(guides)} usable — ranked by the on-target heuristic "
              f"(real Rule Set 3 not set up; see rule_set_3.py to enable it).")

    # --- STEP 2: off-target screen the best on-target candidates ---
    print()
    print("-" * 70)
    print("STEP 2  Off-target screen top candidates")
    print("-" * 70)
    screen_pool = guides[:max(screen_n, 1)]
    references = [(label, seq, locus_mask)] + list(reference_seqs or [])
    ref_names = ", ".join(n for n, _, _ in references)
    print(f"Screening top {len(screen_pool)} on-target candidate(s) against "
          f"{len(references)} reference sequence(s): {ref_names}")
    if not reference_seqs:
        print("(No --reference supplied; screening against the input locus only. "
              "This catches tandem repeats/duplications but NOT genome-wide off-targets — "
              "pass --reference for paralogs/pseudogenes/a chromosome for a real screen.)")

    # Build each reference's search structure ONCE and reuse it for every
    # screened candidate, instead of re-scanning per guide. For references
    # above fast_index_threshold (real chromosome-scale sequence) this
    # switches to a seed-indexed lookup (build once, O(1)-ish per guide) so
    # screening 10-15 candidates against a multi-megabase human reference
    # stays tractable; the exact-off-target safety gate is unaffected by
    # this (see offtarget.py), only lower-risk seed-mismatched hits are
    # left un-enumerated at that scale.
    hits_by_guide: list[list] = [[] for _ in screen_pool]
    for ref_idx, (ref_name, ref_seq, ref_mask) in enumerate(references):
        is_locus = (ref_idx == 0)
        fast = len(ref_seq) > fast_index_threshold
        if fast:
            print(f"  [{ref_name}] {len(ref_seq):,} bp > {fast_index_threshold:,} bp "
                  f"threshold -> building a seed index once, reused for all "
                  f"{len(screen_pool)} candidates (seed-exact matches only at this scale).")
            index = build_index(ref_seq, ref_name=ref_name, seed_len=seed_len, mask=ref_mask)
        else:
            print(f"  [{ref_name}] {len(ref_seq):,} bp -> exhaustive scan")
        for gi, g in enumerate(screen_pool):
            on_pos = g.guide_start_fwd if is_locus else None
            on_strand = g.strand if is_locus else None
            if fast:
                hits_by_guide[gi].extend(scan_with_index(
                    g.protospacer, index, max_mismatches=max_mismatches,
                    on_target_position=on_pos, on_target_strand=on_strand))
            else:
                hits_by_guide[gi].extend(scan_reference(
                    g.protospacer, ref_seq, ref_name=ref_name,
                    max_mismatches=max_mismatches, seed_len=seed_len,
                    on_target_position=on_pos, on_target_strand=on_strand,
                    mask=ref_mask))
    print()

    screened = [(g, summarize_offtargets(hits_by_guide[gi]), hits_by_guide[gi])
                for gi, g in enumerate(screen_pool)]

    if not allow_exact_repeats:
        clean = [t for t in screened if t[1]["exact_offtargets"] == 0]
        excluded = len(screened) - len(clean)
        if clean:
            if excluded:
                print(f"Excluded {excluded} candidate(s) with an exact-match off-target "
                      f"elsewhere in the reference set (pass --allow-exact-repeats to keep them).")
            screened = clean
        else:
            print("WARNING: every screened candidate has an exact-match off-target "
                  "somewhere in the reference set. Keeping all, ranked by specificity, "
                  "so you can see the least-bad option — none are safe to order as-is.")

    if using_rs3:
        screened.sort(key=lambda t: (t[0].rs3_score is not None, t[0].rs3_score,
                                      t[1]["cfd_specificity_score"]), reverse=True)
    else:
        screened.sort(key=lambda t: (t[0].score, t[1]["cfd_specificity_score"]), reverse=True)

    rank_col = f"{'RS3':>7}" if using_rs3 else f"{'score':>6}"
    cds_col = f" {'CDS%':>5}" if (knockout and have_cds) else ""
    print(f"{'#':>2}  {'protospacer':<22} {'PAM':<4} {'str':<4} {'cut':>5} "
          f"{'GC%':>5} {rank_col}{cds_col} {'CFD-spec':>8} {'exact-OT':>9} {'rpt-OT':>7}")
    for i, (g, summ, _) in enumerate(screened[:top_n], 1):
        rank_val = f"{g.rs3_score:>7.2f}" if using_rs3 else f"{g.score:>6.0f}"
        cds_val = ""
        if knockout and have_cds:
            f = cds_frac(g)
            cds_val = f" {f * 100:>5.1f}" if f is not None else f" {'-':>5}"
        print(f"{i:>2}  {g.protospacer:<22} {g.pam:<4} {g.strand:<4} {g.cut_site_fwd:>5} "
              f"{g.gc:>5.0f} {rank_val}{cds_val} {summ['cfd_specificity_score']:>8.1f} "
              f"{summ['exact_offtargets']:>9} {summ['repeat_masked_hits']:>7}")

    # --- STEP 3: pick the top guide ---
    best, best_summary, best_hits = screened[0]
    print()
    print("-" * 70)
    print("STEP 3  Select guide + confirm Cas9 geometry")
    print("-" * 70)
    pam_pos = (best.guide_end_fwd if best.strand == "+"
               else best.guide_start_fwd - 3)
    print(f"Chosen guide  : {best.protospacer}")
    print(f"PAM           : {best.pam}")
    print(f"Strand        : {best.strand}  "
          f"({'sense' if best.strand=='+' else 'antisense'} — Cas9 uses either)")
    print(f"Guide span    : {best.guide_start_fwd}-{best.guide_end_fwd} (fwd coords)")
    print(f"Cut site      : {best.cut_site_fwd} (blunt, 3 bp 5' of PAM)")
    print(f"GC content    : {best.gc:.0f}%")
    if using_rs3:
        print(f"On-target score       : {best.rs3_score:.4f}  "
              f"(REAL Rule Set 3, Doench lab 2022 — {rule_set_3.mode()}, not reimplemented)")
        print(f"  (heuristic on-target score, kept for comparison: {best.score:.0f})")
    else:
        print(f"On-target score       : {best.score:.0f}  (transparent heuristic, not Doench 2016 — "
              f"see rule_set_3.py to enable the real Rule Set 3 score)")
    print(f"CFD specificity score : {best_summary['cfd_specificity_score']}  "
          f"(real Doench et al. 2016 CFD per-hit, aggregated — see cfd_score.py)")
    print(f"  (legacy heuristic specificity score: {best_summary['specificity_score']})")
    print(f"Off-targets found     : {best_summary['total_hits']} total, "
          f"{best_summary['exact_offtargets']} exact, "
          f"{best_summary['seed_intact_near_offtargets']} near-exact w/ intact seed, "
          f"{best_summary['repeat_masked_hits']} in repeat-masked sequence")
    if best_summary["worst_hit"]:
        w = best_summary["worst_hit"]
        print(f"  worst hit: {w.matched_seq} PAM={w.pam} in {w.ref_name!r} "
              f"@{w.position} ({w.strand} strand, {w.mismatches} mismatch(es), "
              f"{w.seed_mismatches} in seed{', repeat-masked' if w.repeat_masked else ''}) "
              f"— inspect before ordering.")
    if best_summary["worst_cfd_hit"] and best_summary["worst_cfd_hit"] is not best_summary["worst_hit"]:
        wc = best_summary["worst_cfd_hit"]
        print(f"  highest-CFD-risk hit: {wc.matched_seq} PAM={wc.pam} in {wc.ref_name!r} "
              f"@{wc.position} (CFD={wc.cfd:.3f})")
    if locus_mask and any(locus_mask[best.guide_start_fwd:best.guide_end_fwd]):
        print("  NOTE: the chosen guide's own target site overlaps repeat-masked "
              "(lowercase/RepeatMasker-flagged) sequence in the input locus — a real "
              "red flag for a repetitive, likely multi-copy site. Consider another "
              "guide, or confirm with a genome-scale off-target tool before ordering.")
    if knockout and have_cds:
        f = cds_frac(best)
        _total = cds_total_len if cds_total_len else (cds_end - cds_start)
        _coding_pos = cds_upstream_bp + (best.cut_site_fwd - cds_start)
        print(f"Position in CDS      : {f * 100:.1f}% "
              f"(coding base {_coding_pos} of {_total})")
        if f < 0.05:
            print("  WARNING: very close to the start codon. A frameshift here can be "
                  "rescued by translation reinitiation at a downstream ATG, so the "
                  "protein may survive. Prefer a cut slightly further in.")
        elif f > 0.60:
            print("  WARNING: late in the coding sequence. A frameshift here still "
                  "leaves most of the protein intact, so it may retain function. "
                  "Prefer an earlier cut.")
        else:
            print("  Good position: early enough to truncate the protein, far enough "
                  "in to avoid reinitiation rescue.")
    print()
    print("  sgRNA selects the site by hybridisation; RNP makes the DSB here.")

    # --- KNOCK-OUT: no donor. Report and verify, then stop. ---
    if knockout:
        print()
        print("-" * 70)
        print("STEP 4  Repair outcome (no donor required)")
        print("-" * 70)
        print("A knock-out relies on NHEJ, which is error-prone by design: repairing "
              "the blunt DSB introduces small insertions or deletions.")
        print("Roughly two thirds of indels shift the reading frame, truncating the "
              "protein downstream of the cut. The remaining in-frame indels may leave "
              "a partly functional protein, which is why knock-outs are screened.")
        print()
        print("What to order : the sgRNA only.")
        print("How to verify : amplify ~200-400 bp around the cut and sequence it "
              "(TIDE, ICE, or amplicon NGS) to measure the indel spectrum.")

        print()
        print("-" * 70)
        print("STEP 5  Verification")
        print("-" * 70)
        checks = [
            ("Guide found with a canonical NGG PAM", best.pam.endswith("GG")),
            ("Cut is 3 bp 5' of the PAM (SpCas9 geometry)",
             abs((best.guide_end_fwd if best.strand == "+"
                  else best.guide_start_fwd) - best.cut_site_fwd) == 3),
            (f"No exact off-target in screened reference set "
             f"({'ALLOWED via --allow-exact-repeats' if allow_exact_repeats else 'gate active'})",
             allow_exact_repeats or best_summary["exact_offtargets"] == 0),
        ]
        if have_cds:
            f = cds_frac(best)
            checks.append(("Cut lands inside the declared CDS", f is not None))
            checks.append(("Cut is in the frameshift-effective 5-60% window",
                            f is not None and 0.05 <= f <= 0.60))
        all_ok = True
        for name, ok in checks:
            print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
            all_ok &= ok
        print()
        print("RESULT:", "ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED")
        if not have_cds:
            print("NOTE: CDS-position checks were skipped because no "
                  "--cds-start/--cds-end was given.")
        return None

    # --- STEP 4: build the knock-in donor ---
    print()
    print("-" * 70)
    print("STEP 4  Build HDR knock-in donor")
    print("-" * 70)
    donor = design_knockin_donor(
        seq,
        cut_site=best.cut_site_fwd,
        insert_seq=cargo,
        insert_name=cargo_name,
        protospacer=best.protospacer,
        pam_pos=best.guide_end_fwd if best.strand == "+" else None,
        arm_len=arm_len,
    )
    print(donor.annotated_map())

    # --- STEP 5: the map ---
    print()
    print("-" * 70)
    print("STEP 5  Annotated donor map")
    print("-" * 70)
    print(donor.ascii_track())

    # --- STEP 6: verification checks ---
    print()
    print("-" * 70)
    print("STEP 6  Verification")
    print("-" * 70)
    la = [f for f in donor.features if f.name == "HA-L"][0]
    ra = [f for f in donor.features if f.name == "HA-R"][0]
    ins = [f for f in donor.features if f.kind == "insert"][0]

    checks = []
    checks.append(("Left arm matches genome upstream of cut",
                   donor.full_donor[la.start:la.end]
                   == seq[max(best.cut_site_fwd - arm_len, 0):best.cut_site_fwd]))
    checks.append(("Right arm matches genome downstream of cut",
                   donor.full_donor[ra.start:ra.end]
                   == seq[best.cut_site_fwd:best.cut_site_fwd + arm_len]))
    checks.append(("Cargo present and intact",
                   donor.full_donor[ins.start:ins.end] == cargo.upper()))
    checks.append(("Donor length = armL + cargo + armR",
                   len(donor.full_donor)
                   == (la.end - la.start) + len(cargo) + (ra.end - ra.start)))
    checks.append(("Cargo separates protospacer from PAM (no re-cut)",
                   cargo.upper() in donor.full_donor
                   and best.protospacer + best.pam not in donor.full_donor))
    checks.append((f"No exact off-target in screened reference set "
                   f"({'ALLOWED via --allow-exact-repeats' if allow_exact_repeats else 'gate active'})",
                   allow_exact_repeats or best_summary["exact_offtargets"] == 0))

    all_ok = True
    for name, ok in checks:
        print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
        all_ok &= ok
    print()
    print("RESULT:", "ALL CHECKS PASSED" if all_ok else "SOME CHECKS FAILED")
    return donor


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="End-to-end CRISPR pipeline: find guides, off-target screen "
                     "them, then either build an HDR knock-in donor or design a "
                     "frameshift knock-out.")
    p.add_argument("fasta", nargs="?",
                    help="Target locus FASTA (single record). Omit to run the built-in demo.")
    p.add_argument("--mode", choices=["knockin", "knockout"], default="knockin",
                    help="knockin (default): insert cargo by HDR, build a donor. "
                         "knockout: disrupt the gene by NHEJ frameshift, no donor.")
    p.add_argument("--cds-start", type=int, default=None, metavar="N",
                    help="0-based offset of the coding sequence START within your "
                         "input FASTA. Knock-out mode uses this to rank guides by "
                         "position in the CDS and to reject cuts outside it.")
    p.add_argument("--cds-end", type=int, default=None, metavar="N",
                    help="0-based, end-exclusive offset of the coding sequence END "
                         "within your input FASTA.")
    p.add_argument("--cds-upstream-bp", type=int, default=0, metavar="N",
                    help="Coding bases that lie BEFORE this window in the full CDS. "
                         "Use with --cds-total-len when your FASTA holds one exon of "
                         "a multi-exon gene. Default 0.")
    p.add_argument("--cds-total-len", type=int, default=None, metavar="N",
                    help="Length of the COMPLETE coding sequence across all exons. "
                         "Without it, position-in-CDS is measured against this window "
                         "alone, which overstates how far in the cut is.")
    p.add_argument("--reference", "-r", nargs="+", default=[], metavar="FASTA",
                    help="FASTA file(s) to off-target-screen candidates against "
                         "(paralogs, pseudogenes, a chromosome arm, ...). Each may "
                         "hold multiple records. Always screened against the input "
                         "locus itself in addition to these.")
    p.add_argument("--cargo-fasta", metavar="FASTA",
                    help="FASTA with the real cargo/insert sequence. Overrides the "
                         "built-in CAR-like placeholder — use this for real designs.")
    p.add_argument("--cargo-name", default="CAR")
    p.add_argument("--arm-len", type=int, default=40, help="Homology arm length (bp).")
    p.add_argument("--top-n", type=int, default=5, help="Rows to display in ranking tables.")
    p.add_argument("--screen-n", type=int, default=15,
                    help="How many top on-target candidates to off-target-screen.")
    p.add_argument("--max-mismatches", type=int, default=6,
                    help="Max mismatches to still count as an off-target hit. Default 6: "
                         "a real, literature-validated off-target (see validate_iguide.py) "
                         "has 5 mismatches and is invisible at 4.")
    p.add_argument("--seed-len", type=int, default=12,
                    help="PAM-proximal seed length (nt) used for specificity weighting.")
    p.add_argument("--allow-exact-repeats", action="store_true",
                    help="Don't auto-exclude candidates with an exact-match off-target "
                         "elsewhere in the screened reference set.")
    p.add_argument("--fast-index-threshold", type=int, default=2_000_000, metavar="BP",
                    help="References longer than this (bp) use a seed-indexed fast scan "
                         "(built once, reused per candidate) instead of an exhaustive "
                         "per-guide scan — needed to screen against real chromosome-scale "
                         "human references. Default 2,000,000.")
    return p


if __name__ == "__main__":
    args = build_argparser().parse_args()

    if args.cargo_fasta:
        cargo_header, CARGO = read_fasta(args.cargo_fasta)
    else:
        # A CAR-like cargo (placeholder cassette, not a real CAR ORF)
        CARGO = ("ATGCTGCTGCTGGTGACCAGCCTGCTGCTGTGCGAGCTGCCCCACCCCGCC"
                 "TTCCTGCTGATCCCCGACATCCAGATGACCCAGAGCCCCAGCAGCCTGAGC")

    # Reference panels use the mask-preserving reader: real human genome
    # downloads (UCSC/Ensembl) commonly soft-mask repeats as lowercase, and
    # that's exactly the sequence off-target screening most needs to see.
    reference_seqs = []
    for path in args.reference:
        reference_seqs.extend(
            (f"{path}:{h}" if h else path, s, m)
            for h, s, m in read_masked_fasta_records(path))

    run_kwargs = dict(
        cargo_name=args.cargo_name, arm_len=args.arm_len, top_n=args.top_n,
        reference_seqs=reference_seqs, screen_n=args.screen_n,
        max_mismatches=args.max_mismatches, seed_len=args.seed_len,
        allow_exact_repeats=args.allow_exact_repeats,
        fast_index_threshold=args.fast_index_threshold,
        mode=args.mode, cds_start=args.cds_start, cds_end=args.cds_end,
        cds_upstream_bp=args.cds_upstream_bp, cds_total_len=args.cds_total_len,
    )

    if args.fasta:
        records = read_masked_fasta_records(args.fasta)
        if len(records) != 1:
            raise SystemExit(
                f"{args.fasta} has {len(records)} records but run_pipeline expects a "
                f"single target locus. Use --reference for multi-sequence off-target panels.")
        header, seq, locus_mask = records[0]
        run(seq, CARGO, label=header or args.fasta, locus_mask=locus_mask, **run_kwargs)
    else:
        # Built-in demo locus: synthetic, AAVS1-style, with a verifiable guide.
        DEMO = (
            "GGGGCCACTAGGGACAGGATTGGTGACAGAAAAGCCCCATCCTTAGGCCTC"
            "CTCCTTCCTAGTCTCCTGATATTGGGTCTAACCCCCACCTCCTGTTAGGCA"
            "GATTCCTTATCTGGTGACACACCCCCATTTCCTGGAGCCATCTCTCTCCTT"
            "GCCAGAACCTCTAAGGTTTGCTTACGATGGAGCCAGAGAGGATCCTGGGAG"
        )
        run(DEMO, CARGO, label="synthetic AAVS1-style demo locus", **run_kwargs)
