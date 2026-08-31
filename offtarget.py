"""
Off-target scanning for SpCas9 guides.

WHAT THIS DOES
--------------
Given a 20 nt protospacer, scans one or more reference sequences for every
NGG/NAG PAM site on BOTH strands and reports how closely the adjacent 20mer
matches the query guide. This is a real mismatch search over whatever
sequence you supply — not a stub — but it is REFERENCE-SCOPED: you decide
what to search (the on-target locus itself, known paralogs/pseudogenes, a
chromosome arm, ...). There is no bundled whole-genome index, so this
cannot replace a genome-wide aligner.

TWO SEARCH MODES (for real, chromosome-scale human loci)
----------------------------------------------------------
- `scan_reference()` — exhaustive. O(len(reference)) per guide, finds every
  hit within max_mismatches ANYWHERE in the 20mer. Correct and thorough,
  but re-scans the whole reference for every guide, so it doesn't scale to
  repeatedly screening many candidates against a multi-megabase reference.
- `build_index()` + `scan_with_index()` — indexed. Builds a hash of the
  PAM-proximal seed (the ~8-12 nt most critical for binding, see below)
  ONCE per reference, then each guide is a near-O(1) dict lookup. This
  scales to real human chromosome-sized references, at one deliberate
  cost: it only finds hits with an EXACT seed match. Crucially, an exact
  20nt match always has an exact seed by construction, so the
  safety-critical "exact off-target elsewhere" gate in run_pipeline.py
  stays fully accurate in this mode — only the lower-risk, seed-mismatched
  partial matches are left un-enumerated at this scale.

`run_pipeline.py` picks between the two automatically based on reference
size (see --fast-index-threshold).

WHY THE DEFAULT IS 6 MISMATCHES, NOT 4
-----------------------------------------
`max_mismatches` defaults to 6. It used to default to 4, until
`validate_iguide.py` showed that was too strict to see a real thing: a
literature-validated off-target (PHACTR1, Nature 2022,
https://doi.org/10.1038/s41586-022-05140-y, Supplementary Table 8, found
by iGUIDE and confirmed by 50,000x-deep sequencing) has 5 mismatches
against its guide. At max_mismatches=4 this tool would report a clean
scan and miss a real, published, wet-lab-confirmed off-target site. Run
`python3 validate_iguide.py` to see the reproduction.

WHY THE SEED REGION MATTERS
----------------------------
Hsu et al. 2013 (Nat Biotechnol 31:827-832) showed the ~8-12 nt of the
protospacer closest to the PAM (the "seed") is far more sensitive to
mismatches than the PAM-distal end: a site with an intact seed and a few
distal mismatches can still be cut, while a single seed mismatch often
abolishes binding. NAG is also a documented (weaker) SpCas9-tolerated PAM,
so it's included here even though sgrna_design.py only *designs* against
canonical NGG sites.

REPEAT-MASKING (real human genome FASTA)
------------------------------------------
Real human genome downloads (UCSC, Ensembl, NCBI) commonly soft-mask
repetitive/low-complexity sequence as lowercase (the RepeatMasker
convention). A guide landing in repeat-masked sequence is a red flag on its
own: by definition that sequence occurs elsewhere in the genome too.
`read_masked_fasta_records()` preserves this signal instead of discarding
it, and hits/candidates overlapping masked sequence are flagged
`repeat_masked`. (Plain-FASTA sources like a raw NCBI efetch usually have
no soft-masking at all — this only helps when your reference has it.)

REAL CFD SCORE (Doench et al. 2016), NOT JUST A HEURISTIC
-------------------------------------------------------------
Every hit's `.cfd` field is the actual published Doench, Fusi et al. 2016
Cutting Frequency Determination score (see cfd_score.py) — a real
per-position mismatch penalty matrix + PAM score table, transcribed
programmatically from CRISPOR's public data files, not reconstructed from
memory. `summarize()` reports `cfd_specificity_score`, an aggregate built
by summing each hit's real CFD score and inverting — the per-hit number is
the published model; the "sum many real scores into one guide-level
number" step is still this module's own convention (Doench et al. don't
define a single canonical guide-level aggregate). `specificity_score()`
is kept as the original from-scratch heuristic for comparison, clearly
labeled as such below — prefer `cfd_specificity_score` when you need a
citable number.

For genome-scale, publication-grade off-target discovery use Cas-OFFinder,
CRISPOR, or FlashFry, and validate experimentally (GUIDE-seq, CIRCLE-seq,
etc) before trusting any guide for real therapeutic or research use. This
toolkit does not attempt on-target activity prediction (Rule Set 2/3,
DeepHF, ...) — those are trained ML models with too many fitted
parameters to responsibly reproduce from a paper alone; sgrna_design.py's
on-target score remains a documented heuristic.
"""

from __future__ import annotations
import re
from dataclasses import dataclass

from cfd_score import cfd_score

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


@dataclass
class OffTargetHit:
    ref_name: str
    strand: str          # "+" or "-", always in forward coords of the reference
    position: int        # 0-based start of the matched 20mer, forward coords
    matched_seq: str
    pam: str
    mismatches: int
    seed_mismatches: int
    is_self: bool = False
    repeat_masked: bool = False
    cfd: float | None = None   # real Doench et al. 2016 CFD score; None if not computable (see _cfd_for)


def _mismatches(a: str, b: str) -> int:
    return sum(x != y for x, y in zip(a, b))


def _cfd_for(query: str, candidate: str, pam: str) -> float | None:
    """cfd_score() only defines the model for 20nt protospacers; return
    None rather than raise for any off-spec guide length so callers using
    a customized sgrna_design guide_len don't crash."""
    if len(query) != 20 or len(candidate) != 20:
        return None
    return cfd_score(query, candidate, pam)


# ---------- FASTA reading that keeps repeat-masking ----------
def read_masked_fasta_records(path: str) -> list[tuple[str, str, list[bool]]]:
    """
    Parse a FASTA file, preserving lowercase soft-masking instead of
    discarding it. Returns (header, UPPERCASE sequence, is_repeat_masked)
    per record, where is_repeat_masked is a same-length list[bool].
    """
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
                parts.append(line)
    if header is not None:
        records.append((header, "".join(parts)))

    out = []
    for h, raw in records:
        bad = set(raw.upper()) - set("ACGTN")
        if bad:
            raise ValueError(f"Non-DNA characters in FASTA record {h!r}: {sorted(bad)}")
        out.append((h, raw.upper(), [c.islower() for c in raw]))
    return out


# ---------- Exhaustive scan (small/moderate references) ----------
def _scan_strand(query: str, strand_seq: str, ref_name: str, strand_label: str,
                  max_mismatches: int, seed_len: int,
                  mask: list[bool] | None = None) -> list[OffTargetHit]:
    guide_len = len(query)
    pam_regex = r"(?=([ACGT][AG]G))"   # NGG or NAG
    hits = []
    for m in re.finditer(pam_regex, strand_seq):
        pam_start = m.start(1)
        g_start = pam_start - guide_len
        if g_start < 0:
            continue
        candidate = strand_seq[g_start:pam_start]
        if len(candidate) != guide_len or set(candidate) - set("ACGT"):
            continue
        mm = _mismatches(query, candidate)
        if mm > max_mismatches:
            continue
        seed_mm = _mismatches(query[-seed_len:], candidate[-seed_len:])
        pam = strand_seq[pam_start:pam_start + 3]
        masked = bool(mask) and any(mask[g_start:pam_start])
        hits.append(OffTargetHit(ref_name, strand_label, g_start, candidate,
                                  pam, mm, seed_mm, repeat_masked=masked,
                                  cfd=_cfd_for(query, candidate, pam)))
    return hits


def scan_reference(protospacer: str, reference: str, ref_name: str = "reference",
                    max_mismatches: int = 6, seed_len: int = 12,
                    on_target_position: int | None = None,
                    on_target_strand: str | None = None,
                    mask: list[bool] | None = None) -> list[OffTargetHit]:
    """
    Exhaustively find every NGG/NAG PAM site on both strands of `reference`
    whose adjacent 20mer is within `max_mismatches` of `protospacer`.
    O(len(reference)) — fine for a locus or a multi-kb/low-Mb panel, not
    meant to be called per-guide against a whole chromosome (use
    build_index()/scan_with_index() for that).

    If on_target_position/on_target_strand are given (forward coords, as
    returned by sgrna_design.find_guides), the matching hit is flagged
    is_self=True so summaries don't count the intended edit site as a risk.
    `mask`, if given, is a forward-coordinate list[bool] from
    read_masked_fasta_records() marking repeat-masked bases.
    """
    protospacer = protospacer.upper()
    reference = reference.upper()

    fwd_hits = _scan_strand(protospacer, reference, ref_name, "+",
                             max_mismatches, seed_len, mask)

    rc = revcomp(reference)
    rc_mask = mask[::-1] if mask is not None else None
    L = len(reference)
    rev_hits = []
    for h in _scan_strand(protospacer, rc, ref_name, "-", max_mismatches, seed_len, rc_mask):
        fwd_pos = L - (h.position + len(protospacer))
        rev_hits.append(OffTargetHit(ref_name, "-", fwd_pos, h.matched_seq,
                                      h.pam, h.mismatches, h.seed_mismatches,
                                      repeat_masked=h.repeat_masked, cfd=h.cfd))

    hits = fwd_hits + rev_hits
    if on_target_position is not None and on_target_strand is not None:
        for h in hits:
            if h.position == on_target_position and h.strand == on_target_strand:
                h.is_self = True

    hits.sort(key=lambda h: (h.mismatches, h.seed_mismatches))
    return hits


# ---------- Indexed scan (real chromosome-scale references) ----------
@dataclass
class ReferenceIndex:
    name: str
    length: int
    seed_len: int
    guide_len: int
    fwd_seq: str
    rc_seq: str
    fwd_mask: list | None
    rc_mask: list | None
    seed_map: dict   # seed(str) -> list[(local_pos, strand)]


def build_index(reference: str, ref_name: str = "reference", seed_len: int = 12,
                 guide_len: int = 20, mask: list[bool] | None = None) -> ReferenceIndex:
    """
    Build once, scan many: index every NGG/NAG PAM site on both strands of
    `reference` by its PAM-proximal seed. O(len(reference)) build cost paid
    ONCE, then scan_with_index() is a near-O(1) lookup per guide — this is
    what makes screening dozens of candidates against a real multi-megabase
    human reference (a chromosome arm, a paralog panel) tractable in pure
    Python.
    """
    reference = reference.upper()
    rc = revcomp(reference)
    rc_mask = mask[::-1] if mask is not None else None

    seed_map: dict[str, list[tuple[int, str]]] = {}
    pam_regex = r"(?=([ACGT][AG]G))"

    def index_strand(strand_seq: str, strand_label: str):
        for m in re.finditer(pam_regex, strand_seq):
            pam_start = m.start(1)
            g_start = pam_start - guide_len
            if g_start < 0:
                continue
            candidate = strand_seq[g_start:pam_start]
            if len(candidate) != guide_len or set(candidate) - set("ACGT"):
                continue
            seed = candidate[-seed_len:]
            seed_map.setdefault(seed, []).append((g_start, strand_label))

    index_strand(reference, "+")
    index_strand(rc, "-")

    return ReferenceIndex(ref_name, len(reference), seed_len, guide_len,
                           reference, rc, mask, rc_mask, seed_map)


def scan_with_index(protospacer: str, index: ReferenceIndex, max_mismatches: int = 6,
                     on_target_position: int | None = None,
                     on_target_strand: str | None = None) -> list[OffTargetHit]:
    """
    Fast path built on a pre-built ReferenceIndex. Only finds hits whose
    PAM-proximal seed matches the query EXACTLY. An exact 20nt match always
    has an exact seed, so the exact-off-target safety gate is unaffected —
    only lower-risk, seed-mismatched partial matches go undetected at this
    scale. See module docstring.
    """
    protospacer = protospacer.upper()
    seed = protospacer[-index.seed_len:]
    L = index.length
    hits = []
    for g_start, strand in index.seed_map.get(seed, []):
        if strand == "+":
            strand_seq, strand_mask = index.fwd_seq, index.fwd_mask
        else:
            strand_seq, strand_mask = index.rc_seq, index.rc_mask
        pam_start = g_start + index.guide_len
        candidate = strand_seq[g_start:pam_start]
        mm = _mismatches(protospacer, candidate)
        if mm > max_mismatches:
            continue
        seed_mm = _mismatches(protospacer[-index.seed_len:], candidate[-index.seed_len:])
        pam = strand_seq[pam_start:pam_start + 3]
        masked = bool(strand_mask) and any(strand_mask[g_start:pam_start])
        fwd_pos = g_start if strand == "+" else L - pam_start

        h = OffTargetHit(index.name, strand, fwd_pos, candidate, pam, mm, seed_mm,
                          repeat_masked=masked, cfd=_cfd_for(protospacer, candidate, pam))
        if (on_target_position is not None and fwd_pos == on_target_position
                and strand == on_target_strand):
            h.is_self = True
        hits.append(h)

    hits.sort(key=lambda h: (h.mismatches, h.seed_mismatches))
    return hits


def specificity_score(hits: list[OffTargetHit]) -> float:
    """
    LEGACY from-scratch heuristic, kept for comparison (see module
    docstring for the real, published alternative: cfd_specificity_score
    in summarize()). Each non-self hit is weighted by mismatch count and
    whether its seed is intact — seed-intact hits (Hsu et al. 2013) are
    the ones most likely to actually get cut, so they dominate the
    penalty. Returns 0-100, higher = more specific.
    """
    burden = 0.0
    for h in hits:
        if h.is_self:
            continue
        seed_intact = h.seed_mismatches == 0
        weight = (1.0 / (1 + h.mismatches)) * (3.0 if seed_intact else 1.0)
        burden += weight
    return round(100.0 / (1.0 + burden), 1)


def cfd_specificity_score(hits: list[OffTargetHit]) -> float:
    """
    Aggregate built from the REAL Doench et al. 2016 CFD score (see
    cfd_score.py): sum each non-self hit's published per-site CFD score,
    then invert. The per-hit input is the real model; summing many real
    scores into one guide-level number is this module's own convention
    (Doench et al. report per-site scores, not a single aggregate).
    Hits where CFD isn't defined (non-20nt guide) are skipped. Returns
    0-100, higher = more specific.
    """
    total = sum(h.cfd for h in hits if not h.is_self and h.cfd is not None)
    return round(100.0 / (1.0 + total), 1)


def summarize(hits: list[OffTargetHit]) -> dict:
    """Roll a hit list up into the numbers a user actually wants to see."""
    real_hits = [h for h in hits if not h.is_self]
    exact = [h for h in real_hits if h.mismatches == 0]
    seed_intact_near = [h for h in real_hits
                         if 0 < h.mismatches <= 3 and h.seed_mismatches == 0]
    repeat_hits = [h for h in real_hits if h.repeat_masked]
    cfd_hits = [h for h in real_hits if h.cfd is not None]
    worst_cfd_hit = max(cfd_hits, key=lambda h: h.cfd, default=None)
    return {
        "total_hits": len(real_hits),
        "exact_offtargets": len(exact),
        "seed_intact_near_offtargets": len(seed_intact_near),
        "repeat_masked_hits": len(repeat_hits),
        "specificity_score": specificity_score(hits),          # legacy heuristic
        "cfd_specificity_score": cfd_specificity_score(hits),  # real Doench et al. 2016 CFD, aggregated
        "worst_hit": real_hits[0] if real_hits else None,
        "worst_cfd_hit": worst_cfd_hit,
    }


# ---------------- Self-test ----------------
def _self_test():
    query = "ACGTACGTACGTACGTACGT"  # 20 nt
    onshot = query + "TGG"                                  # exact, PAM=TGG
    exact_offtarget = "TTTT" + query + "AGG" + "TTTT"        # exact repeat elsewhere
    near_seed_mm = "GGGG" + "ACGTACGTACGT" + "TCGTACGT" + "CGG" + "GGGG"
    reference = onshot + "NNNNNNNNNN" + exact_offtarget + "NNNNNNNNNN" + near_seed_mm

    on_pos = 0
    hits = scan_reference(query, reference, ref_name="synthetic",
                           max_mismatches=4, seed_len=12,
                           on_target_position=on_pos, on_target_strand="+")

    self_hits = [h for h in hits if h.is_self]
    assert len(self_hits) == 1, f"expected exactly 1 self hit, got {len(self_hits)}"

    summary = summarize(hits)
    assert summary["exact_offtargets"] >= 1, "failed to find the planted exact off-target"
    assert summary["total_hits"] >= 2, f"expected >=2 non-self hits, got {summary}"
    assert 0 <= summary["specificity_score"] <= 100

    clean_ref = query + "TGG" + "N" * 200
    clean_hits = scan_reference(query, clean_ref, ref_name="clean",
                                 max_mismatches=4, seed_len=12,
                                 on_target_position=0, on_target_strand="+")
    clean_summary = summarize(clean_hits)
    assert clean_summary["specificity_score"] == 100.0, clean_summary
    assert clean_summary["cfd_specificity_score"] == 100.0, clean_summary

    # --- real CFD score is populated and the exact off-target (0 mismatches)
    # must have cfd == 1.0 (no predicted activity loss) ---
    exact_hit = next(h for h in hits if not h.is_self and h.mismatches == 0)
    assert exact_hit.cfd == 1.0, exact_hit
    assert 0 <= summary["cfd_specificity_score"] <= 100

    # --- indexed fast path must agree with the exhaustive scan on
    # seed-exact hits (the exact off-target above has 0 seed mismatches,
    # so it MUST still be found via the index; the near_seed_mm hit has a
    # seed mismatch, so it's expected to be ABSENT from the indexed result) ---
    index = build_index(reference, ref_name="synthetic", seed_len=12)
    indexed_hits = scan_with_index(query, index, max_mismatches=4,
                                    on_target_position=on_pos, on_target_strand="+")
    indexed_summary = summarize(indexed_hits)
    assert indexed_summary["exact_offtargets"] == summary["exact_offtargets"], (
        "indexed scan must find the same exact off-targets as the exhaustive scan",
        indexed_summary, summary)
    assert indexed_summary["total_hits"] <= summary["total_hits"], (
        "indexed scan (seed-exact only) should never find MORE hits than exhaustive")

    # --- repeat-masking: a lowercase-masked exact off-target should be flagged ---
    masked_ref = onshot + "n" * 10 + exact_offtarget.lower() + "n" * 10
    mask = [c.islower() for c in masked_ref]
    masked_hits = scan_reference(query, masked_ref, ref_name="masked",
                                  max_mismatches=4, seed_len=12,
                                  on_target_position=0, on_target_strand="+",
                                  mask=mask)
    masked_summary = summarize(masked_hits)
    assert masked_summary["repeat_masked_hits"] >= 1, masked_summary

    print("Off-target self-test passed.")
    print(f"  exhaustive summary : {summary}")
    print(f"  indexed summary    : {indexed_summary}")
    print(f"  clean-guide summary: {clean_summary}")
    print(f"  repeat-mask summary: {masked_summary}")


if __name__ == "__main__":
    _self_test()
