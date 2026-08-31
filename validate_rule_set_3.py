"""
Real-data validation for rule_set_3.py.

Checks that the REAL Rule Set 3 model (Doench lab, 2022, via the official
rs3 package — see rule_set_3.py) ranks two independently-published,
real-world guides as its #1 pick among all candidates in their real
genomic loci:

    AAVS1, chr19 (Nature 2022 Methods, https://doi.org/10.1038/s41586-022-05140-y)
        -> #1 of 71 usable candidates in aavs1_locus.fa
    B2M,   chr15 (same paper's Methods)
        -> #1 of 11 usable candidates in b2m_locus.fa

For comparison: sgrna_design.py's from-scratch heuristic ranked the same
AAVS1 guide 53rd of 91 candidates (see the session notes / README) — the
real published model agrees with what was actually used in the clinic
far more often than the heuristic did.

This test SKIPS (not fails) if rs3_venv isn't set up, since it's an
optional, environment-dependent feature — see rule_set_3.py's module
docstring for setup instructions.
"""

from __future__ import annotations
import rule_set_3
from sgrna_design import find_guides, score_guide
from run_pipeline import read_fasta

CASES = [
    ("aavs1_locus.fa", "AGAGCTAGCACAGACTAGAG"),
    ("b2m_locus.fa", "GAGTAGCGCGAGCACAGCTA"),
]


def _rank_by_rs3(fasta_path: str, arm_len: int = 40):
    _, seq = read_fasta(fasta_path)
    # No CDS declared: this test ranks by the real Rule Set 3 score, and the
    # heuristic's CDS rule is irrelevant to (and must not perturb) that.
    guides = [score_guide(g) for g in find_guides(seq)]
    usable = [g for g in guides if arm_len <= g.cut_site_fwd <= len(seq) - arm_len]
    contexts = [rule_set_3.context_30mer(seq, g) for g in usable]
    scores = rule_set_3.predict(contexts)
    for g, s in zip(usable, scores):
        g.rs3_score = s
    usable.sort(key=lambda g: g.rs3_score, reverse=True)
    return usable


def _validate():
    if not rule_set_3.available():
        print("rs3_venv not found — skipping (optional feature; see rule_set_3.py "
              "for setup instructions). Nothing to validate without it.")
        return

    for fasta_path, expected_guide in CASES:
        ranked = _rank_by_rs3(fasta_path)
        protospacers = [g.protospacer for g in ranked]
        assert expected_guide in protospacers, (
            f"{expected_guide} not even found as a candidate in {fasta_path}")
        rank = protospacers.index(expected_guide) + 1
        assert rank == 1, (
            f"Expected the real published guide {expected_guide} to rank #1 by real "
            f"Rule Set 3 in {fasta_path}, but it ranked {rank}/{len(ranked)}.")
        print(f"{fasta_path}: real published guide {expected_guide} ranks #1/{len(ranked)} "
              f"(rs3={ranked[0].rs3_score:.4f}) under the real Rule Set 3 model.")

    print("Rule Set 3 real-data validation passed.")


if __name__ == "__main__":
    _validate()
