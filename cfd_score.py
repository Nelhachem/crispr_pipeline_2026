"""
Real CFD (Cutting Frequency Determination) off-target score.

Doench, Fusi et al. 2016, "Optimized sgRNA design to maximize activity
and minimize off-target effects of CRISPR-Cas9", Nat Biotechnol 34:184-191
(https://doi.org/10.1038/nbt.3437).

Unlike offtarget.py's specificity_score() (a documented from-scratch
heuristic), this IS the published model: a per-position mismatch penalty
matrix (240 entries: 12 possible RNA:DNA mismatch identities x 20
positions) and a PAM score table (16 entries, keyed by the PAM's last 2
nt), both fitted from the paper's experimental cutting-frequency data. A
candidate site's CFD score is the product of the relevant per-position
mismatch penalties times the PAM score: 1.0 = predicted fully active
(no penalty), 0.0 = predicted inactive.

DATA SOURCE
-----------
MISMATCH_SCORES and PAM_SCORES below are transcribed programmatically
(not retyped by hand) from mismatch_score.pkl / pam_scores.pkl in the
CRISPOR project's public repository:
    https://github.com/maximilianh/crisporWebsite/tree/master/CFD_Scoring
CRISPOR (Haeussler et al. 2016) is a widely-used, actively-maintained
CRISPR guide design tool whose CFD implementation is the standard
reference reproduction of Doench et al.'s original supplementary data.
The algorithm below (cfd_score()) is a faithful line-for-line port of
that repository's cfd-score-calculator.py, verified against it: see
validate_iguide.py and this module's self-test for reproductions of real
data points computed independently from the original pickle files.

SCOPE
-----
- Scores single-nucleotide-substitution mismatches only (no bulges/
  indels) -- exactly what offtarget.py's mismatch scanner finds.
- Defined for a 20nt SpCas9 protospacer + PAM; the score only depends on
  the PAM's final 2 nt (position 1, the "N" of NGG, is unconstrained).
- This scores a *specific candidate site* against a *specific guide* --
  it is an off-target risk score, not an on-target activity predictor.
  On-target predictors (Rule Set 1/2/3, DeepHF, CRISPRater, ...) are
  trained gradient-boosted-tree / neural-network models with hundreds to
  thousands of fitted parameters; they cannot be responsibly reproduced
  from published papers alone the way this lookup-table model can, so
  this toolkit does not attempt to reimplement them from scratch.
"""

from __future__ import annotations

_COMPLEMENT = str.maketrans("ACGTUacgtu", "TGCAAtgcaa")


def _complement(base: str) -> str:
    return base.translate(_COMPLEMENT)


MISMATCH_SCORES = {
    'rA:dA,1': 1.0,
    'rA:dA,10': 0.882352941,
    'rA:dA,11': 0.307692308,
    'rA:dA,12': 0.333333333,
    'rA:dA,13': 0.3,
    'rA:dA,14': 0.533333333,
    'rA:dA,15': 0.2,
    'rA:dA,16': 0.0,
    'rA:dA,17': 0.133333333,
    'rA:dA,18': 0.5,
    'rA:dA,19': 0.538461538,
    'rA:dA,2': 0.727272727,
    'rA:dA,20': 0.6,
    'rA:dA,3': 0.705882353,
    'rA:dA,4': 0.636363636,
    'rA:dA,5': 0.363636364,
    'rA:dA,6': 0.7142857140000001,
    'rA:dA,7': 0.4375,
    'rA:dA,8': 0.428571429,
    'rA:dA,9': 0.6,
    'rA:dC,1': 1.0,
    'rA:dC,10': 0.5555555560000001,
    'rA:dC,11': 0.65,
    'rA:dC,12': 0.7222222220000001,
    'rA:dC,13': 0.6521739129999999,
    'rA:dC,14': 0.46666666700000003,
    'rA:dC,15': 0.65,
    'rA:dC,16': 0.192307692,
    'rA:dC,17': 0.176470588,
    'rA:dC,18': 0.4,
    'rA:dC,19': 0.375,
    'rA:dC,2': 0.8,
    'rA:dC,20': 0.764705882,
    'rA:dC,3': 0.611111111,
    'rA:dC,4': 0.625,
    'rA:dC,5': 0.72,
    'rA:dC,6': 0.7142857140000001,
    'rA:dC,7': 0.705882353,
    'rA:dC,8': 0.7333333329999999,
    'rA:dC,9': 0.666666667,
    'rA:dG,1': 0.857142857,
    'rA:dG,10': 0.333333333,
    'rA:dG,11': 0.4,
    'rA:dG,12': 0.263157895,
    'rA:dG,13': 0.21052631600000002,
    'rA:dG,14': 0.214285714,
    'rA:dG,15': 0.272727273,
    'rA:dG,16': 0.0,
    'rA:dG,17': 0.176470588,
    'rA:dG,18': 0.19047619,
    'rA:dG,19': 0.20689655199999998,
    'rA:dG,2': 0.7857142859999999,
    'rA:dG,20': 0.22727272699999998,
    'rA:dG,3': 0.428571429,
    'rA:dG,4': 0.352941176,
    'rA:dG,5': 0.5,
    'rA:dG,6': 0.454545455,
    'rA:dG,7': 0.4375,
    'rA:dG,8': 0.428571429,
    'rA:dG,9': 0.571428571,
    'rC:dA,1': 1.0,
    'rC:dA,10': 0.9411764709999999,
    'rC:dA,11': 0.307692308,
    'rC:dA,12': 0.538461538,
    'rC:dA,13': 0.7,
    'rC:dA,14': 0.7333333329999999,
    'rC:dA,15': 0.066666667,
    'rC:dA,16': 0.307692308,
    'rC:dA,17': 0.46666666700000003,
    'rC:dA,18': 0.642857143,
    'rC:dA,19': 0.46153846200000004,
    'rC:dA,2': 0.9090909090000001,
    'rC:dA,20': 0.3,
    'rC:dA,3': 0.6875,
    'rC:dA,4': 0.8,
    'rC:dA,5': 0.636363636,
    'rC:dA,6': 0.9285714290000001,
    'rC:dA,7': 0.8125,
    'rC:dA,8': 0.875,
    'rC:dA,9': 0.875,
    'rC:dC,1': 0.913043478,
    'rC:dC,10': 0.38888888899999996,
    'rC:dC,11': 0.25,
    'rC:dC,12': 0.444444444,
    'rC:dC,13': 0.13636363599999998,
    'rC:dC,14': 0.0,
    'rC:dC,15': 0.05,
    'rC:dC,16': 0.153846154,
    'rC:dC,17': 0.058823529000000006,
    'rC:dC,18': 0.133333333,
    'rC:dC,19': 0.125,
    'rC:dC,2': 0.695652174,
    'rC:dC,20': 0.058823529000000006,
    'rC:dC,3': 0.5,
    'rC:dC,4': 0.5,
    'rC:dC,5': 0.6,
    'rC:dC,6': 0.5,
    'rC:dC,7': 0.470588235,
    'rC:dC,8': 0.642857143,
    'rC:dC,9': 0.6190476189999999,
    'rC:dT,1': 1.0,
    'rC:dT,10': 0.8666666670000001,
    'rC:dT,11': 0.75,
    'rC:dT,12': 0.7142857140000001,
    'rC:dT,13': 0.384615385,
    'rC:dT,14': 0.35,
    'rC:dT,15': 0.222222222,
    'rC:dT,16': 1.0,
    'rC:dT,17': 0.46666666700000003,
    'rC:dT,18': 0.538461538,
    'rC:dT,19': 0.428571429,
    'rC:dT,2': 0.727272727,
    'rC:dT,20': 0.5,
    'rC:dT,3': 0.8666666670000001,
    'rC:dT,4': 0.842105263,
    'rC:dT,5': 0.571428571,
    'rC:dT,6': 0.9285714290000001,
    'rC:dT,7': 0.75,
    'rC:dT,8': 0.65,
    'rC:dT,9': 0.857142857,
    'rG:dA,1': 1.0,
    'rG:dA,10': 0.8125,
    'rG:dA,11': 0.384615385,
    'rG:dA,12': 0.384615385,
    'rG:dA,13': 0.3,
    'rG:dA,14': 0.26666666699999997,
    'rG:dA,15': 0.14285714300000002,
    'rG:dA,16': 0.0,
    'rG:dA,17': 0.25,
    'rG:dA,18': 0.666666667,
    'rG:dA,19': 0.666666667,
    'rG:dA,2': 0.636363636,
    'rG:dA,20': 0.7,
    'rG:dA,3': 0.5,
    'rG:dA,4': 0.363636364,
    'rG:dA,5': 0.3,
    'rG:dA,6': 0.666666667,
    'rG:dA,7': 0.571428571,
    'rG:dA,8': 0.625,
    'rG:dA,9': 0.533333333,
    'rG:dG,1': 0.7142857140000001,
    'rG:dG,10': 0.4,
    'rG:dG,11': 0.428571429,
    'rG:dG,12': 0.529411765,
    'rG:dG,13': 0.42105263200000004,
    'rG:dG,14': 0.428571429,
    'rG:dG,15': 0.272727273,
    'rG:dG,16': 0.0,
    'rG:dG,17': 0.235294118,
    'rG:dG,18': 0.47619047600000003,
    'rG:dG,19': 0.448275862,
    'rG:dG,2': 0.692307692,
    'rG:dG,20': 0.428571429,
    'rG:dG,3': 0.384615385,
    'rG:dG,4': 0.529411765,
    'rG:dG,5': 0.7857142859999999,
    'rG:dG,6': 0.681818182,
    'rG:dG,7': 0.6875,
    'rG:dG,8': 0.615384615,
    'rG:dG,9': 0.538461538,
    'rG:dT,1': 0.9,
    'rG:dT,10': 0.933333333,
    'rG:dT,11': 1.0,
    'rG:dT,12': 0.933333333,
    'rG:dT,13': 0.923076923,
    'rG:dT,14': 0.75,
    'rG:dT,15': 0.9411764709999999,
    'rG:dT,16': 1.0,
    'rG:dT,17': 0.933333333,
    'rG:dT,18': 0.692307692,
    'rG:dT,19': 0.7142857140000001,
    'rG:dT,2': 0.846153846,
    'rG:dT,20': 0.9375,
    'rG:dT,3': 0.75,
    'rG:dT,4': 0.9,
    'rG:dT,5': 0.8666666670000001,
    'rG:dT,6': 1.0,
    'rG:dT,7': 1.0,
    'rG:dT,8': 1.0,
    'rG:dT,9': 0.642857143,
    'rU:dC,1': 0.956521739,
    'rU:dC,10': 0.5,
    'rU:dC,11': 0.4,
    'rU:dC,12': 0.5,
    'rU:dC,13': 0.260869565,
    'rU:dC,14': 0.0,
    'rU:dC,15': 0.05,
    'rU:dC,16': 0.346153846,
    'rU:dC,17': 0.117647059,
    'rU:dC,18': 0.333333333,
    'rU:dC,19': 0.25,
    'rU:dC,2': 0.84,
    'rU:dC,20': 0.176470588,
    'rU:dC,3': 0.5,
    'rU:dC,4': 0.625,
    'rU:dC,5': 0.64,
    'rU:dC,6': 0.571428571,
    'rU:dC,7': 0.588235294,
    'rU:dC,8': 0.7333333329999999,
    'rU:dC,9': 0.6190476189999999,
    'rU:dG,1': 0.857142857,
    'rU:dG,10': 0.533333333,
    'rU:dG,11': 0.666666667,
    'rU:dG,12': 0.947368421,
    'rU:dG,13': 0.7894736840000001,
    'rU:dG,14': 0.28571428600000004,
    'rU:dG,15': 0.272727273,
    'rU:dG,16': 0.666666667,
    'rU:dG,17': 0.705882353,
    'rU:dG,18': 0.428571429,
    'rU:dG,19': 0.275862069,
    'rU:dG,2': 0.857142857,
    'rU:dG,20': 0.090909091,
    'rU:dG,3': 0.428571429,
    'rU:dG,4': 0.647058824,
    'rU:dG,5': 1.0,
    'rU:dG,6': 0.9090909090000001,
    'rU:dG,7': 0.6875,
    'rU:dG,8': 1.0,
    'rU:dG,9': 0.923076923,
    'rU:dT,1': 1.0,
    'rU:dT,10': 0.857142857,
    'rU:dT,11': 0.75,
    'rU:dT,12': 0.8,
    'rU:dT,13': 0.692307692,
    'rU:dT,14': 0.6190476189999999,
    'rU:dT,15': 0.578947368,
    'rU:dT,16': 0.9090909090000001,
    'rU:dT,17': 0.533333333,
    'rU:dT,18': 0.666666667,
    'rU:dT,19': 0.28571428600000004,
    'rU:dT,2': 0.846153846,
    'rU:dT,20': 0.5625,
    'rU:dT,3': 0.7142857140000001,
    'rU:dT,4': 0.47619047600000003,
    'rU:dT,5': 0.5,
    'rU:dT,6': 0.8666666670000001,
    'rU:dT,7': 0.875,
    'rU:dT,8': 0.8,
    'rU:dT,9': 0.9285714290000001,
}


PAM_SCORES = {
    'AA': 0.0,
    'AC': 0.0,
    'AG': 0.25925925899999996,
    'AT': 0.0,
    'CA': 0.0,
    'CC': 0.0,
    'CG': 0.107142857,
    'CT': 0.0,
    'GA': 0.06944444400000001,
    'GC': 0.022222222000000003,
    'GG': 1.0,
    'GT': 0.016129031999999998,
    'TA': 0.0,
    'TC': 0.0,
    'TG': 0.038961038999999996,
    'TT': 0.0,
}


def cfd_score(guide_protospacer: str, candidate_protospacer: str, candidate_pam: str) -> float:
    """
    Real Doench et al. 2016 CFD score for one candidate site relative to
    one guide. Both protospacers must be 20nt DNA, same orientation --
    i.e. candidate_protospacer/candidate_pam are exactly
    offtarget.OffTargetHit.matched_seq/.pam (already in the reference's
    forward-coordinate, strand-resolved frame).

    Returns a score in [0, 1]; 1.0 = no predicted activity loss, 0.0 =
    predicted inactive at this site.
    """
    guide = guide_protospacer.upper().replace("T", "U")
    candidate = candidate_protospacer.upper().replace("T", "U")
    if len(guide) != 20 or len(candidate) != 20:
        raise ValueError("cfd_score() requires 20nt protospacers")

    score = 1.0
    for i, (g, c) in enumerate(zip(guide, candidate), start=1):
        if g == c:
            continue
        key = f"r{g}:d{_complement(c)},{i}"
        score *= MISMATCH_SCORES[key]
    score *= PAM_SCORES.get(candidate_pam[-2:].upper(), 0.0)
    return score


# ---------------- Self-test ----------------
def _self_test():
    # Reference values computed independently, directly from the original
    # CRISPOR pickle files (mismatch_score.pkl / pam_scores.pkl), not from
    # this module -- see the session notes for the exact commands used.
    cases = [
        # (name, guide, candidate, pam, expected)
        ("perfect match, PAM=AGG",
         "GACAGAAAAGCCCCATCCTT", "GACAGAAAAGCCCCATCCTT", "AGG", 1.0),
        ("real PHACTR1 off-target of the PD1 guide (5 mismatches, PAM=AGG); "
         "see validate_iguide.py for the source",
         "CGACTGGCCAGGGCGCCTGT", "TTCATGGCCAGGGAGCCTGT", "AGG", 0.080383),
        ("single mismatch at position 20, PAM=TGG",
         "GACAGAAAAGCCCCATCCTT", "GACAGAAAAGCCCCATCCTA", "TGG", 0.5625),
        ("perfect protospacer match but a non-canonical NAG-style PAM",
         "GACAGAAAAGCCCCATCCTT", "GACAGAAAAGCCCCATCCTT", "AAG", 0.259259),
    ]
    for name, guide, candidate, pam, expected in cases:
        got = cfd_score(guide, candidate, pam)
        assert abs(got - expected) < 1e-5, f"{name}: expected {expected}, got {got}"

    assert len(MISMATCH_SCORES) == 240, len(MISMATCH_SCORES)
    assert len(PAM_SCORES) == 16, len(PAM_SCORES)

    print("CFD self-test passed (matches the original CRISPOR/Doench et al. data).")
    for name, guide, candidate, pam, expected in cases:
        print(f"  {name}: {cfd_score(guide, candidate, pam):.6f}")


if __name__ == "__main__":
    _self_test()
