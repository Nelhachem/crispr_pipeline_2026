"""
Real-data validation for offtarget.py.

Reproduces a literature-validated off-target site from Supplementary
Table 8 ("Off-target sites identified by iGUIDE") of:

    Nature (2022). https://doi.org/10.1038/s41586-022-05140-y

That table reports a real off-target of the trial's PD1-targeting sgRNA,
found by genome-wide unbiased sequencing (iGUIDE) and independently
confirmed by 50,000x-deep targeted sequencing:

    on-target  : PDCD1,   chr2:+:241858824,  616 reads (GRCh38)
    off-target : PHACTR1, chr6:-:13230036,    89 reads, 5 mismatches

This script scans a 101 bp excerpt of the real genome around that site
(GRCh38 chr6:13,229,980-13,230,080, accession NC_000006.12, fetched once
from NCBI and embedded below so this needs no network access to run) with
the paper's real PD1 sgRNA sequence, and checks that offtarget.py
independently recovers the exact same site: same sequence, same PAM,
same mismatch count, same strand.

It also checks the negative: at the OLD default (max_mismatches=4) this
5-mismatch site is invisible. That gap is why run_pipeline.py's default
was raised to 6 (see offtarget.py's module docstring).
"""

from __future__ import annotations
from offtarget import scan_reference

# GRCh38 chr6:13,229,980-13,230,080 (NC_000006.12), via NCBI efetch.
PHACTR1_WINDOW = (
    "CTGTGGGAGCCCAGGCAGATATGTAAGCCTTAATTGACTCATTTCTGTCTCCTACAGGCTCCCTGGCC"
    "ATGAAGGTCTGCAGGAAGGACTCCTTAGCCATC"
)

# PD1-targeting sgRNA used in the trial (same paper's Methods).
PD1_GUIDE = "CGACTGGCCAGGGCGCCTGT"

# Ground truth, decoded from Supplementary Table 8's alignment string.
EXPECTED_OFFTARGET_SEQ = "TTCATGGCCAGGGAGCCTGT"
EXPECTED_PAM = "AGG"
EXPECTED_MISMATCHES = 5
EXPECTED_STRAND = "-"
# Real Doench et al. 2016 CFD score for this exact site (see cfd_score.py),
# computed independently from the original CRISPOR/Doench pickle data.
EXPECTED_CFD = 0.080383


def _validate():
    hits = scan_reference(PD1_GUIDE, PHACTR1_WINDOW, ref_name="PHACTR1_chr6_window",
                           max_mismatches=6, seed_len=12)
    matches = [h for h in hits if h.matched_seq == EXPECTED_OFFTARGET_SEQ]
    assert matches, (
        "Failed to recover the literature-validated PHACTR1 off-target site "
        "(Supplementary Table 8, https://doi.org/10.1038/s41586-022-05140-y) "
        "from its real genomic sequence — the mismatch scan is broken.")
    hit = matches[0]
    assert hit.pam == EXPECTED_PAM, hit
    assert hit.mismatches == EXPECTED_MISMATCHES, hit
    assert hit.strand == EXPECTED_STRAND, hit
    assert hit.cfd is not None and abs(hit.cfd - EXPECTED_CFD) < 1e-5, hit

    # The whole point: the OLD default (max_mismatches=4) must NOT see this
    # real, validated off-target — that's the evidence for raising it to 6.
    old_default_hits = scan_reference(PD1_GUIDE, PHACTR1_WINDOW, ref_name="w",
                                       max_mismatches=4, seed_len=12)
    old_default_matches = [h for h in old_default_hits if h.matched_seq == EXPECTED_OFFTARGET_SEQ]
    assert not old_default_matches, (
        "This site has 5 mismatches; expected it to be invisible at "
        "max_mismatches=4 — if this fails, the earlier finding no longer holds.")

    print("Real-data validation passed.")
    print(f"  Recovered the exact PHACTR1 off-target from Supplementary Table 8:")
    print(f"    {hit.matched_seq} PAM={hit.pam} strand={hit.strand} "
          f"mismatches={hit.mismatches}  (real: {EXPECTED_MISMATCHES} mismatches, "
          f"89 reads vs 616 on-target)")
    print(f"  Real Doench et al. 2016 CFD score for this site: {hit.cfd:.6f} "
          f"(low, consistent with it being a low-abundance site that needed "
          f"50,000x sequencing to catch).")
    print("  Confirmed invisible at the old default (max_mismatches=4), "
          "found at the current default (max_mismatches=6).")


if __name__ == "__main__":
    _validate()
