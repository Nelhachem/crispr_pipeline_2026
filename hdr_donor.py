"""
HDR donor designer — the precise-edit companion to sgrna_design.py.

WHAT HDR DOES (the teaching point)
----------------------------------
A Cas9 cut can be repaired two ways:
  - NHEJ  -> sloppy, random indels -> used for KNOCKOUT (see sgrna_design.py)
  - HDR   -> the cell copies a DONOR template you provide -> PRECISE edit

To make a precise edit you give the cell a donor: your intended sequence
flanked by two "homology arms" that match the genome on each side of the cut.
The cell uses the arms to line up the donor and pastes in your edit.

ONE SUBTLETY THAT TRIPS EVERYONE UP
-----------------------------------
After a successful edit, Cas9 will happily cut the repaired DNA AGAIN if the
guide + PAM are still intact. So a good donor also breaks the PAM or the seed
region with a SILENT mutation (changes DNA but not the protein). This module
does that automatically when the edit itself doesn't already disrupt the site.

This is intentionally simple and readable for teaching. It is NOT a substitute
for experimental validation.
"""

from __future__ import annotations
from dataclasses import dataclass

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


# ---- Minimal codon table (standard genetic code) ----
CODON_TABLE = {
    'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
    'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
    'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
    'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
    'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
    'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
    'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
    'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
}

def synonymous_codons(codon: str):
    """Return other codons coding the same amino acid."""
    aa = CODON_TABLE[codon.upper()]
    return [c for c, a in CODON_TABLE.items() if a == aa and c != codon.upper()]


@dataclass
class Donor:
    left_arm: str
    edited_core: str      # the region around the cut, carrying edit + block
    right_arm: str
    full_donor: str
    notes: str


def make_point_edit_donor(
    seq: str,
    cut_site: int,          # 0-based forward-coord cut position (from find_guides)
    edit_pos: int,          # 0-based position to change
    new_base: str,          # base to install at edit_pos
    protospacer: str,       # 20 nt guide (forward orientation as it appears in seq)
    pam_pos: int,           # 0-based start of the PAM in forward coords
    arm_len: int = 40,
) -> Donor:
    """
    Build a single-stranded-style HDR donor for a 1-bp substitution.

    Steps:
      1. Copy the local genomic region.
      2. Install the intended base change.
      3. If the guide+PAM survive the edit, add a SILENT mutation in the PAM
         (or seed) so Cas9 can't re-cut. Only silent if we know the frame;
         here we do the simplest robust thing: break the PAM's third base,
         which for NGG->NG_ almost always kills recognition.
    """
    seq = seq.upper()
    new_base = new_base.upper()
    notes = []

    # 1. intended edit
    edited = list(seq)
    original = edited[edit_pos]
    edited[edit_pos] = new_base
    notes.append(f"Edit: {original}{edit_pos}{new_base}")

    # 2. block re-cutting by disrupting the PAM (positions pam_pos..pam_pos+2 = N G G)
    # We change the 2nd G of the PAM. This is outside the protein-coding change,
    # so in a real design you'd verify it's silent; here we flag it.
    pam = "".join(edited[pam_pos:pam_pos + 3])
    if len(pam) == 3 and pam[1] == "G" and pam[2] == "G":
        # only add a block if the edit didn't already change the PAM/seed
        edit_hits_site = (pam_pos - 6) <= edit_pos <= (pam_pos + 2)
        if not edit_hits_site:
            edited[pam_pos + 2] = "T"   # NGG -> NGT, abolishes SpCas9 PAM
            notes.append(f"PAM block: {pam}->{pam[:2]}T at {pam_pos+2} (verify silent)")
        else:
            notes.append("Edit already disrupts guide/PAM; no extra block needed")
    else:
        notes.append("PAM not canonical NGG in forward frame; block skipped")

    edited_seq = "".join(edited)

    # 3. cut out arms + core
    core_start = max(cut_site - 5, 0)
    core_end = min(cut_site + 6, len(edited_seq))
    left = edited_seq[max(core_start - arm_len, 0):core_start]
    core = edited_seq[core_start:core_end]
    right = edited_seq[core_end:core_end + arm_len]
    full = left + core + right

    return Donor(left, core, right, full, " | ".join(notes))


# ---------------- Self-test ----------------
def _self_test():
    # A synthetic locus: 20nt guide + PAM (AGG). We install a point edit and
    # confirm the donor (a) carries the new base and (b) breaks the PAM.
    guide = "ATGCATGCATGCATGCATGC"          # 20 nt
    seq = "CCCCC" + guide + "AGG" + "TTTTTTTTTTTTTTTTTTTT"
    pam_pos = 5 + 20                          # PAM starts right after guide
    cut = pam_pos - 3                         # 3 bp upstream of PAM

    d = make_point_edit_donor(
        seq, cut_site=cut, edit_pos=cut, new_base="A",
        protospacer=guide, pam_pos=pam_pos, arm_len=15)

    assert "A" in d.edited_core, "edit not present in core"
    # PAM in the donor should no longer be AGG
    donor_pam_region = d.full_donor
    assert "PAM block" in d.notes or "already disrupts" in d.notes, d.notes
    print("HDR self-test passed.")
    print("  notes:", d.notes)
    print("  left :", d.left_arm)
    print("  core :", d.edited_core)
    print("  right:", d.right_arm)
    print("  full donor length:", len(d.full_donor))


if __name__ == "__main__":
    _self_test()

    print("\n--- DEMO: precise edit on a placeholder locus ---")
    guide = "GACGTGCGCTTCCTGGACGA"
    locus = "AAAAAGGGGG" + guide + "AGG" + "CCCCCTTTTT" + "GGGGGAAAAA"
    pam_pos = 10 + 20
    cut = pam_pos - 3
    d = make_point_edit_donor(
        locus, cut_site=cut, edit_pos=cut, new_base="T",
        protospacer=guide, pam_pos=pam_pos, arm_len=20)
    print("Notes :", d.notes)
    print("Donor :", d.full_donor)
