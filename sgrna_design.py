"""
Simple, teachable sgRNA designer for SpCas9 (NGG PAM) — knockout-oriented.
Built for a CAR-T knockout workflow (e.g. ABCC3). Readable on purpose.

KEY BIOLOGY THE CODE ENCODES
----------------------------
1. Cas9 ignores gene strand. It binds any 20 nt protospacer + NGG PAM on
   EITHER DNA strand. So guides can be "+" (sense) or "-" (antisense) —
   both are legitimate.
2. The blunt DSB is ALWAYS 3 bp 5' of the PAM, on whichever strand the
   guide sits. We report the cut site in forward-reference coordinates so
   every guide is comparable regardless of strand.
3. For knockout we favour cuts early in the CDS (frameshift kills all
   downstream protein) but not at the very 5' end (reinitiation rescue).
"""

from __future__ import annotations
import re
from dataclasses import dataclass

# ---------- 1. Sequence utilities ----------
COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]

def gc_content(seq: str) -> float:
    seq = seq.upper()
    return 100 * (seq.count("G") + seq.count("C")) / len(seq)


# ---------- 2. Guide model ----------
@dataclass
class Guide:
    protospacer: str   # 20 nt, written 5'->3' as it would appear in the sgRNA
    pam: str           # 3 nt NGG
    strand: str        # "+" (guide on forward strand) or "-" (reverse)
    guide_start_fwd: int   # 0-based start of protospacer in forward coords
    guide_end_fwd: int     # 0-based end (exclusive)
    cut_site_fwd: int      # 0-based position of the blunt cut in forward coords
    gc: float
    score: float


# ---------- 3. Find candidate guides on both strands ----------
def find_guides(seq: str, guide_len: int = 20, pam_regex: str = r"(?=([ACGT]GG))"):
    """
    Scan forward and reverse strands for [20 nt][NGG].
    We use a lookahead so overlapping PAMs are all found.
    Coordinates are always mapped back to the FORWARD strand of `seq`,
    so a '+' guide and a '-' guide can be compared on the same axis.
    """
    seq = seq.upper()
    L = len(seq)
    guides = []

    # --- forward strand ---
    for m in re.finditer(pam_regex, seq):
        pam_start = m.start(1)
        g_start = pam_start - guide_len
        if g_start < 0:
            continue
        proto = seq[g_start:pam_start]
        pam = seq[pam_start:pam_start + 3]
        if set(proto) - set("ACGT"):
            continue
        # cut is 3 bp 5' of PAM => between pam_start-3 and pam_start-4... 
        # blunt cut sits 3 nt into the protospacer from the PAM side.
        cut = pam_start - 3
        guides.append(Guide(proto, pam, "+", g_start, pam_start, cut,
                            gc_content(proto), 0.0))

    # --- reverse strand ---
    rc = revcomp(seq)
    for m in re.finditer(pam_regex, rc):
        pam_start = m.start(1)
        g_start = pam_start - guide_len
        if g_start < 0:
            continue
        proto = rc[g_start:pam_start]          # 5'->3' on the reverse strand
        pam = rc[pam_start:pam_start + 3]
        if set(proto) - set("ACGT"):
            continue
        cut_rc = pam_start - 3                  # cut in reverse-strand coords
        # map reverse-strand indices back to forward coords
        g_start_fwd = L - pam_start
        g_end_fwd   = L - g_start
        cut_fwd     = L - cut_rc
        guides.append(Guide(proto, pam, "-", g_start_fwd, g_end_fwd, cut_fwd,
                            gc_content(proto), 0.0))

    return guides


# ---------- 4. Transparent heuristic scoring ----------
# NOT Doench 2016 — swappable teaching rules. Explainable line by line.
def cds_fraction(g: Guide, cds_start: int, cds_end: int,
                  cds_upstream_bp: int = 0,
                  cds_total_len: int | None = None) -> float | None:
    """
    How far into the CODING SEQUENCE this guide cuts, 0.0 to 1.0.
    None if the cut falls outside the coding window entirely.

    cds_start/cds_end are 0-based, end-exclusive offsets into the input
    sequence (your FASTA), not genome coordinates.

    Most genes are multi-exon, and a FASTA window usually holds only part of
    the CDS. Two optional arguments make the fraction refer to the WHOLE
    coding sequence rather than just the part in this window:

      cds_upstream_bp : coding bases lying before this window in the full CDS
      cds_total_len   : length of the complete CDS across all exons

    Without them the fraction is computed against the window alone, which
    over-states how far in the cut is whenever the window is a single exon.
    """
    window_len = cds_end - cds_start
    if window_len <= 0:
        return None
    offset_in_window = g.cut_site_fwd - cds_start
    if not 0 <= offset_in_window <= window_len:
        return None
    total = cds_total_len if cds_total_len else window_len
    if total <= 0:
        return None
    return (cds_upstream_bp + offset_in_window) / total


def score_guide(g: Guide, cds_start: int | None = None,
                 cds_end: int | None = None) -> Guide:
    """
    Transparent heuristic on-target score. Used only as a FALLBACK when the
    real Rule Set 3 model is unavailable (see rule_set_3.py).

    The knockout position rule is applied ONLY when you declare where the
    coding sequence actually is. Earlier versions derived it from the length
    of the input FASTA, which meant the same guide scored differently
    depending purely on how you cropped the window — a real bug, since the
    window is an arbitrary choice and the biology is not.
    """
    s = 100.0
    if g.gc < 40 or g.gc > 70:          # extreme GC hurts activity/specificity
        s -= 20
    if "TTTT" in g.protospacer:          # poly-T = PolIII terminator
        s -= 40
    if g.protospacer.count("G") > 12 or g.protospacer.count("C") > 12:
        s -= 10
    if g.protospacer[0] == "G":          # U6 likes a 5' G
        s += 5
    # Knockout preference: a frameshift early in the CDS truncates everything
    # downstream, but cutting in the very 5' end can be rescued by translation
    # reinitiation at a downstream ATG. Requires a declared CDS.
    if cds_start is not None and cds_end is not None:
        frac = cds_fraction(g, cds_start, cds_end)
        if frac is not None:
            if 0.05 <= frac <= 0.40:
                s += 15
            elif frac < 0.05:
                s -= 10
    g.score = max(s, 0.0)
    return g


# ---------- 5. Pipeline ----------
def design(seq: str, top_n: int = 10, cds_start: int | None = None,
            cds_end: int | None = None):
    guides = [score_guide(g, cds_start, cds_end) for g in find_guides(seq)]
    guides.sort(key=lambda x: (x.score, -x.cut_site_fwd), reverse=True)
    return guides[:top_n]


# ---------- 6. Self-test: prove the strand/cut logic is correct ----------
def _self_test():
    # Construct a sequence with ONE known forward PAM and ONE known reverse PAM.
    # Forward guide: 20 A's then AGG (PAM). Cut must be 3 bp before PAM.
    fwd = "A"*20 + "AGG" + "T"*10
    gs = find_guides(fwd, pam_regex=r"(?=([ACGT]GG))")
    fwd_hit = [g for g in gs if g.strand == "+" and g.pam.endswith("GG")]
    assert any(g.protospacer == "A"*20 for g in fwd_hit), "forward protospacer wrong"
    fg = [g for g in fwd_hit if g.protospacer == "A"*20][0]
    # PAM starts at index 20, cut = 17
    assert fg.cut_site_fwd == 17, f"forward cut wrong: {fg.cut_site_fwd}"

    # Reverse: put a CCN on the forward strand => NGG on reverse.
    # "CCT" on forward => reverse strand reads "AGG" as PAM.
    rev = "G"*10 + "CCT" + "C"*20
    gs2 = find_guides(rev, pam_regex=r"(?=([ACGT]GG))")
    rev_hit = [g for g in gs2 if g.strand == "-"]
    assert rev_hit, "no reverse guide found"
    print("Self-test passed: forward cut@17, reverse guides found:", len(rev_hit))

if __name__ == "__main__":
    _self_test()

    # ---- Demo on a PLACEHOLDER sequence (replace with real ABCC3 exon FASTA) ----
    # >>> This is NOT real ABCC3 sequence. Paste the verified exon here. <<<
    PLACEHOLDER_ABCC3_EXON = (
        "ATGGCCTGCGCCTTCCAGGGCCTGCTGCTGCTGGGCACCCTGCTGCGCGGC"
        "GCCTGGGCCGAGCCCGGCGGCGGCAGCGAGCTGGACGTGCGCTTCCTGGAC"
        "GAGGGCACCCTGCGCCTGGGCGGCTTCTGCGACCTGCTGCGCGCCGTGGGC"
        "CAGGGCCTGCAGCCCGGCGACCTGCTGGCCGTGGTGGGCCCCGTGGGCTGC"
    )
    print("\n--- DEMO (placeholder sequence, not real ABCC3) ---")
    header = "protospacer(5'-3')"
    print(f"{'#':>2} {header:<22} {'PAM':<4} {'str':<3} "
          f"{'cut_fwd':>7} {'GC%':>5} {'score':>5}")
    # The whole placeholder IS the coding sequence here, so the CDS window is
    # the whole string. In real use you must say where the CDS actually sits.
    for i, g in enumerate(design(PLACEHOLDER_ABCC3_EXON, top_n=8,
                                 cds_start=0,
                                 cds_end=len(PLACEHOLDER_ABCC3_EXON)), 1):
        print(f"{i:>2} {g.protospacer:<22} {g.pam:<4} {g.strand:<3} "
              f"{g.cut_site_fwd:>7} {g.gc:>5.0f} {g.score:>5.0f}")
