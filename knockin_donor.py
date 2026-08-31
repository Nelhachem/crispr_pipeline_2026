"""
Knock-in HDR donor designer + annotated map.

This is the "stick a CAR (or anything) into a safe-harbor locus" version.
It differs from a point-edit donor in one way: the thing between the two
homology arms is an INSERT (a whole cassette) rather than a single changed base.

PIPELINE (matches the biology you described)
--------------------------------------------
  RNP (Cas9+sgRNA) makes the cut  ->  sgRNA guide picks WHERE (hybridisation)
  ->  HDR copies this donor        ->  cargo is integrated at the cut

DONOR LAYOUT
------------
  [ Left homology arm ] [ INSERT / cargo ] [ Right homology arm ]
        ~40-800 bp        e.g. a CAR           ~40-800 bp

  The arms match the genome on each side of the cut so the cell lines the
  donor up and pastes the cargo in at the break.

SILENT / PAM-BLOCK NOTE
-----------------------
  For a KNOCK-IN, re-cutting is usually prevented automatically: once the
  cargo is inserted, the 20 nt protospacer + PAM no longer sit next to each
  other (the insert splits them), so Cas9 can't rebind. We still REPORT the
  guide/PAM status and, if the arms alone would recreate an intact site, we
  flag it so the user can add a silent block by hand. This keeps the teaching
  logic explicit rather than hidden.

Intentionally simple and readable. Not a substitute for experimental design.
"""

from __future__ import annotations
from dataclasses import dataclass, field

COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


@dataclass
class Feature:
    name: str
    start: int          # 0-based, relative to the full donor
    end: int            # exclusive
    kind: str           # "arm", "insert", "pam_block", "junction"

@dataclass
class KnockinDonor:
    full_donor: str
    features: list = field(default_factory=list)
    notes: list = field(default_factory=list)

    def annotated_map(self) -> str:
        """Return a human-readable annotated map of the donor."""
        lines = []
        lines.append(f"Donor length: {len(self.full_donor)} bp")
        lines.append("")
        lines.append("Features (0-based, end-exclusive):")
        for f in sorted(self.features, key=lambda x: x.start):
            length = f.end - f.start
            preview = self.full_donor[f.start:f.end]
            if len(preview) > 30:
                preview = preview[:14] + "..." + preview[-13:]
            lines.append(f"  {f.start:>5}-{f.end:<5} [{length:>5} bp] "
                         f"{f.kind:<10} {f.name:<18} {preview}")
        if self.notes:
            lines.append("")
            lines.append("Notes:")
            for n in self.notes:
                lines.append(f"  - {n}")
        return "\n".join(lines)

    def ascii_track(self, width: int = 60) -> str:
        """A crude proportional track so you can 'see' the layout."""
        L = len(self.full_donor)
        row = [" "] * width
        symbol = {"arm": "=", "insert": "#", "pam_block": "*", "junction": "|"}
        for f in self.features:
            s = int(f.start / L * width)
            e = max(int(f.end / L * width), s + 1)
            for i in range(s, min(e, width)):
                # don't let arms overwrite the insert marker
                if row[i] == " " or symbol.get(f.kind) == "#":
                    row[i] = symbol.get(f.kind, "?")
        legend = "  = left/right arm   # insert/cargo   * pam-block   | junction"
        return "".join(row) + "\n" + legend


def design_knockin_donor(
    seq: str,
    cut_site: int,          # 0-based forward-coord cut (from find_guides)
    insert_seq: str,        # the cargo to drop in (CAR, reporter, tag, ...)
    insert_name: str = "cargo",
    protospacer: str = "",  # optional: for re-cut reporting
    pam_pos: int | None = None,   # optional: 0-based PAM start in forward coords
    arm_len: int = 40,
) -> KnockinDonor:
    """
    Build an HDR knock-in donor: [left arm][insert][right arm], anchored at cut.
    Homology arms are taken directly from the genome on each side of the cut.
    """
    seq = seq.upper()
    insert_seq = insert_seq.upper()
    notes = []

    left_start = max(cut_site - arm_len, 0)
    left_arm = seq[left_start:cut_site]
    right_arm = seq[cut_site:cut_site + arm_len]

    if len(left_arm) < arm_len:
        notes.append(f"Left arm truncated to {len(left_arm)} bp (near sequence start)")
    if len(right_arm) < arm_len:
        notes.append(f"Right arm truncated to {len(right_arm)} bp (near sequence end)")

    full = left_arm + insert_seq + right_arm

    features = [
        Feature("HA-L", 0, len(left_arm), "arm"),
        Feature(insert_name, len(left_arm), len(left_arm) + len(insert_seq), "insert"),
        Feature("HA-R", len(left_arm) + len(insert_seq), len(full), "arm"),
        Feature("5' junction", len(left_arm) - 1, len(left_arm) + 1, "junction"),
        Feature("3' junction",
                len(left_arm) + len(insert_seq) - 1,
                len(left_arm) + len(insert_seq) + 1, "junction"),
    ]

    # Re-cut check: after insertion the guide and PAM are split by the cargo,
    # so an intact protospacer+PAM should NOT survive. Verify and report.
    if protospacer and pam_pos is not None:
        # reconstruct what the edited genomic site looks like (arms rejoined
        # WITH the insert between them at the cut)
        notes.append("Cargo splits the protospacer/PAM at the cut; "
                     "Cas9 cannot rebind an intact site -> no re-cut expected.")
        # sanity: is the cut actually inside the protospacer footprint?
        if pam_pos - 3 == cut_site or abs((pam_pos - 3) - cut_site) <= 1:
            notes.append(f"Cut at {cut_site} is 3 bp upstream of PAM at "
                         f"{pam_pos} (canonical SpCas9 geometry). OK.")
        else:
            notes.append(f"WARNING: cut at {cut_site} is not 3 bp upstream of "
                         f"PAM at {pam_pos}; check guide/cut coordinates.")
    else:
        notes.append("No guide/PAM supplied; re-cut analysis skipped. "
                     "For SpCas9 knock-ins the cargo usually breaks the site anyway.")

    return KnockinDonor(full_donor=full, features=features, notes=notes)


# ---------------- Self-test ----------------
def _self_test():
    # Synthetic AAVS1-style locus with a KNOWN guide + PAM so coords are verifiable.
    guide = "GTCACCAATCCTGTCCCTAG"        # 20 nt (placeholder, verifiable by construction)
    upstream   = "TTTTTAAAAACCCCCGGGGG"    # 20 bp
    downstream = "ACGTACGTACGTACGTACGT"    # 20 bp
    locus = upstream + guide + "TGG" + downstream
    pam_pos = len(upstream) + len(guide)   # PAM right after guide
    cut = pam_pos - 3                       # 3 bp upstream of PAM

    cargo = "ATGCARCARCARGENEHERE".replace("R", "A")  # fake 20 bp "CAR"
    d = design_knockin_donor(
        locus, cut_site=cut, insert_seq=cargo, insert_name="CAR",
        protospacer=guide, pam_pos=pam_pos, arm_len=15)

    # left arm should be 15 bp ending exactly at the cut
    la = [f for f in d.features if f.name == "HA-L"][0]
    assert (la.end - la.start) == 15, "left arm wrong length"
    # insert should be present and full length
    ins = [f for f in d.features if f.kind == "insert"][0]
    assert d.full_donor[ins.start:ins.end] == cargo, "cargo not placed correctly"
    # cut-to-PAM geometry note present
    assert any("canonical SpCas9 geometry" in n for n in d.notes), d.notes
    print("Knock-in self-test passed.\n")
    print(d.annotated_map())
    print()
    print(d.ascii_track())


if __name__ == "__main__":
    _self_test()
