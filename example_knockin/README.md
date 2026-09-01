# Worked example: knocking EGFP into AAVS1

A complete, runnable knock-in design using only real sequence. Every file
here was fetched from NCBI — nothing is synthetic.

| File | What it is | Source |
|---|---|---|
| `aavs1_target_2kb.fa` | The locus to edit: chr19:55,115,000-55,117,000 (2,001 bp), GRCh38 | NCBI `NC_000019.10` |
| `egfp_cargo.fa` | The insert: the real EGFP CDS, 720 bp, in frame (ATG -> TAA) | NCBI `U55762` (pEGFP-N1), CDS 679-1398 |
| `chr19_flanking_panel.fa` | Off-target screen: 198 kb of chr19 flanking the target, as 2 records | NCBI `NC_000019.10` |
| `chr19_BAD_contains_target_200kb.fa` | A deliberately wrong panel, kept to demonstrate a common mistake | NCBI `NC_000019.10` |

## Run it

```bash
cd example_knockin

python3 ../run_pipeline.py aavs1_target_2kb.fa \
    --cargo-fasta egfp_cargo.fa \
    --cargo-name EGFP \
    --arm-len 500 \
    --reference chr19_flanking_panel.fa \
    --screen-n 20 \
    --top-n 5
```

Or in Docker, from this directory:

```bash
docker run --rm -v "$(pwd):/data" crispr-pipeline:arm64 \
    aavs1_target_2kb.fa --cargo-fasta egfp_cargo.fa --cargo-name EGFP \
    --arm-len 500 --reference chr19_flanking_panel.fa --screen-n 20
```

## What you should see

```
Found 427 candidate guides.
Dropped 195 guide(s) too close to the sequence ends (need 500 bp flanking ...).
232 usable — ranked by the REAL Rule Set 3 score (Doench lab, 2022 ...)

 #  protospacer            PAM  str    cut   GC%     RS3 CFD-spec  exact-OT  rpt-OT
 1  GAGATGGCTCCAGGAAATGG   GGG  +      643    55    1.33     76.2         0       0
 ...
Donor length: 1720 bp      (500 bp arm + 720 bp EGFP + 500 bp arm)
RESULT: ALL CHECKS PASSED
```

Row 3 of that table is `AGAGCTAGCACAGACTAGAG` — the AAVS1 guide actually
used in the Nature 2022 CAR-T trial. Nothing was tuned to place it there.

## The mistake this example is designed to teach

Swap in the deliberately wrong panel and re-run:

```bash
python3 ../run_pipeline.py aavs1_target_2kb.fa \
    --cargo-fasta egfp_cargo.fa --cargo-name EGFP --arm-len 500 \
    --reference chr19_BAD_contains_target_200kb.fa --screen-n 20
```

Every candidate now reports `exact-OT 1` and the run ends in
`RESULT: SOME CHECKS FAILED`.

Nothing is broken. That 200 kb window **physically contains** the 2 kb
target, so each guide's own on-target site is found a second time inside a
`--reference` file — and a perfect match in a reference is, correctly,
treated as a real off-target somewhere else in the genome. The pipeline
only excuses the on-target site within your input locus.

The effect on the numbers is large:

| | Contaminated panel | Correct panel |
|---|---|---|
| `exact-OT` | 1 (false alarm) | 0 |
| `CFD-spec` | 43.3 | 76.2 |
| Result | SOME CHECKS FAILED | ALL CHECKS PASSED |

**Rule: build reference panels from regions flanking your target, never
spanning it.** Here that meant two records — chr19:55,015,000-55,114,999
and chr19:55,117,001-55,215,000.

## Bringing a shortlist over from CHOPCHOP

`chopchop_shortlist.tsv` in this folder is a CHOPCHOP-format results table
(the real columns: `Rank`, `Target sequence`, `Genomic location`, `Strand`,
`GC content (%)`, `Self-complementarity`, `MM0`-`MM3`, `Efficiency`).

Score exactly those guides — no rediscovery:

```bash
python3 ../score_guides.py aavs1_target_2kb.fa \
    --guides-file chopchop_shortlist.tsv \
    --reference chr19_flanking_panel.fa
```

```
 #  protospacer            PAM  str     cut   GC%     RS3 CFD-spec  exact-OT
 1  GAGATGGCTCCAGGAAATGG   GGG  +       643    55    1.33     76.2         0
 2  CCGGAGAGGACCCAGACACG   GGG  +      1332    70    1.30     86.8         0
 3  AGAGCTAGCACAGACTAGAG   AGG  +       997    50    1.21     99.4         0
 4  TATAAGGTGGTCCCAGCTCG   GGG  +       854    55    1.19    100.0         0
```

Two things to notice:

1. The deliberately bogus 5th row in the TSV (a poly-T guide at a fake
   coordinate) is reported as **NOT FOUND** rather than silently dropped.
2. CHOPCHOP ranked `AGAGCTAGCACAGACTAGAG` first; Rule Set 3 ranks it third.
   It has the best CFD-spec of the four, and it is the guide actually used
   in the Nature 2022 CAR-T trial. Weigh all three columns — plus
   CHOPCHOP's genome-wide MM0-MM3, which nothing here reproduces.

Then build the donor for whichever guide you choose:

```bash
python3 ../score_guides.py aavs1_target_2kb.fa \
    --guides AGAGCTAGCACAGACTAGAGAGG \
    --build-donor AGAGCTAGCACAGACTAGAG \
    --cargo-fasta egfp_cargo.fa --cargo-name EGFP --arm-len 500
```

That prints the full 1,720 bp donor sequence, ready to order.
