# Simple CRISPR design toolkit

A small, readable sgRNA + HDR donor designer for teaching. Inspired by CHOPCHOP's
core logic, kept deliberately simple so it can be walked through line by line.

**No dependencies.** Pure Python 3.7+. Nothing to install.

## Files

| File | What it does |
|---|---|
| `sgrna_design.py` | Finds SpCas9 guides (NGG PAM) on both strands, scores on-target activity (heuristic fallback) |
| `rule_set_3.py` | Bridges to the REAL, published Rule Set 3 on-target model (Doench lab, 2022) — see below |
| `offtarget.py` | Off-target search: scans reference sequence(s) for near-matches |
| `cfd_score.py` | The real, published Doench et al. 2016 CFD off-target score (not a heuristic) |
| `validate_iguide.py` | Real-data regression test: reproduces a literature-validated off-target site + its exact CFD score |
| `validate_rule_set_3.py` | Real-data regression test: the real Rule Set 3 model ranks two published guides #1 |
| `knockin_donor.py` | Builds an HDR donor for inserting cargo (a CAR, tag, reporter) |
| `hdr_donor.py` | Builds an HDR donor for a precise point edit (with silent PAM block) |
| `run_pipeline.py` | Runs the whole thing end to end: find guides -> off-target screen -> pick best -> build donor -> verify |
| `score_guides.py` | Scores a **specific** guide list (e.g. a CHOPCHOP shortlist) instead of discovering its own — see below |
| `Dockerfile` | Builds a ~1GB image with real Rule Set 3 pre-installed — see "Running it as a Docker image" below |

## Quick start

```bash
# Run the built-in demo (synthetic locus)
python3 run_pipeline.py

# Run on a real locus bundled in this repo (GRCh38, fetched from NCBI):
#   aavs1_locus.fa - chr19:55,115,796-55,116,196, the AAVS1 safe-harbor site
#   b2m_locus.fa   - chr15:44,711,477-44,711,653, B2M exon 1 (+40bp arms)
python3 run_pipeline.py aavs1_locus.fa
python3 run_pipeline.py b2m_locus.fa

# Run on your own sequence
python3 run_pipeline.py my_locus.fa

# Real design: screen candidates against a reference panel (paralogs,
# pseudogenes, a chromosome arm) and supply your actual cargo sequence
python3 run_pipeline.py my_locus.fa \
    --reference paralogs.fa pseudogenes.fa \
    --cargo-fasta my_cargo.fa --cargo-name MYINSERT

# Full option list
python3 run_pipeline.py --help
```

## Running it as a Docker image

The easiest way to hand this to someone else, or run it on a machine you
don't want to configure by hand: build the image once and it ships with
the real Rule Set 3 model already installed — no CMake fights, no Python
version juggling (that's all resolved at image-build time, on the image's
own Python 3.10; see the Dockerfile and `rule_set_3.py`'s docstring for
why that's necessary at all).

```bash
# Build once (~1 min if wheels are cached, ~3-4 min cold; ~1GB image).
# Every self-test — including the real-data regressions — runs as part
# of the build, so a broken image is refused before it's ever tagged.
docker build -t crispr-pipeline .

# Run the built-in demo
docker run --rm crispr-pipeline

# Run one of the bundled real loci
docker run --rm crispr-pipeline /app/aavs1_locus.fa

# Run on YOUR locus of interest: mount the directory containing it to
# /data (that's the container's working directory), then pass a path
# relative to it. Any --reference / --cargo-fasta files just need to be
# in the same mounted directory.
docker run --rm -v "$(pwd):/data" crispr-pipeline my_locus.fa
docker run --rm -v "$(pwd):/data" crispr-pipeline my_locus.fa \
    --reference paralogs.fa --cargo-fasta my_cargo.fa --cargo-name MYINSERT

# Full option list
docker run --rm crispr-pipeline --help
```

Verified end to end, not just written: built the image, confirmed all 8
self-tests pass during the build (`rule_set_3.py`/`validate_rule_set_3.py`
reproduce the same real AAVS1-ranks-#1/71 and B2M-ranks-#1/11 result
inside the container), then ran it against a locus mounted from an
arbitrary host directory the way a real user would.

Two things worth knowing:
- Inside the container, `rs3` (real Rule Set 3) is installed **directly**
  into the image's Python 3.10 — no subprocess/venv bridge needed there,
  unlike local dev on this project's own (newer, incompatible) system
  Python. `rule_set_3.py` detects which mode it's in and both work.
- `build-essential`/`cmake` are deliberately **not** in the final image:
  verified that prebuilt manylinux wheels exist for `rs3`'s pinned
  scikit-learn/lightgbm on both amd64 and arm64, so nothing compiles from
  source in Docker (only local dev on an unsupported Python needed that).
  If you ever bump `rs3`'s pinned versions and lose wheel coverage, the
  Dockerfile comments say what to add back.

### Exporting the image for a collaborator

If they don't want to (or can't) build it themselves, export a portable
tarball instead — no source code, no Docker Hub account, no internet
access needed on their end beyond `docker load`.

```bash
# Build for their CPU architecture. Apple Silicon Mac = arm64;
# almost everything else (Intel/AMD Mac, Windows, most cloud/Linux) = amd64.
docker build --platform linux/arm64 -t crispr-pipeline:arm64 .
docker build --platform linux/amd64 -t crispr-pipeline:amd64 .   # builds via emulation on Apple Silicon

# Export + compress (~215-225MB each)
docker save crispr-pipeline:arm64 | gzip > crispr-pipeline-arm64.tar.gz
docker save crispr-pipeline:amd64 | gzip > crispr-pipeline-amd64.tar.gz
```

Send them the one file matching their machine. They load and run it with:

```bash
gunzip -c crispr-pipeline-arm64.tar.gz | docker load   # or the amd64 file
docker run --rm crispr-pipeline:arm64                   # built-in demo
docker run --rm -v "$(pwd):/data" crispr-pipeline:arm64 my_locus.fa
```

Verified, not assumed: built both architectures (amd64 via QEMU emulation
on this Mac), removed the local image tags entirely, reloaded each
**purely from its exported `.tar.gz`**, and ran both the built-in demo and
a mounted real locus against the freshly-loaded image — the same
experience a collaborator would actually have. This also caught a real
bug: the image's default (no-args) behavior was showing `--help` instead
of running the demo, which didn't match `python3 run_pipeline.py`'s local
behavior or this README's own claim above — fixed in the Dockerfile
(`CMD []` instead of `CMD ["--help"]`) and reverified before either
tarball was produced.

`aavs1_locus.fa` and `b2m_locus.fa` are both real GRCh38 sequence pulled
from NCBI (not synthetic) and both contain literature-validated guides —
see "Validated against real published guides" below for what was checked
and where each guide sits.

Each module also self-tests when run directly:

```bash
python3 sgrna_design.py         # verifies strand + cut-site geometry
python3 cfd_score.py            # verifies the real CFD table against known reference values
python3 offtarget.py            # verifies off-target search + specificity scoring
python3 validate_iguide.py      # reproduces a real, literature-validated off-target site + its CFD score
python3 rule_set_3.py           # (optional) proves the real Rule Set 3 bridge works, if rs3_venv is set up
python3 validate_rule_set_3.py  # (optional) the real model ranks 2 published guides #1 — see below
python3 knockin_donor.py        # verifies donor assembly
python3 hdr_donor.py            # verifies point-edit donor + PAM block
```

## Knock-out and knock-in

The pipeline does both, selected with `--mode` (default `knockin`).

```bash
# Knock-in: insert cargo by HDR, get a donor sequence back
python3 run_pipeline.py my_locus.fa --cargo-fasta cargo.fa --arm-len 500

# Knock-out: disrupt the gene by NHEJ frameshift. No donor is designed.
python3 run_pipeline.py b2m_locus.fa --mode knockout \
    --cds-start 70 --cds-end 137 --cds-total-len 360
```

They differ in more than the donor:

| | Knock-in | Knock-out |
|---|---|---|
| Repair pathway | HDR, needs a donor template | NHEJ, no template |
| Guides filtered by | Room for homology arms at both ends | Landing inside the coding sequence |
| Extra ranking signal | none | Position within the CDS |
| Output | Donor sequence to order | sgRNA only, plus a genotyping plan |
| Re-cut prevention | Automatic, the cargo splits guide from PAM | Not applicable |

**Declaring the CDS.** Knock-out ranking rewards cuts in roughly the first
5 to 60 percent of the coding sequence: far enough in that translation
cannot reinitiate at a downstream ATG, early enough that a frameshift
truncates most of the protein. That needs to know where the CDS actually
is, which a bare FASTA cannot tell you, so you declare it:

- `--cds-start` / `--cds-end` are 0-based offsets **into your FASTA**, not
  genome coordinates.
- `--cds-total-len` is the length of the complete CDS across all exons, and
  `--cds-upstream-bp` is how many coding bases precede this window.

Those last two matter more than they look. Most genes are multi-exon and a
FASTA window usually holds one exon, so without them the position is
measured against the window alone. In the B2M example above the real
published guide sits at **6.1 percent** of the 360 bp CDS; omit
`--cds-total-len` and the same guide reports **32.8 percent**, because it is
then being measured against just the 67 bp of exon 1. The pipeline prints a
warning whenever it is making that assumption.

If you supply no CDS at all, knock-out mode still runs, but it says plainly
that it cannot rank by position and skips those checks.

## Using this alongside CHOPCHOP or CRISPOR (recommended)

This toolkit has no genome index, so it cannot do genome-wide off-target
search — that is its single biggest limitation. Tools like
[CHOPCHOP](https://chopchop.cbu.uib.no/) and CRISPOR can, across 700+
pre-indexed organisms. What they do *not* offer is Rule Set 3, the 2022
Doench-lab on-target model (CHOPCHOP's newest is Doench 2016).

So the two compose well:

| Step | Tool | Why |
|---|---|---|
| 1. Find guides, screen genome-wide | CHOPCHOP / CRISPOR | Real genome index; gives MM0-MM3 off-target counts |
| 2. Re-score the shortlist | `score_guides.py` here | Rule Set 3 (2022) + real CFD specificity |
| 3. Build the donor | `score_guides.py --build-donor` | Emits a ready-to-order HDR donor sequence |

`score_guides.py` exists specifically to close this handoff. Unlike
`run_pipeline.py`, which discovers its own guides, it scores **exactly the
guides you give it** — so a shortlist from another tool can be carried
across without eyeballing two tables.

```bash
# 1. Export your CHOPCHOP results table (it has a "Target sequence" column).
#    Fetch the same locus as FASTA — coordinates come from CHOPCHOP's
#    "Genomic location" column.

# 2. Score that shortlist here. 23 nt sequences (protospacer+PAM) are
#    trimmed for you; a plain one-per-line list works too.
python3 score_guides.py my_locus.fa \
    --guides-file chopchop_results.tsv \
    --reference flanking_panel.fa

# 3. Build the donor for whichever guide wins on the combined evidence.
python3 score_guides.py my_locus.fa \
    --guides AGAGCTAGCACAGACTAGAGAGG \
    --build-donor AGAGCTAGCACAGACTAGAG \
    --cargo-fasta egfp_cargo.fa --cargo-name EGFP --arm-len 500
```

**How to weigh the two sets of numbers.** They answer different questions,
so use each for what it is good at:

- **CHOPCHOP MM0-MM3** — is this guide unique *in the genome*? Treat a
  non-zero MM0 (beyond the on-target site itself) or a high MM1/MM2 as
  disqualifying. Nothing here can tell you this.
- **RS3** — how active is it likely to be? Primary ranking key among
  guides that already passed the genome-wide check.
- **CFD-spec** — how risky are the near-matches in the region you care
  about? Good for breaking ties, and for spotting local repeats.

Expect the two tools to disagree on ranking; that is the point of running
both. In `example_knockin/` CHOPCHOP-style ranking put
`AGAGCTAGCACAGACTAGAG` first while Rule Set 3 put it third — but it had the
cleanest CFD-spec (99.4) and is the guide actually used in the Nature 2022
trial. Disagreement is information, not a fault.

A worked end-to-end version of this lives in
[`example_knockin/`](example_knockin/).

## Off-target screening

`run_pipeline.py` now screens the top on-target candidates (`--screen-n`,
default 15) before picking one: it scans every NGG/NAG PAM site on both
strands of the input locus, plus any `--reference` FASTA files you supply,
for near-matches to each candidate protospacer. Candidates with an
exact-match off-target elsewhere are excluded by default (override with
`--allow-exact-repeats`); the rest are ranked by on-target score with a
**CFD specificity score** shown alongside (`CFD-spec` column).

This is a **real search** — it finds actual matching sites in whatever
sequence(s) you give it — but it is reference-scoped, not genome-wide: with
no `--reference` supplied it only screens the input locus itself (catches
tandem repeats/duplications, not genome-wide off-targets). For a real
design, point `--reference` at known paralogs/pseudogenes or a chromosome
region; for therapeutic-grade specificity, use a genome-scale tool
(Cas-OFFinder, CRISPOR, FlashFry) and validate experimentally (GUIDE-seq,
CIRCLE-seq) before trusting any guide. See the module docstring in
`offtarget.py` for the scoring rationale and its limits.

## The real published off-target score (CFD, Doench et al. 2016)

Earlier versions of this scoring were a from-scratch heuristic. They no
longer are, for off-target risk: `cfd_score.py` implements the actual
**Cutting Frequency Determination (CFD) score** from Doench, Fusi et al.
2016 (*Nat Biotechnol* 34:184-191, https://doi.org/10.1038/nbt.3437) — a
per-position mismatch penalty matrix (240 entries) and a PAM score table
(16 entries), fitted from real experimental cutting-frequency data. These
tables were transcribed **programmatically**, not retyped by hand, from
[CRISPOR's public data files](https://github.com/maximilianh/crisporWebsite/tree/master/CFD_Scoring)
(Haeussler et al.) — the standard open-source reproduction of the original
supplementary data — and the scoring function is a line-for-line port of
CRISPOR's reference implementation, verified against it.

Every `OffTargetHit` now carries a real `.cfd` score (0-1; 1.0 = no
predicted activity loss). `summarize()` reports `cfd_specificity_score`,
a guide-level aggregate (sum the real per-site CFD scores, then invert) —
the per-site number is the published model; summing many of them into one
guide-level number is still this toolkit's own convention, since Doench
et al. report per-site scores rather than a single canonical aggregate.
The old from-scratch heuristic is kept as `specificity_score()` for
comparison, clearly labeled legacy.

## The real published on-target score (Rule Set 3, Doench lab 2022)

Off-target risk uses real published data (above); on-target activity now
can too. Rule Set 3 is a trained LightGBM model with hundreds of fitted
parameters — unlike CFD, there's no responsible way to reimplement it
from the paper alone (there is no small public coefficient table; it's a
gradient-boosted tree ensemble). So `rule_set_3.py` doesn't reimplement
anything: it calls the **real, unmodified, official** package —
[rs3](https://github.com/gpp-rnd/rs3) (Doench lab, Broad Institute
Genetic Perturbation Platform) — in a dedicated subprocess.

**Why a subprocess:** `rs3` pins old scikit-learn/LightGBM/numpy versions
that don't have prebuilt wheels for a current Python, so getting it
running took a separate Python 3.10 environment. Rather than force that
onto everyone using this toolkit, `run_pipeline.py` looks for a
pre-built venv at `./rs3_venv` and uses it automatically when present;
if it's not there, on-target ranking transparently falls back to
`sgrna_design.py`'s heuristic — nothing breaks either way, and
`run_pipeline.py` always prints which one it used.

**One-time optional setup** (macOS/Homebrew; ~5 min, ~150MB):

```bash
brew install python@3.10 cmake libomp
/opt/homebrew/bin/python3.10 -m venv rs3_venv
CMAKE_POLICY_VERSION_MINIMUM=3.5 CMAKE_ARGS="-DCMAKE_POLICY_VERSION_MINIMUM=3.5" \
    rs3_venv/bin/python -m pip install rs3
```

(The `CMAKE_POLICY_VERSION_MINIMUM` flag exists because `rs3`'s pinned
LightGBM version predates CMake 4's removal of old
`cmake_minimum_required()` policies — an artifact of an old ML package
meeting a new build toolchain, unrelated to this project.)

`rs3_venv/` is a local build artifact tied to this machine (compiled
LightGBM binaries) — don't commit it or copy it elsewhere; rerun the
three commands above on any other machine that wants real Rule Set 3
scoring.

**Once set up, the result is real and independently striking:** on
`aavs1_locus.fa`, the actual Rule Set 3 model ranks the real AAVS1 guide
used in the Nature 2022 CAR-T trial **#1 of 71 candidates** — the same
guide the from-scratch heuristic had ranked 53rd of 91. On
`b2m_locus.fa`, it ranks the real B2M guide from the same paper **#1 of
11**. Both are pinned as a regression test in `validate_rule_set_3.py`.
This is strong, if informal, evidence that the real published model
tracks what researchers actually chose far better than a hand-written
heuristic — exactly the reason to prefer it when it's available.

### Scaling to real human reference sequences

Two things matter once `--reference` points at a real, non-trivial genomic
file (a chromosome arm, a multi-record panel of paralogs/pseudogenes)
instead of a toy sequence:

- **Speed.** Each reference is indexed once and reused for every screened
  candidate, instead of re-scanning per guide. References under
  `--fast-index-threshold` (default 2,000,000 bp) use an exhaustive
  per-guide scan; larger ones automatically switch to a seed-indexed lookup
  (`build_index()`/`scan_with_index()` in `offtarget.py`) that's built once
  and then near-O(1) per guide. A synthetic 3 Mb reference screens 15
  candidates in well under a second in pure Python. The trade-off at that
  scale: only seed-exact hits are found — but an exact 20nt match always
  has an exact seed, so the exact-off-target exclusion gate is exactly as
  accurate as at small scale; only lower-risk, seed-mismatched partial
  matches go unenumerated.
- **Repeat-masking.** Real human genome downloads (UCSC, Ensembl, NCBI)
  commonly soft-mask repetitive/low-complexity sequence as lowercase
  (the RepeatMasker convention). `--reference` files (and the input locus)
  are read with this preserved (`read_masked_fasta_records()`); any
  off-target hit — or the chosen guide's own site — landing in masked
  sequence is flagged (`rpt-OT` column, and a note if the on-target site
  itself is masked). A guide sitting in repetitive sequence is a red flag
  on its own: by definition it occurs elsewhere in the genome too. (Note:
  a plain NCBI `efetch` FASTA, like the one used to build `aavs1_locus.fa`
  in this repo, typically has no soft-masking at all — pull from UCSC if
  you want this signal.)

## Validated against a real off-target (not just synthetic tests)

`validate_iguide.py` reproduces a real off-target site from Supplementary
Table 8 of a Nature paper on CRISPR-engineered CAR-T cells
(https://doi.org/10.1038/s41586-022-05140-y): a PD1-targeting sgRNA has a
genuine, 50,000x-deep-sequencing-confirmed off-target inside *PHACTR1*
(chr6, GRCh38), found originally by iGUIDE (an unbiased genome-wide assay).
Given the real genomic sequence around that site, `offtarget.py`
independently recovers the exact same 20nt sequence, PAM, mismatch count
(5), and strand reported in the paper — and `cfd_score.py` computes a real
CFD score of **0.080** for it: low, consistent with it being the
low-abundance (89 vs 616 reads) site that needed 50,000x sequencing to
catch in the first place.

That off-target has **5 mismatches**. This is why `--max-mismatches`
defaults to **6**, not 4 — at the old default this real, published,
wet-lab-confirmed off-target would have been invisible to a screen. If
you tighten `--max-mismatches` for speed on a very large reference, know
that you may be trading away exactly this kind of catch.

### Validated against real published guides

Ad hoc checks (not scripted into the test suite, run manually against
real GRCh38 sequence pulled from NCBI):

| Locus | Guide source | What matched |
|---|---|---|
| AAVS1, chr19 | [Nature 2022](https://doi.org/10.1038/s41586-022-05140-y) Methods | Exact protospacer/PAM/coordinates |
| AAVS1, chr19 (different site) | [IDT Alt-R HDR positive-control sheet](https://sfvideo.blob.core.windows.net/sitefinity/docs/default-source/flyer/using-crispr-cas9-hdr-positive-controls-product-sheet.pdf) | `knockin_donor.py`, given matching 40bp arms, reconstructed IDT's commercial donor **byte-for-byte** |
| TRAC, chr14 | [STAR Protocols 2022](https://doi.org/10.1016/j.xpro.2021.101031) (guide originally from Ren et al. 2017) | Exact protospacer/PAM(CGG)/strand |
| B2M, chr15 | Nature 2022 Methods | Exact protospacer/PAM(AGG); cut site confirmed via NCBI's gene table to land inside **exon 1**'s coding sequence (44,711,547-44,711,613), at the paper's reported coordinate |

These are single-source checks (only `validate_iguide.py`'s PHACTR1 case
is cross-checked against primary supplementary data and pinned as a
regression test) — treat them as spot checks on real sequence, not a
substitute for wet-lab validation of any specific guide.

`b2m_locus.fa` is deliberately scoped tight to exon 1 (+40bp flanking on
each side, matching `--arm-len`'s default) rather than the whole 6.6kb
B2M gene. On the full gene, `run_pipeline.py`'s on-target heuristic has no
concept of exon/intron boundaries and happily picked its top candidate
from an intron — a real limitation worth knowing about: for exon-targeted
edits, scope your input FASTA to the exon (+ enough flanking for arms),
don't hand it a whole gene and assume the top-ranked guide lands somewhere
functionally meaningful.

## FASTA format

A plain single-record FASTA. Non-ACGTN characters are rejected.

```
>my_locus description here
GGGGCCACTAGGGACAGGATTGGTGACAGAAAAGCCCCATCCTTAGGCCTC
CTCCTTCCTAGTCTCCTGATATTGGGTCTAACCCCCACCTCCTGTTAGGCA
```

Give at least ~100 bp of flanking sequence on each side of your intended cut,
or the homology arms get truncated (the pipeline warns you and drops guides
that are too close to the ends).

## The biology it encodes

1. **RNP makes the cut.** Cas9 protein + sgRNA.
2. **The sgRNA picks where** by hybridising to a 20 nt protospacer.
3. **A PAM (NGG) must sit immediately 3' of the protospacer** — no PAM, no cut.
4. **The blunt cut is always 3 bp 5' of the PAM.**
5. **Cas9 reads either DNA strand** — guides can be sense or antisense. All
   coordinates are reported on the forward strand so they're comparable.
6. **HDR copies your donor** into the break: `[left arm][cargo][right arm]`.
7. **Re-cutting is prevented** — for knock-ins the cargo splits the protospacer
   from its PAM automatically; for point edits a silent PAM mutation is added.

## Important caveats

- **Both scores can now be the real, published models**: off-target
  specificity always uses the real Doench et al. 2016 CFD score
  (`cfd_score.py`); on-target activity uses the real Rule Set 3 model
  (Doench lab 2022, `rule_set_3.py`) whenever `rs3_venv` is set up (see
  above), otherwise it falls back to a documented from-scratch heuristic
  and says so explicitly in the output. Check `run_pipeline.py`'s STEP 1
  output line to see which one a given run actually used.
- Off-target screening is **reference-scoped, not genome-wide** — it finds
  real matches in whatever sequence(s) you supply via `--reference`, but does
  not ship a whole-genome index. For genome-scale specificity use an aligner
  (Bowtie/BWA) or a dedicated tool (Cas-OFFinder, CRISPOR, FlashFry), and
  validate experimentally (e.g. GUIDE-seq) before trusting any guide.
- The built-in demo locus is **synthetic**, not a real genomic sequence. Supply
  verified sequence from Ensembl or NCBI for real designs.
