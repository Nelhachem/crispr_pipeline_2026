# CRISPR knock-in design pipeline.
#
# Python 3.10 is used deliberately: the real Rule Set 3 on-target model
# (rs3, Doench lab 2022 -- see rule_set_3.py) pins old scikit-learn/
# lightgbm/numpy versions that don't have prebuilt wheels for newer
# Pythons. Inside this image rs3 is installed directly (no subprocess/venv
# bridge needed -- that trick is only for running this toolkit's other
# files under a newer system Python; see rule_set_3.py's docstring).
FROM python:3.10-slim

LABEL org.opencontainers.image.title="crispr-pipeline" \
      org.opencontainers.image.description="Pure-Python CRISPR knock-in design pipeline: real Doench 2016 CFD off-target scoring, real Rule Set 3 on-target scoring, HDR donor building." \
      org.opencontainers.image.licenses="MIT"

# libgomp1: LightGBM's runtime OpenMP dependency.
#
# build-essential/cmake are deliberately NOT installed: verified (by
# actually building this image) that manylinux wheels exist for rs3's
# pinned scikit-learn 1.0.2 / lightgbm 3.3.5 on both amd64 and arm64, so
# nothing needs to compile from source here -- unlike local dev on this
# project's own newer system Python (see rule_set_3.py's docstring for
# why THAT environment needed CMake). If you change rs3's pinned versions
# and hit a "Please install CMake" build error, add back
# `build-essential cmake` here and the CMAKE_POLICY_VERSION_MINIMUM env
# vars from rule_set_3.py's docstring.
RUN apt-get update && apt-get install -y --no-install-recommends \
        libgomp1 \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# Real Rule Set 3 (Doench lab, 2022) -- installed once, at build time, so
# every container run starts with it ready. See rule_set_3.py for what
# this model is and why it isn't reimplemented from scratch.
RUN pip install --no-cache-dir rs3

# The toolkit itself has no other dependencies.
COPY sgrna_design.py offtarget.py cfd_score.py rule_set_3.py \
     knockin_donor.py hdr_donor.py run_pipeline.py score_guides.py \
     validate_iguide.py validate_rule_set_3.py ./
COPY aavs1_locus.fa b2m_locus.fa ./

# Sanity-check the image at build time: every self-test must pass,
# including the real-data regressions, so a broken image never ships.
RUN python3 sgrna_design.py && \
    python3 cfd_score.py && \
    python3 offtarget.py && \
    python3 validate_iguide.py && \
    python3 rule_set_3.py && \
    python3 validate_rule_set_3.py && \
    python3 knockin_donor.py && \
    python3 hdr_donor.py

# Mount your own locus (and any --reference/--cargo-fasta files) into
# /data and pass paths relative to it, e.g.:
#   docker run --rm -v "$(pwd):/data" crispr-pipeline my_locus.fa
# No args at all matches local `python3 run_pipeline.py`: runs the
# built-in demo. Use `docker run --rm crispr-pipeline --help` for options.
WORKDIR /data
ENTRYPOINT ["python3", "/app/run_pipeline.py"]
CMD []
