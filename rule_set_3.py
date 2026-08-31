"""
Real Rule Set 3 on-target scoring (Doench lab, 2022) via the official
`rs3` package: https://github.com/gpp-rnd/rs3

THIS IS NOT A REIMPLEMENTATION. Rule Set 3 is a trained LightGBM model
with hundreds of fitted parameters — there is no responsible way to
reproduce it from the paper alone. Instead, this module calls the real,
unmodified, official package, so the rest of this toolkit stays
dependency-free by default.

TWO WAYS THIS RUNS
--------------------
1. IN-PROCESS: if `import rs3` succeeds directly (e.g. inside the Docker
   image, which is built on Python 3.10 with rs3 installed system-wide —
   see Dockerfile), this module just calls it. No subprocess, no venv.
2. SUBPROCESS BRIDGE: for local dev on a newer/incompatible system Python
   (this repo was developed on Python 3.13, which `rs3`'s old pinned
   scikit-learn/lightgbm/numpy don't have wheels for), a separate Python
   3.10 virtualenv at ./rs3_venv is used instead, via a subprocess.

Either way, if neither is available, on-target scoring transparently
falls back to sgrna_design.py's heuristic — nothing breaks, and
run_pipeline.py always prints which mode was actually used.

SETTING UP rs3_venv FOR LOCAL DEV (one-time, optional; not needed in Docker)
--------------------------------------------------------------------------------
    brew install python@3.10 cmake libomp
    /opt/homebrew/bin/python3.10 -m venv rs3_venv
    CMAKE_POLICY_VERSION_MINIMUM=3.5 CMAKE_ARGS="-DCMAKE_POLICY_VERSION_MINIMUM=3.5" \\
        rs3_venv/bin/python -m pip install rs3

(The CMake flags exist because rs3's pinned LightGBM version predates
CMake 4's dropped support for old cmake_minimum_required() policies —
unrelated to this toolkit, purely an artifact of an old ML package
meeting a new build toolchain. The Dockerfile hits the same issue and
sets the same flags when installing rs3 during the image build.)

INPUT FORMAT
-------------
Rule Set 3 (sequence model) scores 30nt context sequences: 4nt upstream
of the protospacer + the 20nt protospacer + the 3nt PAM + 3nt downstream.
`context_30mer()` extracts this correctly for either strand from a
sgrna_design.Guide and its source sequence.
"""

from __future__ import annotations
import json
import subprocess
from pathlib import Path

from sgrna_design import revcomp

_HERE = Path(__file__).resolve().parent
VENV_PYTHON = _HERE / "rs3_venv" / "bin" / "python"
_RUNNER = _HERE / "rs3_venv_runner.py"


def _has_inprocess_rs3() -> bool:
    try:
        import rs3.seq  # noqa: F401
        return True
    except ImportError:
        return False


_INPROCESS = _has_inprocess_rs3()


def available() -> bool:
    """True if real Rule Set 3 is usable, in-process or via rs3_venv."""
    return _INPROCESS or VENV_PYTHON.exists()


def mode() -> str | None:
    """Which real Rule Set 3 path is active, for user-facing messages."""
    if _INPROCESS:
        return "in-process (rs3 installed directly)"
    if VENV_PYTHON.exists():
        return "subprocess bridge (rs3_venv)"
    return None


def context_30mer(seq: str, guide) -> str | None:
    """
    The 30nt context Rule Set 3 expects: 4nt upstream + 20nt protospacer +
    3nt PAM + 3nt downstream, in the guide's own 5'->3' reading direction.
    Returns None if `seq` doesn't have enough flanking sequence.
    """
    L = len(seq)
    if guide.strand == "+":
        lo, hi = guide.guide_start_fwd - 4, guide.guide_end_fwd + 6
        if lo < 0 or hi > L:
            return None
        return seq[lo:hi]
    else:
        lo, hi = guide.guide_start_fwd - 6, guide.guide_end_fwd + 4
        if lo < 0 or hi > L:
            return None
        return revcomp(seq[lo:hi])


def predict(context_sequences: list[str]) -> list[float] | None:
    """
    Real Rule Set 3 scores for a batch of 30nt context sequences (batched,
    not per-guide — the model has real load time). Returns None if real
    Rule Set 3 isn't available at all. Raises RuntimeError if it's
    available but scoring fails.
    """
    if _INPROCESS:
        import io
        import contextlib
        from rs3.seq import predict_seq
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):   # rs3 prints progress info
            scores = predict_seq(context_sequences)
        return [float(x) for x in scores]

    if not VENV_PYTHON.exists():
        return None
    proc = subprocess.run(
        [str(VENV_PYTHON), str(_RUNNER)],
        input=json.dumps(context_sequences),
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise RuntimeError(f"rs3 scoring failed:\n{proc.stderr[-2000:]}")
    return json.loads(proc.stdout)


# ---------------- Self-test (only runs meaningfully if real rs3 is reachable) ----------------
def _self_test():
    if not available():
        print("Real Rule Set 3 not available (no in-process rs3, no rs3_venv) — "
              "skipping. This is an optional feature; see this module's docstring "
              "to set it up. sgrna_design.py's heuristic score is used in its place.")
        return

    mode = "in-process (rs3 installed directly)" if _INPROCESS else "subprocess bridge (rs3_venv)"
    # Real AAVS1 guide (Nature 2022 Methods) in its real genomic context —
    # same locus used by validate_iguide.py / aavs1_locus.fa.
    context = "TGGAAGAGCTAGCACAGACTAGAGAGGTAA"
    assert context[4:24] == "AGAGCTAGCACAGACTAGAG"
    scores = predict([context])
    assert scores is not None and len(scores) == 1
    assert isinstance(scores[0], float)
    print(f"Rule Set 3 self-test passed [{mode}]. Real score for the AAVS1 guide "
          f"in its real 30nt context: {scores[0]:.4f}")


if __name__ == "__main__":
    _self_test()
