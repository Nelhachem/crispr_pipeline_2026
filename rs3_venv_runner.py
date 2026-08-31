"""
Runs INSIDE rs3_venv only. Reads a JSON list of 30nt context sequences from
stdin, calls the real rs3 package's predict_seq(), writes a JSON list of
scores to stdout. Kept separate from the main (dependency-free) codebase --
see rule_set_3.py for the caller.
"""
import sys
import json
import io
import contextlib

def main():
    context_seqs = json.load(sys.stdin)
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):   # rs3 prints progress info; keep it off real stdout
        from rs3.seq import predict_seq
        scores = predict_seq(context_seqs)
    sys.stdout.write(json.dumps([float(x) for x in scores]))

if __name__ == "__main__":
    main()
