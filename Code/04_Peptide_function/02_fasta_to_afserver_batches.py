#!/usr/bin/env python3
import argparse, json, re, sys, math
from pathlib import Path

AA_RE = re.compile(r"[A-Z*]")            # uppercase amino acids + optional '*'
GAP_WS_RE = re.compile(r"[\s\-]")        # strip whitespace and '-' gaps
SANITIZE_NAME_RE = re.compile(r"[^A-Za-z0-9._-]+")

def parse_fasta(path):
    """Yield (header, sequence) from a FASTA file."""
    header, buf = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line: 
                continue
            if line.startswith(">"):
                if header is not None:
                    seq = GAP_WS_RE.sub("", "".join(buf)).upper()
                    if not seq: raise SystemExit(f"Empty sequence for header {header}")
                    yield header, seq
                header = line[1:].strip() or "seq"
                buf = []
            else:
                buf.append(line.strip())
    if header is not None:
        seq = GAP_WS_RE.sub("", "".join(buf)).upper()
        if not seq: raise SystemExit(f"Empty sequence for header {header}")
        yield header, seq

def validate_seq(seq, hdr):
    bad = [c for c in seq if not AA_RE.fullmatch(c)]
    if bad:
        raise SystemExit(f"Unsupported characters in {hdr}: {''.join(sorted(set(bad)))}")

def header_to_job_name(h):
    """
    Build a clean job name from a FASTA header.
    Example: 'Q96GM8-2|TOE1|DCR=1027|frame=2|alt' -> 'Q96GM8-2_TOE1_DCR-1027_frame-2_alt'
    """
    parts = h.split("|")
    parts = [p.replace("=", "-") for p in parts]
    name = "_".join(parts)
    name = SANITIZE_NAME_RE.sub("_", name).strip("_")
    return name or "job"

def make_server_job(name, seq, seed, count, use_template):
    return {
        "name": name,
        "modelSeeds": [str(seed)],
        "sequences": [{
            "proteinChain": {
                "sequence": seq,
                "count": count,
                "useStructureTemplate": bool(use_template)
            }
        }],
        "dialect": "alphafoldserver",
        "version": 1
    }

def main():
    ap = argparse.ArgumentParser(
        description="Split multi-FASTA into AlphaFold Server JSON batches (one job per protein)."
    )
    ap.add_argument("fasta", help="Input multi-FASTA")
    ap.add_argument("--outdir", default="batches", help="Output directory (default: batches)")
    ap.add_argument("--batch-size", type=int, default=20, help="Jobs per JSON (default: 20)")
    ap.add_argument("--prefix", default="batch_", help="Filename prefix (default: batch_)")
    ap.add_argument("--start-index", type=int, default=1, help="First batch index (default: 1)")
    ap.add_argument("--seed", default="256845973", help="modelSeeds value (default: 256845973)")
    ap.add_argument("--count", type=int, default=1, help="proteinChain.count (default: 1)")
    ap.add_argument("--no-template", action="store_true",
                    help="Disable templates (useStructureTemplate=false)")
    args = ap.parse_args()

    Path(args.outdir).mkdir(parents=True, exist_ok=True)

    # Load and validate all entries once to know how many batches we’ll write.
    entries = []
    for hdr, seq in parse_fasta(args.fasta):
        validate_seq(seq, hdr)
        name = header_to_job_name(hdr)
        entries.append((name, seq, hdr))
    if not entries:
        raise SystemExit("No sequences found.")

    n = len(entries)
    bsize = max(1, args.batch_size)
    nbatches = math.ceil(n / bsize)
    pad = max(2, len(str(args.start_index + nbatches - 1)))  # zero-pad for nice sorting

    idx = args.start_index
    written = []
    for i in range(0, n, bsize):
        batch = entries[i : i + bsize]
        jobs = [make_server_job(name, seq, args.seed, args.count, not args.no_template)
                for (name, seq, _hdr) in batch]
        fname = f"{args.prefix}{str(idx).zfill(pad)}.json"
        outpath = Path(args.outdir) / fname
        with open(outpath, "w") as w:
            json.dump(jobs, w, indent=2)
        written.append((fname, len(jobs)))
        idx += 1

    # Summary
    total_jobs = sum(m for _, m in written)
    print(f"Input proteins: {n}")
    for fn, m in written:
        print(f"Wrote {fn}  ({m} jobs)")
    print(f"Total jobs written: {total_jobs} in {len(written)} batch file(s).")

if __name__ == "__main__":
    sys.exit(main())
