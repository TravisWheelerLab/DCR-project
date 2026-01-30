#!/usr/bin/env python3
"""
Robust DSSP batch runner for AlphaFold .cif files.

Usage:
  python batch_dssp.py ./CIFs --out ./dssp_out --format mmcif --sleep 0.8 --retries 6
  # add --mkdssp /path/to/mkdssp to enable local fallback
"""
import argparse, time, random, subprocess, shutil
from pathlib import Path
import requests

API_URL = "https://pdb-redo.eu/dssp/do"

def mk_session(total_retries=6, backoff_factor=0.8):
    s = requests.Session()
    # Manual retry loop below; keep the session for keep-alive
    s.headers.update({"User-Agent": "dssp-batch/1.0"})
    s.backoff_factor = backoff_factor
    s.total_retries = total_retries
    return s

def post_dssp(session, cif_path: Path, out_format: str, timeout=180):
    # Send multipart/form-data; field name MUST be 'data'
    files = {
        "data": (cif_path.name, cif_path.read_bytes(), "chemical/x-mmcif"),
    }
    data = {"format": out_format}
    return session.post(API_URL, data=data, files=files, timeout=timeout)

def run_mkdssp(mkdssp_bin: str, cif_path: Path, out_path: Path, out_format: str):
    # mkdssp writes classic DSSP; for "annotated mmCIF" we can’t replicate server’s enrichment.
    if out_format == "mmcif":
        raise RuntimeError("Local mkdssp fallback only supports classic DSSP text output.")
    cmd = [mkdssp_bin, "-i", str(cif_path), "-o", str(out_path)]
    subprocess.run(cmd, check=True)

def process_one(cif_path: Path, out_path: Path, out_format: str, session, mkdssp_bin: str | None,
                retries: int, base_sleep: float, timeout: int):
    # Try server first with exponential backoff + jitter; fallback to mkdssp if requested.
    for attempt in range(1, retries + 1):
        try:
            r = post_dssp(session, cif_path, out_format, timeout=timeout)
            if r.status_code == 200 and r.text.strip():
                out_path.write_text(r.text)
                return "server"
            # Treat typical transient codes as retryable
            if r.status_code in (429, 500, 502, 503, 504):
                raise RuntimeError(f"HTTP {r.status_code}")
            # Non-retryable
            raise RuntimeError(f"HTTP {r.status_code}: {r.text[:200]}")
        except Exception as e:
            if attempt == retries and mkdssp_bin:
                # Local fallback (only classic DSSP)
                if out_format != "dssp":
                    # Write sidecar note and still try dssp text locally
                    local_out = out_path.with_suffix(".dssp")
                else:
                    local_out = out_path
                run_mkdssp(mkdssp_bin, cif_path, local_out, "dssp")
                return "local"
            # backoff then retry
            sleep_s = (base_sleep * (2 ** (attempt - 1))) + random.uniform(0, 0.5)
            time.sleep(sleep_s)
    return "fail"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cif_dir", type=Path)
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument("--glob", default="*.cif")
    ap.add_argument("--format", choices=["dssp", "mmcif"], default="dssp")
    ap.add_argument("--sleep", type=float, default=0.5, help="Base pause between jobs (sec)")
    ap.add_argument("--retries", type=int, default=6, help="Server retries before fallback")
    ap.add_argument("--timeout", type=int, default=180)
    ap.add_argument("--mkdssp", type=str, default=None, help="Path to local mkdssp for fallback")
    args = ap.parse_args()

    files = sorted(args.cif_dir.glob(args.glob))
    if not files:
        raise SystemExit(f"No files matched: {args.cif_dir / args.glob}")

    out_dir = args.out or args.cif_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    ext = ".dssp" if args.format == "dssp" else ".annotated.cif"

    mkdssp_bin = args.mkdssp
    if mkdssp_bin is None:
        # Try to find mkdssp on PATH for convenience
        mkdssp_bin = shutil.which("mkdssp")

    sess = mk_session(total_retries=args.retries)

    server_ok = local_ok = fail = 0
    for i, f in enumerate(files, 1):
        out_path = out_dir / (f.stem + ext)
        print(f"[{i}/{len(files)}] {f.name} -> {out_path.name}")
        where = process_one(f, out_path, args.format, sess, mkdssp_bin,
                            retries=args.retries, base_sleep=args.sleep, timeout=args.timeout)
        if where == "server":
            server_ok += 1
        elif where == "local":
            local_ok += 1
        else:
            fail += 1
        time.sleep(max(0.0, args.sleep))  # polite baseline gap
    print(f"Done. Server: {server_ok}, Local: {local_ok}, Failed: {fail}")

if __name__ == "__main__":
    main()

