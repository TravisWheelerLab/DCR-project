#!/usr/bin/env python3
"""stage_selection_inputs.py -- collect and verify the 11 intermediate files read by selection_v3.py.

Purpose  The alternate-frame selection test (Figure 2B(iii)) starts from intermediates built during
         the July-August 2026 revision session (liftOver of DCR / latent-null intervals to mm39,
         GENCODE v48 CDS intersections, protein-derived DCR frames). Their generating code is kept
         verbatim under ../archive_prior_session/ (as executed, not refactored). This script finds
         the files under one or more search directories (exported artifact folders may prefix file
         names with 'v<8 hex>_'), checks size + md5 against selection_inputs_manifest.tsv, and
         symlinks (default) or copies them into --dest under their canonical names.
Inputs   --search DIR [DIR ...]   directories searched recursively
         --manifest               default: selection_inputs_manifest.tsv next to this script
Outputs  --dest/<file> for every manifest entry (default $WORK/selection_inputs)
Usage    python stage_selection_inputs.py --search "$PRIOR_EXPORT" --dest "$WORK/selection_inputs"
Exit status is 1 if any file is missing or fails the md5 check.
"""
import argparse, hashlib, os, re, shutil, sys

HERE = os.path.dirname(os.path.abspath(__file__))


def md5(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--search", nargs="+", required=True)
    ap.add_argument("--dest", default=os.path.join(os.environ.get("WORK", "work"), "selection_inputs"))
    ap.add_argument("--manifest", default=os.path.join(HERE, "selection_inputs_manifest.tsv"))
    ap.add_argument("--copy", action="store_true", help="copy instead of symlink")
    a = ap.parse_args()

    want = {}
    for line in open(a.manifest):
        if line.startswith("file\t"):
            continue
        f, size, m, _ = line.rstrip("\n").split("\t", 3)
        want[f] = (int(size), m)

    cand = {f: [] for f in want}
    pre = re.compile(r"^v[0-9a-f]{8}_")
    for root in a.search:
        for dp, _, fns in os.walk(root):
            for fn in fns:
                base = pre.sub("", fn)
                if base in cand:
                    cand[base].append(os.path.join(dp, fn))

    os.makedirs(a.dest, exist_ok=True)
    bad = 0
    for f, (size, m) in sorted(want.items()):
        hit = None
        for p in cand[f]:
            if os.path.getsize(p) == size and md5(p) == m:
                hit = p
                break
        dst = os.path.join(a.dest, f)
        if hit is None:
            print(f"MISSING  {f}  ({len(cand[f])} candidate(s), none with md5 {m})")
            bad += 1
            continue
        if os.path.lexists(dst):
            os.remove(dst)
        if a.copy:
            shutil.copy2(hit, dst)
        else:
            os.symlink(os.path.abspath(hit), dst)
        print(f"ok       {f}")
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
