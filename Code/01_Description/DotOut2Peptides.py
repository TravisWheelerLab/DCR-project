#!/usr/bin/env python3
"""
Make an isoform-level multi-FASTA from a directory of DCR '.out' reports.

- Parses both 'Group N' and 'Groups 1,2,3' formats.
- Maps each Group to its isoform IDs from lines like:
    > Group 1: homo_sapiens|GENE/ID|1
               homo_sapiens|GENE/ID|1.1
- Extracts the DCR peptide for each Group and replicates it per isoform in that Group.
- Adds the Group's frame to the header when available.
- Filters by identity from:
  "Percents identity for alignments of translated dual-coding nucleotides to group aminos"

Output header format:
    ><DCRID>_iso<ISOID>[_f<frame>]
  Example: >203_iso1.1_f3

Usage
-----
python dcr_out_dir_to_iso_multifasta.py INDIR -o all_DCR_iso.fasta
python dcr_out_dir_to_iso_multifasta.py INDIR -o all.fasta --recursive --minid 85
python dcr_out_dir_to_iso_multifasta.py INDIR -o all.fasta --drop-missing
"""

import argparse, glob, re, sys
from pathlib import Path
from typing import Dict, List, Tuple, Optional

IDENTITY_HDR = "Percents identity for alignments of translated dual-coding nucleotides to group aminos"
OVERLAID_HDR = "Overlaid alignment of the dual-coding regions for each group"
HEADER_MAXLEN = 50

def parse_dcr_id(text: str, fallback: str) -> str:
    m = re.search(r"DCE\s+Index\s*:\s*(\d+)", text)
    return m.group(1).strip() if m else fallback

def parse_group_identities_and_frames(text: str) -> Dict[str, Tuple[float, Optional[int]]]:
    """ Return {gid: (percent_identity, frame)}; handles 'Group' and 'Groups' lines. """
    out: Dict[str, Tuple[float, Optional[int]]] = {}
    start = text.find(IDENTITY_HDR)
    if start == -1:
        return out
    end = text.find(OVERLAID_HDR, start)
    block = text[start:end] if end != -1 else text[start:]
    pat = re.compile(r"Group[s]?\s+([0-9, ]+)\s*:=\s*([0-9]+(?:\.[0-9]+)?)\s*%(?:\s*\[frame:(\d)\])?", re.I)
    for ln in block.splitlines():
        m = pat.search(ln)
        if not m:
            continue
        ids_str, pct_str, frame_str = m.groups()
        pct = float(pct_str)
        frame = int(frame_str) if frame_str else None
        ids = [s.strip() for s in ids_str.replace(" ", "").split(",") if s.strip()]
        for gid in ids:
            out[gid] = (pct, frame)
    return out

def parse_overlaid_group_sequences(text: str) -> Dict[str, str]:
    """ Parse DCR AA sequences per Group from the 'Overlaid alignment' block. """
    start = text.find(OVERLAID_HDR)
    block = text[start:] if start != -1 else text
    seqs: Dict[str, List[str]] = {}
    line_pat = re.compile(r"^\s*Group[s]?\s+([0-9, ]+)\s+(.+)$")
    for ln in block.splitlines():
        m = line_pat.match(ln)
        if not m:
            continue
        tail = m.group(2)
        if ":=" in tail:  # skip identity echoes
            continue
        aa = "".join(c for c in tail if c.isalpha())
        if not aa:
            continue
        ids_str = m.group(1)
        ids = [s.strip() for s in ids_str.replace(" ", "").split(",") if s.strip()]
        for gid in ids:
            seqs.setdefault(gid, []).append(aa)
    return {gid: "".join(parts) for gid, parts in seqs.items() if "".join(parts)}

def parse_group_to_isoforms(text: str) -> Dict[str, List[str]]:
    """
    Map Group N -> [isoform IDs], accepting decimals like '1.1'.
    Example:
    > Group 1: homo_sapiens|gene/id|1
               homo_sapiens|gene/id|14
    """
    lines = text.splitlines()
    group2isos: Dict[str, List[str]] = {}
    i = 0
    while i < len(lines):
        line = lines[i]
        m = re.match(r"^\s*>\s*Group\s+(\d+)\s*:\s*(.+)$", line)
        if not m:
            i += 1
            continue
        gid = m.group(1)
        payload = m.group(2)
        def iso_from_payload(s: str) -> Optional[str]:
            return s.strip().split("|")[-1] if "|" in s else None
        isos = []
        first = iso_from_payload(payload)
        if first:
            isos.append(first.strip())
        j = i + 1
        while j < len(lines):
            nxt = lines[j]
            if nxt.strip() == "":
                break
            if nxt.lstrip().startswith(("Genome Ranges:", "Protein Seq.", "Coding Nucl.", "Translation",
                                        "Chromosome", "Species", "> Group", "Group ")):
                break
            more = iso_from_payload(nxt)
            if more:
                isos.append(more.strip())
                j += 1
            else:
                break
        group2isos[gid] = isos
        i = j
    return group2isos

def wrap(seq: str, width: int = 60) -> str:
    return "\n".join(seq[i:i+width] for i in range(0, len(seq), width)) or ""

def unique_and_trim_header(header: str, seen: Dict[str, int]) -> str:
    def enforce(h: str) -> str:
        return h if len(h) <= HEADER_MAXLEN else h[:HEADER_MAXLEN]
    h = enforce(header)
    if h not in seen:
        seen[h] = 1
        return h
    while True:
        seen[h] += 1
        cand = enforce(f"{header}#{seen[h]}")
        if cand not in seen:
            seen[cand] = 1
            return cand

def main():
    import argparse
    ap = argparse.ArgumentParser(description="Isoform-level multi-FASTA from DCR '.out' files.")
    ap.add_argument("indir", type=Path, help="Directory with .out files")
    ap.add_argument("-o", "--outfile", type=Path, required=True, help="Output multi-FASTA path")
    ap.add_argument("--pattern", default="*.out", help="Glob pattern (default: *.out)")
    ap.add_argument("--recursive", action="store_true", help="Recurse into subdirectories")
    ap.add_argument("--minid", type=float, default=80.0, help="Minimum identity to keep (default: 80.0)")
    ap.add_argument("--drop-missing", action="store_true", help="Drop groups without identity entry")
    ap.add_argument("--no-frame", action="store_true", help="Do NOT append frame to headers")
    ap.add_argument("--wrap", type=int, default=60, help="FASTA wrap width (default: 60)")
    args = ap.parse_args()

    files = ([Path(p) for p in glob.glob(str(args.indir / "**" / args.pattern), recursive=True)]
             if args.recursive else sorted(args.indir.glob(args.pattern)))
    files = [f for f in files if f.is_file()]
    if not files:
        print(f"[ERROR] No input files in {args.indir} matching {args.pattern}", file=sys.stderr); sys.exit(2)

    entries = []
    for f in files:
        text = f.read_text(errors="ignore")
        dcr = parse_dcr_id(text, f.stem)
        id_frames = parse_group_identities_and_frames(text)
        seqs = parse_overlaid_group_sequences(text)
        g2i = parse_group_to_isoforms(text)

        for gid, seq in seqs.items():
            # identity & frame
            if gid in id_frames:
                pct, frame = id_frames[gid]
                if pct < args.minid:
                    continue
            else:
                frame = None
                if args.drop_missing:
                    continue
            # expand to isoforms; if none recorded for the group, skip (isoform-only output)
            isos = g2i.get(gid, [])
            if not isos:
                continue
            for iso in isos:
                iso_clean = re.sub(r"[^A-Za-z0-9._-]", "_", iso)  # allow decimals like 1.1
                header = f"{dcr}_iso{iso_clean}"
                if not args.no_frame and frame is not None:
                    header += f"_f{frame}"
                entries.append((header, seq))

    if not entries:
        args.outfile.write_text("")
        print(f"[INFO] No entries passed filters; wrote empty file: {args.outfile}", file=sys.stderr)
        sys.exit(0)

    args.outfile.parent.mkdir(parents=True, exist_ok=True)
    seen = {}
    with args.outfile.open("w") as out:
        for header, seq in entries:
            out.write(f">{unique_and_trim_header(header, seen)}\n{wrap(seq, args.wrap)}\n")
    print(f"[OK] Wrote {len(entries)} sequences from {len(files)} files to {args.outfile}", file=sys.stderr)

if __name__ == "__main__":
    main()
