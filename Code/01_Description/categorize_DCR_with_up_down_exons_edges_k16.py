
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced architecture classifier (string windows) USING EDGE COMPARISON (k=16).

We refine LW/RW connection categories using:
  - Exon identity (from lowercase flank of window, as before)
  - DCR edge identity via uppercase block edges (prefix/suffix up to 16 bp)

New refined statuses (same labels as before, but edges instead of junctions):
  - Connected-shared                  (same exon + same EDGE across frames)
  - Connected-distinct-same-exon      (same exon + different EDGE)
  - Connected-distinct-diff-exon      (different exon)

Terminals (N-Terminal / C-Terminal) are preserved as-is.

Also emits L/R edge tokens and per-row edge sharing status vs other frames.

INPUT
-----
TSV/CSV with columns:
  Index, Gene, Groups, Frames, LeftWindow, RightWindow, TransAliPctID, LW, RW
Where LeftWindow = <lowercase_exon><UPPERCASE_DCR> and
      RightWindow = <UPPERCASE_DCR><lowercase_exon>.

USAGE
-----
python categorize_DCR_with_up_down_exons_edges_k16.py \
  --input your_table.tsv \
  --out-detailed DCR_classification_detailed.tsv \
  --out-upset DCR_classification_upset.tsv \
  [--sep \\t or ,]

The UpSet file remains based on LW_status_new / RW_status_new.
"""

import argparse
import os
import re
from collections import Counter, defaultdict

import pandas as pd

TERMINALS = {"N-Terminal", "C-Terminal"}
EDGE_K = 16  # default

def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True, help="Input TSV/CSV table with architecture calls")
    ap.add_argument("--sep", default=None, help="Override delimiter (auto by suffix if omitted)")
    ap.add_argument("--out-detailed", default="DCR_classification_detailed.tsv",
                    help="Detailed per-row table with refined statuses and edge tokens")
    ap.add_argument("--out-upset", default="DCR_classification_upset.tsv",
                    help="UpSet-ready wide table collapsed at Index")
    return ap.parse_args()

def detect_sep(path, override):
    if override is not None:
        return override
    ext = os.path.splitext(path)[1].lower()
    if ext in (".tsv", ".txt"):
        return "\t"
    return ","

def is_terminal(val: str) -> bool:
    if val is None:
        return True
    s = str(val)
    if s in TERMINALS:
        return True
    if s.startswith("Insufficient_Nucls"):
        return True
    return False

def exon_token_from_window(window: str, end: str) -> str | None:
    """
    For LW: lowercase run BEFORE first uppercase block.
    For RW: lowercase run AFTER last uppercase block.
    """
    if window is None:
        return None
    s = str(window)
    if is_terminal(s):
        return None
    if s.startswith("Insufficient_Nucls_"):
        return None
    if end == "LW":
        m = re.match(r"^([a-z]+)[A-Z]+", s)
        return m.group(1) if m else None
    else:
        m = re.match(r"[A-Z]+([a-z]+)$", s)
        return m.group(1) if m else None

def edge_token_from_window(window: str, end: str, k: int = EDGE_K) -> str | None:
    """
    Extract edge token from uppercase DCR block:
      - LW: take the FIRST uppercase run, return its PREFIX of length <= k
      - RW: take the LAST  uppercase run, return its SUFFIX of length <= k
    If fewer than k uppercase letters are available, return all available.
    """
    if window is None:
        return None
    s = str(window)
    if is_terminal(s):
        return None
    runs = list(re.finditer(r"[A-Z]+", s))
    if not runs:
        return None
    if end == "LW":
        u = runs[0].group(0)
        n = min(k, len(u))
        return u[:n]
    else:
        u = runs[-1].group(0)
        n = min(k, len(u))
        return u[-n:]

def refined_status_from_exon_and_edge(orig_status, exon_token, edge_token,
                                      other_exons, other_edges):
    """
    Decide refined 'Connected-*' label using exon (same/diff) and edge (shared/distinct).
    Terminals kept as-is.
    If exon/edge token missing or no other frames present, fall back to original status.
    """
    if is_terminal(orig_status):
        return orig_status
    if exon_token is None or edge_token is None:
        return orig_status
    if not other_exons and not other_edges:
        return orig_status

    same_exon = exon_token in other_exons
    same_edge = edge_token in other_edges

    if same_exon and same_edge:
        return "Connected-shared"
    elif same_exon and not same_edge:
        return "Connected-distinct-same-exon"
    else:
        return "Connected-distinct-diff-exon"

def main():
    args = parse_args()
    sep = detect_sep(args.input, args.sep)
    df = pd.read_csv(args.input, sep=sep, dtype=str)

    # Preserve originals for comparison
    df["LW_status_original"] = df["LW"]
    df["RW_status_original"] = df["RW"]

    # Outputs to build
    lw_new, rw_new = [], []
    Ltok_list, Rtok_list = [], []
    Lstat_list, Rstat_list = [], []

    # Process per Index (DCR) to compare across frames
    for idx_val, sub in df.groupby("Index", sort=False):
        # Precompute tokens for the group
        sub = sub.copy()
        sub["lw_exon"] = sub["LeftWindow"].apply(lambda x: exon_token_from_window(x, "LW"))
        sub["rw_exon"] = sub["RightWindow"].apply(lambda x: exon_token_from_window(x, "RW"))
        sub["L_edge_token"] = sub["LeftWindow"].apply(lambda x: edge_token_from_window(x, "LW", EDGE_K))
        sub["R_edge_token"] = sub["RightWindow"].apply(lambda x: edge_token_from_window(x, "RW", EDGE_K))

        for i, row in sub.iterrows():
            frame = row["Frames"]

            # Other frames' tokens (non-terminal only)
            others = sub[sub["Frames"] != frame]
            other_lw_exons = set(t for t in others["lw_exon"] if t is not None)
            other_rw_exons = set(t for t in others["rw_exon"] if t is not None)
            other_L_edges = set(t for t in others["L_edge_token"] if t is not None)
            other_R_edges = set(t for t in others["R_edge_token"] if t is not None)

            # Decide refined LW/RW statuses using exon + EDGE
            lw_refined = refined_status_from_exon_and_edge(
                row["LW"], row["lw_exon"], row["L_edge_token"], other_lw_exons, other_L_edges
            )
            rw_refined = refined_status_from_exon_and_edge(
                row["RW"], row["rw_exon"], row["R_edge_token"], other_rw_exons, other_R_edges
            )

            # Edge-only shared/distinct flags (helpful for QA)
            def edge_stat(tok, otherset):
                if tok is None:
                    return "NA"
                return "shared" if tok in otherset else "distinct"

            Lstat = edge_stat(row["L_edge_token"], other_L_edges)
            Rstat = edge_stat(row["R_edge_token"], other_R_edges)

            # Append
            lw_new.append(lw_refined)
            rw_new.append(rw_refined)
            Ltok_list.append(row["L_edge_token"])
            Rtok_list.append(row["R_edge_token"])
            Lstat_list.append(Lstat)
            Rstat_list.append(Rstat)

    # Attach to df in original row order
    df["LW_status_new"] = lw_new
    df["RW_status_new"] = rw_new
    df["L_edge_token"] = Ltok_list
    df["R_edge_token"] = Rtok_list
    df["L_edge_status"] = Lstat_list
    df["R_edge_status"] = Rstat_list

    # Save detailed
    df.to_csv(args.out_detailed, sep="\t", index=False)

    # Build UpSet-wide table per Index using refined LW/RW statuses
    per_index_patterns = defaultdict(set)
    for _, row in df.iterrows():
        i = row["Index"]
        per_index_patterns[i].add(f"LW:{row['LW_status_new']}")
        per_index_patterns[i].add(f"RW:{row['RW_status_new']}")

    all_cols = sorted({p for pats in per_index_patterns.values() for p in pats})
    freq = Counter()
    for pats in per_index_patterns.values():
        for p in pats:
            freq[p] += 1
    all_cols_sorted = sorted(all_cols, key=lambda c: (-freq[c], c))

    upset_rows = []
    for index_val, pats in per_index_patterns.items():
        rowdict = {"Index": index_val, "PatternString": ",".join(sorted(pats))}
        for c in all_cols_sorted:
            rowdict[c] = 1 if c in pats else 0
        upset_rows.append(rowdict)

    upset_df = pd.DataFrame(upset_rows)
    upset_df = upset_df[["Index", "PatternString"] + all_cols_sorted]
    upset_df.to_csv(args.out_upset, sep="\t", index=False)

    print(f"Wrote detailed table: {args.out_detailed}")
    print(f"Wrote UpSet table:    {args.out_upset}")
    print(f"Edge k-mer used: {EDGE_K} | Pattern columns: {', '.join(all_cols_sorted)}")

if __name__ == "__main__":
    main()
