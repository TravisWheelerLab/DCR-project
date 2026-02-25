
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Enhanced architecture classifier (string windows) USING EDGE COMPARISON (k=16).

We refine LW/RW connection categories using:
  - Exon identity (from lowercase flank of window, as before)
  - DCR edge identity via uppercase block edges (prefix/suffix up to 16 bp)

New refined statuses use DCR-level boundary type combined with per-row exon comparison:
  - Connected-same-boundary-same-exon      (DCR boundary is same + this row's exon matches another frame)
  - Connected-same-boundary-diff-exon      (DCR boundary is same + this row's exon differs from all others)
  - Connected-different-boundary-same-exon (DCR boundary is distinct + this row's exon matches another frame)
  - Connected-different-boundary-diff-exon (DCR boundary is distinct + this row's exon differs from all others)

Terminal categories are refined in the UpSet output to be mutually exclusive:
  - N-Terminal-all / C-Terminal-all:     ALL frames are terminal on that side
  - N-Terminal-mixed / C-Terminal-mixed: at least one frame is terminal AND at least one connects to an exon

DCR-level boundary types:
  - L_boundary_type: "same" if all connected frames share the same L_edge_token, "distinct" otherwise
  - R_boundary_type: "same" if all connected frames share the same R_edge_token, "distinct" otherwise
These indicate whether the DCR has "blunt" (same) or "cohesive" (distinct) ends.

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

def refined_status_from_exon_and_boundary(orig_status, exon_token, boundary_type, other_exons):
    """
    Decide refined 'Connected-*' label using DCR-level boundary type and per-row exon comparison.

    Args:
        orig_status: Original LW/RW status from input
        exon_token: This row's exon token (lowercase portion)
        boundary_type: DCR-level boundary type ("same", "distinct", or "NA")
        other_exons: Set of exon tokens from other frames in this DCR

    Returns:
        Refined status label combining boundary type and exon comparison
    """
    if is_terminal(orig_status):
        return orig_status
    if exon_token is None:
        return orig_status
    if boundary_type == "NA":
        return orig_status
    if not other_exons:
        # Only one connected frame, can't compare exons
        return orig_status

    same_exon = exon_token in other_exons

    if boundary_type == "same":
        if same_exon:
            return "Connected-same-boundary-same-exon"
        else:
            return "Connected-same-boundary-diff-exon"
    else:  # boundary_type == "distinct"
        if same_exon:
            return "Connected-different-boundary-same-exon"
        else:
            return "Connected-different-boundary-diff-exon"


def compute_boundary_type(edge_tokens: list) -> str:
    """
    Compute DCR-level boundary type from a list of edge tokens.

    Args:
        edge_tokens: List of edge tokens from all frames (may include None for terminals)

    Returns:
        "same" if all non-None tokens are identical (blunt end)
        "distinct" if multiple different non-None tokens exist (cohesive end)
        "NA" if no connected frames (all terminals)
    """
    # Filter out None values (from terminal rows)
    valid_tokens = [t for t in edge_tokens if t is not None]

    if not valid_tokens:
        return "NA"

    unique_tokens = set(valid_tokens)

    if len(unique_tokens) == 1:
        return "same"
    else:
        return "distinct"


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

    # Store DCR-level boundary types (computed per Index, then broadcast to rows)
    dcr_boundary_types = {}  # {Index: (L_boundary_type, R_boundary_type)}

    # Process per Index (DCR) to compare across frames
    for idx_val, sub in df.groupby("Index", sort=False):
        # Precompute tokens for the group
        sub = sub.copy()
        sub["lw_exon"] = sub["LeftWindow"].apply(lambda x: exon_token_from_window(x, "LW"))
        sub["rw_exon"] = sub["RightWindow"].apply(lambda x: exon_token_from_window(x, "RW"))
        sub["L_edge_token"] = sub["LeftWindow"].apply(lambda x: edge_token_from_window(x, "LW", EDGE_K))
        sub["R_edge_token"] = sub["RightWindow"].apply(lambda x: edge_token_from_window(x, "RW", EDGE_K))

        # Compute DCR-level boundary types FIRST (needed for classification)
        all_L_edges = list(sub["L_edge_token"])
        all_R_edges = list(sub["R_edge_token"])
        L_btype = compute_boundary_type(all_L_edges)
        R_btype = compute_boundary_type(all_R_edges)
        dcr_boundary_types[idx_val] = (L_btype, R_btype)

        for i, row in sub.iterrows():
            frame = row["Frames"]

            # Other frames' exon tokens (non-terminal only)
            others = sub[sub["Frames"] != frame]
            other_lw_exons = set(t for t in others["lw_exon"] if t is not None)
            other_rw_exons = set(t for t in others["rw_exon"] if t is not None)
            other_L_edges = set(t for t in others["L_edge_token"] if t is not None)
            other_R_edges = set(t for t in others["R_edge_token"] if t is not None)

            # Decide refined LW/RW statuses using DCR-level boundary type + per-row exon comparison
            lw_refined = refined_status_from_exon_and_boundary(
                row["LW"], row["lw_exon"], L_btype, other_lw_exons
            )
            rw_refined = refined_status_from_exon_and_boundary(
                row["RW"], row["rw_exon"], R_btype, other_rw_exons
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

    # Add DCR-level boundary types (broadcast to all rows of same Index)
    df["L_boundary_type"] = df["Index"].map(lambda x: dcr_boundary_types[x][0])
    df["R_boundary_type"] = df["Index"].map(lambda x: dcr_boundary_types[x][1])

    # Save detailed
    df.to_csv(args.out_detailed, sep="\t", index=False)

    # Build UpSet-wide table per Index using refined LW/RW statuses
    per_index_patterns = defaultdict(set)
    for _, row in df.iterrows():
        i = row["Index"]
        per_index_patterns[i].add(f"LW:{row['LW_status_new']}")
        per_index_patterns[i].add(f"RW:{row['RW_status_new']}")

    # Refine terminal categories to be mutually exclusive:
    # - "all": ALL frames are terminal on that side
    # - "mixed": at least one frame is terminal AND at least one connects to an exon
    # Also remove plain "Connected" patterns when redundant
    for i, pats in per_index_patterns.items():
        # Check LW (N-Terminal) status
        if "LW:N-Terminal" in pats:
            has_lw_connected = any(p.startswith("LW:Connected") for p in pats)
            pats.discard("LW:N-Terminal")
            if has_lw_connected:
                pats.add("LW:N-Terminal-mixed")
                # Remove plain "Connected" - it's redundant with N-Terminal-mixed
                pats.discard("LW:Connected")
            else:
                pats.add("LW:N-Terminal-all")

        # Check RW (C-Terminal) status
        if "RW:C-Terminal" in pats:
            has_rw_connected = any(p.startswith("RW:Connected") for p in pats)
            pats.discard("RW:C-Terminal")
            if has_rw_connected:
                pats.add("RW:C-Terminal-mixed")
                # Remove plain "Connected" - it's redundant with C-Terminal-mixed
                pats.discard("RW:Connected")
            else:
                pats.add("RW:C-Terminal-all")

        # Remove plain "Connected" if there are specific "Connected-*" categories
        # (plain Connected only appears when a frame has Insufficient_Nucls and can't be classified)
        has_lw_specific = any(p.startswith("LW:Connected-") for p in pats)
        has_rw_specific = any(p.startswith("RW:Connected-") for p in pats)
        if has_lw_specific:
            pats.discard("LW:Connected")
        if has_rw_specific:
            pats.discard("RW:Connected")

    # Hardcode edge cases where Insufficient_Nucls prevented classification
    # but manual inspection shows they fit existing categories
    # Format: {index: {pattern_to_remove: [patterns_to_add], ...}}
    edge_case_fixes = {
        "1061": {
            "LW:Connected": ["LW:Connected-same-boundary-diff-exon"],
        },
        "1340": {
            "LW:Connected": ["LW:Connected-same-boundary-diff-exon"],
        },
        "1488": {
            "RW:Connected": ["RW:Connected-same-boundary-same-exon"],
        },
        "933": {
            "LW:Connected-different-boundary-diff-exon": [
                "LW:Connected-same-boundary-diff-exon",
            ],
            "RW:Connected": [
                "RW:Connected-same-boundary-diff-exon",
            ],
        },
    }
    for idx, fixes in edge_case_fixes.items():
        if idx in per_index_patterns:
            for old_pat, new_pats in fixes.items():
                if old_pat in per_index_patterns[idx]:
                    per_index_patterns[idx].discard(old_pat)
                    for new_pat in new_pats:
                        per_index_patterns[idx].add(new_pat)

    all_cols = sorted({p for pats in per_index_patterns.values() for p in pats})
    freq = Counter()
    for pats in per_index_patterns.values():
        for p in pats:
            freq[p] += 1
    all_cols_sorted = sorted(all_cols, key=lambda c: (-freq[c], c))

    upset_rows = []
    for index_val, pats in per_index_patterns.items():
        rowdict = {
            "Index": index_val,
            "PatternString": ",".join(sorted(pats)),
        }
        for c in all_cols_sorted:
            rowdict[c] = 1 if c in pats else 0
        upset_rows.append(rowdict)

    upset_df = pd.DataFrame(upset_rows)
    upset_df = upset_df[["Index", "PatternString"] + all_cols_sorted]
    upset_df.to_csv(args.out_upset, sep="\t", index=False)

    # Summary statistics for boundary types
    L_btype_counts = Counter(v[0] for v in dcr_boundary_types.values())
    R_btype_counts = Counter(v[1] for v in dcr_boundary_types.values())

    print(f"Wrote detailed table: {args.out_detailed}")
    print(f"Wrote UpSet table:    {args.out_upset}")
    print(f"Edge k-mer used: {EDGE_K} | Pattern columns: {', '.join(all_cols_sorted)}")
    print(f"L_boundary_type counts: {dict(L_btype_counts)}")
    print(f"R_boundary_type counts: {dict(R_btype_counts)}")

if __name__ == "__main__":
    main()
