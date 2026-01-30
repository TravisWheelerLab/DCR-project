#!/usr/local/bin/python 
import pandas as pd
from pybedtools import BedTool as bt
import os
import argparse
import numpy as np

parser = argparse.ArgumentParser(description="Categorize DCR with up/down exons")
parser.add_argument("-i", "--input", required=True, help="Input TSV file")
args = parser.parse_args()

intable    = pd.read_csv(args.input, sep = "\t")
DCRindices = pd.unique(intable['Index'])

def find_end_patterns():
    endPats = []

    for DCR in DCRindices:
        DCRtable = intable[intable['Index'] == DCR]
        #print(DCRtable)
        frames = pd.unique(DCRtable['Frames'])

        # Build mappings: sequence -> set of frames it appears in
        lw_frame_map = DCRtable.groupby('LeftWindow')['Frames'].apply(set).to_dict()
        rw_frame_map = DCRtable.groupby('RightWindow')['Frames'].apply(set).to_dict()

        for frame in frames:
            subtable = DCRtable[DCRtable['Frames'] == frame]

            frameLW = []
            frameRW = []

            for _, row in subtable.iterrows():
                # Process LW
                if row['LW'] == 'Connected':
                    seq = row['LeftWindow']
                    if len(lw_frame_map.get(seq, set())) > 1:
                        frameLW.append('Connected_shared')
                    else:
                        frameLW.append('Connected_distinct')
                else:
                    frameLW.append(row['LW'])

                # Process RW
                if row['RW'] == 'Connected':
                    seq = row['RightWindow']
                    if len(rw_frame_map.get(seq, set())) > 1:
                        frameRW.append('Connected_shared')
                    else:
                        frameRW.append('Connected_distinct')
                else:
                    frameRW.append(row['RW'])

            frameLW = list(pd.unique(frameLW))
            frameRW = list(pd.unique(frameRW))

            endpattern = [frameLW, frameRW]
            endpattern_tuple = (tuple(frameLW), tuple(frameRW))

            if endpattern_tuple not in [tuple(map(tuple, ep)) for ep in endPats]:
                endPats.append(endpattern)

    return endPats


def count_end_patterns():
    '''
    Counts the frequency of each unique end pattern in endPats 
    for each DCR and returns a DataFrame.
    '''
    pattern_counts = []

    for DCR in DCRindices:
        DCRtable = intable[intable['Index'] == DCR]
        frames = pd.unique(DCRtable['Frames'])

        # Build mappings: sequence -> set of frames it appears in
        lw_frame_map = DCRtable.groupby('LeftWindow')['Frames'].apply(set).to_dict()
        rw_frame_map = DCRtable.groupby('RightWindow')['Frames'].apply(set).to_dict()

        # Initialize a count list with zeros, one for each pattern in endPats
        counts = [0] * len(endPats)

        for frame in frames:
            subtable = DCRtable[DCRtable['Frames'] == frame]

            frameLW = []
            frameRW = []

            for _, row in subtable.iterrows():
                # LW
                if row['LW'] == 'Connected':
                    seq = row['LeftWindow']
                    if len(lw_frame_map.get(seq, set())) > 1:
                        frameLW.append('Connected_shared')
                    else:
                        frameLW.append('Connected_distinct')
                else:
                    frameLW.append(row['LW'])

                # RW
                if row['RW'] == 'Connected':
                    seq = row['RightWindow']
                    if len(rw_frame_map.get(seq, set())) > 1:
                        frameRW.append('Connected_shared')
                    else:
                        frameRW.append('Connected_distinct')
                else:
                    frameRW.append(row['RW'])

            # Remove duplicates from frame-level list
            frameLW = list(pd.unique(frameLW))
            frameRW = list(pd.unique(frameRW))

            # Match against known patterns
            for i, pat in enumerate(endPats):
                if (tuple(frameLW), tuple(frameRW)) == (tuple(pat[0]), tuple(pat[1])):
                    counts[i] += 1
                    break

        pattern_counts.append([DCR] + counts)

    # Create DataFrame with DCRs as rows and patterns as columns
    columns = ["DCR"] + pattern_labels
    df_counts = pd.DataFrame(pattern_counts, columns=columns)

    return df_counts

# Get list of unique endind patterns
endPats = find_end_patterns()
# Generate unique labels for each pattern
pattern_labels = [f"P{i}" for i in range(len(endPats))]  # P0, P1, P2, ...
pattern_dict = dict(zip(pattern_labels, endPats))  # Dictionary for reference
# Compute counts
df_result = count_end_patterns()
# Convert pattern correspondence table to DataFrame
df_patterns = pd.DataFrame(list(pattern_dict.items()), columns=["Pattern_Label", "Pattern"])

# Clean and standardize before value_counts
df_patterns_only = df_result.drop(columns=["DCR"]).copy()

# Coerce all to int — even if they already are
df_patterns_only = df_patterns_only.astype("int64")

# Count unique pattern rows
df_unique_series = df_patterns_only.value_counts().reset_index()
df_unique_series.columns = pattern_labels + ["DCR_Count"]

# Sanity check
total_counted = df_unique_series["DCR_Count"].sum()
expected_total = df_result.shape[0]
print(f"✅ Total DCRs grouped: {total_counted} / {expected_total}")

# # Group DCRs by unique pattern distributions and count occurrences
# df_unique_series = df_result.drop(columns=["DCR"]).value_counts().reset_index()
# df_unique_series.columns = pattern_labels + ["DCR_Count"]

# Save results
df_result.to_csv("output/01_Description/DCR_pattern_matrix.csv", index=False)  # Matrix file
df_patterns.to_csv("output/01_Description/pattern_correspondence.csv", index=False)  # Correspondence table
df_unique_series.to_csv("output/01_Description/unique_series_counts.csv", index=False)  # Unique series counts

# Print results
print("\nDCR Pattern Matrix:")
print(df_result)
print("\nPattern Correspondence Table:")
print(df_patterns)
print("\nUnique Series Count Table:")
print(df_unique_series)