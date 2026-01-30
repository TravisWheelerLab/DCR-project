import pandas as pd

# Read the table (tabs or whitespace-delimited)
# If you already have a DataFrame `df`, skip this read_csv.
df = pd.read_csv(
    "output/01_Description/isoforms_to_DCR_metadata.tsv",
    sep=r"\s+",
    header=None,
    names=["isoform", "frame", "ignore", "DCR"],
    engine="python"
)

# Extract isoform number after the last '|'
df["isoform_num"] = df["isoform"].str.extract(r'\|([^|]+)$')
df["frame"] = df["frame"].astype(int)
df["DCR"] = df["DCR"].astype(int)

# Canonical frame(s) per DCR: frames used by isoform exactly "1"
canon = (
    df[df["isoform_num"] == "1"]
    .groupby("DCR")["frame"]
    .agg(lambda x: set(x))
    .rename("canonical_frames")
    .reset_index()
)

# DCRs with no isoform 1 present
dcr_no_iso1 = sorted(set(df["DCR"]) - set(canon["DCR"]))

# Attach canonical frame sets (NaN for DCRs lacking iso1)
out = df.merge(canon, on="DCR", how="left")

def label_row(r):
    if pd.isna(r["canonical_frames"]):
        return "undetermined"
    return "canonical" if r["frame"] in r["canonical_frames"] else "alternative"

out["frame_status"] = out.apply(label_row, axis=1)

# Make a clean string column without braces; "NA" when undetermined
out["canonical_frame"] = out["canonical_frames"].apply(
    lambda x: "NA" if pd.isna(x) else ",".join(str(i) for i in sorted(x))
)

# Keep exactly 7 columns, with canonical_frame as plain text
out = out[["isoform", "frame", "ignore", "DCR", "isoform_num", "canonical_frame", "frame_status"]]

print(out.to_string(index=False))
print("DCR with no isoform 1:", dcr_no_iso1)

# Optional: write to TSV
out.to_csv("output/01_Description/isoforms_to_DCR_metadata_canon.tsv", sep="\t", index=False)
# out["frame_status"] = out.apply(label_row, axis=1)

# Write output DataFrame to TSV
# out.to_csv("output/01_Description/isoforms_to_DCR_metadata_canon.tsv", sep="\t", index=False)
# print(out)
# print("DCR with no isoform 1:", dcr_no_iso1)
