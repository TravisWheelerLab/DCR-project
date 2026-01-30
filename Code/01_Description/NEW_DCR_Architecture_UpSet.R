library(ComplexUpset)
library(dplyr)
library(ggplot2)

# Read table
df <- readr::read_tsv("/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/output/01_Description/NEW_DCR_architectures_exon_edges_upset.tsv", show_col_types = FALSE)

# Identify the binary pattern columns (everything except Index and PatternString)
pattern_cols <- setdiff(colnames(df), c("Index", "PatternString"))

# Ensure numeric -> logical
df <- df |>
  mutate(across(all_of(pattern_cols), ~ .x %in% c(1, "1", TRUE)))

# Rename for the paper:
# - LW -> Upstream, RW -> Downstream
# - same-boundary -> same-DCR-boundary
# - different-boundary -> different-DCR-boundary
# - Connected$ (alone, rare edge case) -> Connected-unclassified
names(df) <- names(df) |>
  gsub(pattern = "LW:", replacement = "Upstream:") |>
  gsub(pattern = "RW:", replacement = "Downstream:") |>
  gsub(pattern = "Connected-same-boundary", replacement = "Connected-same-DCR-boundary") |>
  gsub(pattern = "Connected-different-boundary", replacement = "Connected-different-DCR-boundary") |>
  gsub(pattern = ":Connected$", replacement = ":Connected-unclassified")

pattern_cols <- pattern_cols |>
  gsub(pattern = "LW:", replacement = "Upstream:") |>
  gsub(pattern = "RW:", replacement = "Downstream:") |>
  gsub(pattern = "Connected-same-boundary", replacement = "Connected-same-DCR-boundary") |>
  gsub(pattern = "Connected-different-boundary", replacement = "Connected-different-DCR-boundary") |>
  gsub(pattern = ":Connected$", replacement = ":Connected-unclassified")

# Remove all-zero columns (just in case)
keep_cols <- pattern_cols[colSums(df[, pattern_cols]) > 0]
if (length(keep_cols) == 0) stop("No active pattern columns found in the input.")
df <- df[, c("Index", "PatternString", keep_cols), drop = FALSE]

# Save architecture table
DCRarc <- df
write.table(df, "/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/output/01_Description/DCR_architecture.tsv", quote = FALSE, col.names = TRUE, row.names = FALSE, sep = "\t")

# Sort patterns by total frequency
freq <- colSums(df[, keep_cols])
keep_cols <- keep_cols[order(freq, decreasing = TRUE)]

# Print the column order for reference
message("Pattern columns (sorted by frequency):")
print(keep_cols)

# Build the UpSet plot
# Group columns: Upstream first, then Downstream
# Within each group: ordered by decreasing frequency

upstream_cols <- keep_cols[grepl("^Upstream:", keep_cols)]
downstream_cols <- keep_cols[grepl("^Downstream:", keep_cols)]

# Order each group by frequency (decreasing) - keep_cols is already sorted by freq
# so we just need to preserve that order within each group
upstream_ordered <- upstream_cols[order(freq[upstream_cols], decreasing = TRUE)]
downstream_ordered <- downstream_cols[order(freq[downstream_cols], decreasing = TRUE)]

ordered_cols <- c(upstream_ordered, downstream_ordered)

p <- ComplexUpset::upset(
  df,
  rev(ordered_cols),
  name = "DCR Architecture",
  width_ratio = 0.2,
  sort_sets = FALSE,
  stripes = "white",
  min_size = 1,
  base_annotations = list(
    "Intersection size" = intersection_size(
      text = list(size = 2.5)
    )
  )
)

p <- p & theme(panel.grid.major.x = element_blank())

# Display the plot
print(p)

# Save plot as PDF
output_file <- "/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/output/01_Description/DCR_architecture_upset.pdf"
ggplot2::ggsave(output_file, p, width = 14, height = 8, device = cairo_pdf, bg = "white")
message("Saved UpSet plot to: ", output_file)

# Also save as PNG
output_png <- "/Users/clementgoubert/Library/CloudStorage/GoogleDrive-goubert.clement@gmail.com/My Drive/Wheeler_Lab/2025/06-17_DCR_pipeline/DCR-pipeline/output/01_Description/DCR_architecture_upset.png"
ggplot2::ggsave(output_png, p, width = 14, height = 8, dpi = 300, bg = "white")
message("Saved UpSet plot to: ", output_png)

# Print frequency table summary
freq_tbl <- tibble::tibble(pattern = keep_cols, count = as.integer(freq[keep_cols]))
print(freq_tbl)
