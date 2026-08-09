#!/usr/bin/env Rscript

# Create iTOL color-strip annotation files from metadata for nifH trees.
# Annotation fields:
#   1. sample type
#   2. site
#   3. family
#   4. native status
#   5. tribe
#
# The script is local/RStudio friendly. It does not require extra packages.

project_dir <- "/Users/rosa/Desktop/ALLWork/Madison/Project/ryan-nitfix/intbio-nitfix-git/symbiosis_sorted_all_sample_consensus_sequences_50percent"

metadata_file <- file.path(project_dir, "result", "tables", "intbio_metadata_draft4.csv")
tree_dir <- file.path(project_dir, "result", "trees")
out_dir <- file.path(project_dir, "result", "annotations")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

tree_files <- c(
  strict_pct80_Nle20 = file.path(tree_dir, "nifH_consensus50_strict_single_dominant.treefile"),
  iupac_pct80_Nle20 = file.path(tree_dir, "nifH_consensus50_iupac_all_pass.treefile"),
  strict_pct60_Nle40 = file.path(tree_dir, "nifH_consensus50_strict_single_dominant_pct60_depth10_Nle40.treefile")
)

annotation_fields <- list(
  sample_type = "Type",
  site = "Site",
  family = "Family",
  native = "native",
  tribe = "Tribe"
)

field_labels <- c(
  sample_type = "Sample type",
  site = "Site",
  family = "Host family",
  native = "Native status",
  tribe = "Host tribe"
)

clean_value <- function(x) {
  x <- trimws(as.character(x))
  x[x == "" | is.na(x) | x %in% c("#N/A", "NA", "N/A", "na", "n/a")] <- "unknown"
  x
}

sample_from_tip <- function(tip) {
  sub("\\|.*$", "", tip)
}

normalize_sample_id <- function(x) {
  # Keep exact sample IDs as the main matching rule, but normalize known OAES
  # underscore/hyphen formatting so OAES_19-2-Rh and OAES-19-2-Rh can match.
  x <- trimws(as.character(x))
  x <- sub("^OAES_19", "OAES-19", x)
  x <- sub("^OAES_24", "OAES-24", x)
  x
}

read_newick_tips <- function(tree_file) {
  txt <- paste(readLines(tree_file, warn = FALSE), collapse = "")
  txt <- gsub("[\n\r\t ]+", "", txt)
  # IQ-TREE writes these labels unquoted, for example:
  # CLBJ-14-1-No|nifH|consensus50|strict_single_dominant:0.000002
  # Extract only real leaf labels containing the expected nifH marker.
  # This avoids accidentally keeping Newick parentheses or internal node
  # support labels as part of the iTOL IDs.
  m <- gregexpr("[^(),:;]+\\|nifH\\|consensus50\\|[^(),:;]+", txt, perl = TRUE)
  tips <- regmatches(txt, m)[[1]]
  unique(tips)
}

make_colors <- function(values) {
  values <- sort(unique(values))
  n <- length(values)
  if (n <= 8) {
    base_cols <- c(
      "#2B8CBE", "#F39C12", "#7A5195", "#E74C3C",
      "#2CA25F", "#6C757D", "#D81B60", "#8D6E63"
    )
    cols <- base_cols[seq_len(n)]
  } else {
    cols <- grDevices::hcl.colors(n, palette = "Dark 3")
  }
  names(cols) <- values
  cols
}

write_itol_colorstrip <- function(dat, field_key, value_col, tree_label) {
  values <- clean_value(dat[[value_col]])
  colors <- make_colors(values)
  out_file <- file.path(out_dir, paste0(tree_label, "_itol_", field_key, "_colorstrip.txt"))

  con <- file(out_file, open = "wt")
  on.exit(close(con), add = TRUE)

  writeLines("DATASET_COLORSTRIP", con)
  writeLines("SEPARATOR TAB", con)
  writeLines(paste("DATASET_LABEL", field_labels[[field_key]], sep = "\t"), con)
  writeLines("COLOR\t#000000", con)
  writeLines("STRIP_WIDTH\t35", con)
  writeLines("MARGIN\t5", con)
  writeLines("SHOW_INTERNAL\t0", con)
  writeLines("DATA", con)

  for (i in seq_len(nrow(dat))) {
    tip <- dat$tree_tip[i]
    value <- values[i]
    writeLines(paste(tip, colors[[value]], value, sep = "\t"), con)
  }

  out_file
}

metadata <- read.csv(metadata_file, stringsAsFactors = FALSE, check.names = FALSE)
metadata$metadata_sample <- metadata[["Sample ID"]]
metadata$sample_norm <- normalize_sample_id(metadata$metadata_sample)

# If the metadata contains duplicate normalized sample IDs, keep the first row
# for annotation. This avoids many-to-many joins in tree-tip annotation.
metadata_unique <- metadata[!duplicated(metadata$sample_norm), ]

all_summary <- data.frame()

for (tree_label in names(tree_files)) {
  tree_file <- tree_files[[tree_label]]
  if (!file.exists(tree_file)) {
    message("Skipping missing tree: ", tree_file)
    next
  }

  tips <- read_newick_tips(tree_file)
  tip_df <- data.frame(
    tree_tip = tips,
    sample = sample_from_tip(tips),
    sample_norm = normalize_sample_id(sample_from_tip(tips)),
    stringsAsFactors = FALSE
  )

  annotated <- merge(
    tip_df,
    metadata_unique,
    by = "sample_norm",
    all.x = TRUE,
    sort = FALSE
  )

  # Restore tree order after merge.
  annotated <- annotated[match(tip_df$tree_tip, annotated$tree_tip), ]

  for (field_key in names(annotation_fields)) {
    value_col <- annotation_fields[[field_key]]
    if (!value_col %in% names(annotated)) {
      warning("Missing metadata column: ", value_col)
      next
    }
    write_itol_colorstrip(annotated, field_key, value_col, tree_label)
  }

  summary_one <- data.frame(
    tree = tree_label,
    tree_file = tree_file,
    tips = nrow(tip_df),
    tips_with_metadata = sum(!is.na(annotated$metadata_sample)),
    tips_missing_metadata = sum(is.na(annotated$metadata_sample)),
    stringsAsFactors = FALSE
  )
  all_summary <- rbind(all_summary, summary_one)

  missing_file <- file.path(out_dir, paste0(tree_label, "_tips_missing_metadata.tsv"))
  missing <- annotated[is.na(annotated$metadata_sample), c("tree_tip", "sample")]
  write.table(missing, missing_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

summary_file <- file.path(out_dir, "metadata_itol_annotation_summary.tsv")
write.table(all_summary, summary_file, sep = "\t", quote = FALSE, row.names = FALSE)

message("Done.")
message("Annotation folder: ", out_dir)
message("Summary: ", summary_file)
