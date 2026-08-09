#!/usr/bin/env Rscript

# Step 27: Create iTOL BLAST/taxon annotations for the relaxed strict
# pct60/depth10/Nle40/nm5000 nifH tree.
#
# Tree:
#   result/trees/nifH_consensus50_strict_single_dominant_pct60_depth10_Nle40_nm5000.treefile
#
# Required BLAST best-hit table:
#   result/tables/strict_pct60_nifH_best_reference_hit_with_taxon.tsv
#
# The BLAST/taxon labels mean "closest known nifH reference hit", not confirmed
# bacterial species identity.

project_dir <- "/Users/rosa/Desktop/ALLWork/Madison/Project/ryan-nitfix/intbio-nitfix-git/symbiosis_sorted_all_sample_consensus_sequences_50percent"

tree_file <- file.path(project_dir, "result/trees/nifH_consensus50_strict_single_dominant_pct60_depth10_Nle40_nm5000.treefile")
blast_file <- file.path(project_dir, "result/tables/strict_pct60_nifH_best_reference_hit_with_taxon.tsv")
out_dir <- file.path(project_dir, "result/itol_annotations_blast_taxon_strict_pct60_Nle40_nm5000")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(paste0(label, " does not exist:\n", path), call. = FALSE)
  }
}

stop_if_missing(tree_file, "Strict pct60/Nle40 nm5000 tree file")
stop_if_missing(blast_file, "Strict pct60/Nle40 BLAST best-hit table")

read_tsv_base <- function(path) {
  read.delim(
    path,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}

read_tree_tips <- function(path) {
  if (requireNamespace("ape", quietly = TRUE)) {
    return(ape::read.tree(path)$tip.label)
  }

  txt <- paste(readLines(path, warn = FALSE), collapse = "")
  txt <- gsub(":[0-9.eE+-]+", "", txt)
  pieces <- unlist(strsplit(txt, "[(),;]+"))
  pieces <- trimws(pieces)
  pieces <- pieces[nchar(pieces) > 0]
  pieces[grepl("\\|nifH\\|consensus50\\|", pieces)]
}

tree_tips <- read_tree_tips(tree_file)
blast <- read_tsv_base(blast_file)

required_cols <- c(
  "sample_id",
  "best_reference_id",
  "best_reference_taxon",
  "percent_identity",
  "query_coverage_percent"
)
missing_cols <- setdiff(required_cols, names(blast))
if (length(missing_cols) > 0) {
  stop(paste("BLAST table is missing required columns:", paste(missing_cols, collapse = ", ")), call. = FALSE)
}

if (!"reference_coverage_percent" %in% names(blast)) {
  blast$reference_coverage_percent <- NA_real_
}

blast$percent_identity <- as.numeric(blast$percent_identity)
blast$query_coverage_percent <- as.numeric(blast$query_coverage_percent)
blast$reference_coverage_percent <- as.numeric(blast$reference_coverage_percent)

extract_genus <- function(taxon) {
  taxon <- gsub("^nifH\\|", "", taxon)
  taxon <- gsub("^Candidatus_", "Candidatus-", taxon)
  genus <- sub("_.*$", "", taxon)
  genus <- gsub("^Candidatus-", "Candidatus_", genus)
  genus[is.na(taxon) | taxon == ""] <- NA
  genus
}

blast$closest_genus <- extract_genus(blast$best_reference_taxon)
blast$blast_confidence <- ifelse(
  blast$percent_identity >= 90 & blast$query_coverage_percent >= 80,
  "high_identity90_qcov80",
  ifelse(
    blast$percent_identity >= 85 & blast$query_coverage_percent >= 70,
    "moderate_identity85_qcov70",
    "low_or_uncertain"
  )
)

ann <- merge(
  data.frame(sample_id = tree_tips, stringsAsFactors = FALSE),
  blast,
  by = "sample_id",
  all.x = TRUE,
  sort = FALSE
)

# Restore original tree order.
ann <- ann[match(tree_tips, ann$sample_id), ]

ann$in_blast_table <- ifelse(is.na(ann$best_reference_taxon), "no", "yes")
ann$tip_label_closest_taxon <- ifelse(
  is.na(ann$best_reference_taxon),
  ann$sample_id,
  paste0(
    ann$sample_id,
    " | closest: ", ann$best_reference_taxon,
    " | ", sprintf("%.1f", ann$percent_identity), "% id",
    " | ", sprintf("%.1f", ann$query_coverage_percent), "% qcov"
  )
)

make_colors <- function(values) {
  vals <- sort(unique(values[!is.na(values) & values != ""]))
  cols <- grDevices::hcl.colors(length(vals), palette = "Dark 3")
  names(cols) <- vals
  cols
}

write_colorstrip <- function(path, dataset_label, values, color_map) {
  values2 <- values
  values2[is.na(values2) | values2 == ""] <- "Unassigned"
  if (!"Unassigned" %in% names(color_map)) {
    color_map <- c(color_map, Unassigned = "#bdbdbd")
  }

  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  writeLines("DATASET_COLORSTRIP", con)
  writeLines("SEPARATOR TAB", con)
  writeLines(paste0("DATASET_LABEL\t", dataset_label), con)
  writeLines("COLOR\t#000000", con)
  writeLines(paste0("LEGEND_TITLE\t", dataset_label), con)
  writeLines(paste0("LEGEND_SHAPES\t", paste(rep("1", length(color_map)), collapse = "\t")), con)
  writeLines(paste0("LEGEND_COLORS\t", paste(unname(color_map), collapse = "\t")), con)
  writeLines(paste0("LEGEND_LABELS\t", paste(names(color_map), collapse = "\t")), con)
  writeLines("STRIP_WIDTH\t25", con)
  writeLines("MARGIN\t5", con)
  writeLines("SHOW_INTERNAL\t0", con)
  writeLines("DATA", con)

  for (i in seq_along(tree_tips)) {
    val <- values2[i]
    writeLines(paste(tree_tips[i], color_map[[val]], val, sep = "\t"), con)
  }
}

write_gradient <- function(path, dataset_label, values, color_min, color_max) {
  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  writeLines("DATASET_GRADIENT", con)
  writeLines("SEPARATOR TAB", con)
  writeLines(paste0("DATASET_LABEL\t", dataset_label), con)
  writeLines("COLOR\t#000000", con)
  writeLines(paste0("COLOR_MIN\t", color_min), con)
  writeLines(paste0("COLOR_MAX\t", color_max), con)
  writeLines("DATA", con)

  for (i in seq_along(tree_tips)) {
    if (!is.na(values[i])) {
      writeLines(paste(tree_tips[i], sprintf("%.4f", values[i]), sep = "\t"), con)
    }
  }
}

write_labels <- function(path) {
  con <- file(path, "w")
  on.exit(close(con), add = TRUE)
  writeLines("LABELS", con)
  writeLines("SEPARATOR TAB", con)
  writeLines("DATA", con)
  for (i in seq_along(tree_tips)) {
    writeLines(paste(ann$sample_id[i], ann$tip_label_closest_taxon[i], sep = "\t"), con)
  }
}

taxon_colors <- make_colors(ann$best_reference_taxon)
genus_colors <- make_colors(ann$closest_genus)
confidence_colors <- c(
  high_identity90_qcov80 = "#1a9850",
  moderate_identity85_qcov70 = "#fee08b",
  low_or_uncertain = "#d73027",
  Unassigned = "#bdbdbd"
)

write_colorstrip(
  file.path(out_dir, "strict_pct60_nm5000_blast_closest_taxon_colorstrip.txt"),
  "Closest BLAST taxon",
  ann$best_reference_taxon,
  taxon_colors
)
write_colorstrip(
  file.path(out_dir, "strict_pct60_nm5000_blast_closest_genus_colorstrip.txt"),
  "Closest BLAST genus",
  ann$closest_genus,
  genus_colors
)
write_colorstrip(
  file.path(out_dir, "strict_pct60_nm5000_blast_confidence_colorstrip.txt"),
  "BLAST assignment confidence",
  ann$blast_confidence,
  confidence_colors
)
write_gradient(
  file.path(out_dir, "strict_pct60_nm5000_blast_percent_identity_gradient.txt"),
  "BLAST percent identity",
  ann$percent_identity,
  "#fff7bc",
  "#7f0000"
)
write_gradient(
  file.path(out_dir, "strict_pct60_nm5000_blast_query_coverage_gradient.txt"),
  "BLAST query coverage",
  ann$query_coverage_percent,
  "#edf8fb",
  "#006d2c"
)
write_gradient(
  file.path(out_dir, "strict_pct60_nm5000_blast_reference_coverage_gradient.txt"),
  "BLAST reference coverage",
  ann$reference_coverage_percent,
  "#f7fcfd",
  "#084081"
)
write_labels(file.path(out_dir, "strict_pct60_nm5000_blast_taxon_tip_labels.txt"))

write.table(
  ann,
  file.path(out_dir, "strict_pct60_nm5000_blast_taxon_annotation_table.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

summary_table <- data.frame(
  metric = c(
    "tree_file",
    "blast_file",
    "tree_tips",
    "tips_with_blast_assignment",
    "tips_without_blast_assignment",
    "unique_closest_taxa",
    "unique_closest_genera",
    "mean_percent_identity",
    "min_percent_identity",
    "max_percent_identity",
    "mean_query_coverage_percent",
    "min_query_coverage_percent",
    "max_query_coverage_percent"
  ),
  value = c(
    tree_file,
    blast_file,
    length(tree_tips),
    sum(ann$in_blast_table == "yes"),
    sum(ann$in_blast_table == "no"),
    length(unique(ann$best_reference_taxon[!is.na(ann$best_reference_taxon)])),
    length(unique(ann$closest_genus[!is.na(ann$closest_genus)])),
    sprintf("%.4f", mean(ann$percent_identity, na.rm = TRUE)),
    sprintf("%.4f", min(ann$percent_identity, na.rm = TRUE)),
    sprintf("%.4f", max(ann$percent_identity, na.rm = TRUE)),
    sprintf("%.4f", mean(ann$query_coverage_percent, na.rm = TRUE)),
    sprintf("%.4f", min(ann$query_coverage_percent, na.rm = TRUE)),
    sprintf("%.4f", max(ann$query_coverage_percent, na.rm = TRUE))
  )
)

write.table(
  summary_table,
  file.path(out_dir, "strict_pct60_nm5000_blast_taxon_annotation_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

cat("Done.\n")
cat("Tree tips:", length(tree_tips), "\n")
cat("Tips with BLAST assignment:", sum(ann$in_blast_table == "yes"), "\n")
cat("Tips without BLAST assignment:", sum(ann$in_blast_table == "no"), "\n")
cat("Output folder:\n", out_dir, "\n", sep = "")
cat("\nSummary:\n")
print(summary_table, row.names = FALSE)
