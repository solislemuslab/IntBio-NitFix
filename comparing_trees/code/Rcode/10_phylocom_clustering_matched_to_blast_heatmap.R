#!/usr/bin/env Rscript

# Step 10: Phylocom-style clustering tests matched to the BLAST heatmap rows.
#
# Purpose
# -------
# This script tests whether tips with the same label are closer together in each
# functional-gene tree than expected by random labels.
#
# Most important design choice:
#   For closest BLAST genus, this script starts from the same deduplicated BLAST
#   assignment table used for the BLAST genus heatmap:
#     11b_blast_assignment_rows_deduplicated_strict_iupac.tsv
#
# Therefore, for a cell such as:
#   tree_set = iupac, gene = nifH, closest BLAST genus = Mesorhizobium
# the BLAST count n and the Phylocom ntaxa should be the same number.
#
# Method
# ------
# phylocomr::ph_comstruct() is used with:
#   null_model = 0
#   randomizations = 999 by default
#   abundance = FALSE
#
# Main statistic:
#   NRI = -1 * (observed MPD - mean randomized MPD) / sd randomized MPD
#
# Interpretation:
#   Positive NRI: same-label tips are closer than expected by random labels.
#   NRI near 0: same-label tips are similar to random labels.
#   Negative NRI: same-label tips are more spread out than expected by random.
#   FDR(MPD) < 0.05 and NRI > 0: significant phylogenetic clustering.
#
# Note:
#   This is a group-level test, not a per-sample test. For example, all 597
#   nifH Mesorhizobium-assigned tips are tested together as one group.
#
# References:
#   Webb et al. 2008 Phylocom. DOI: 10.1093/bioinformatics/btn358
#   phylocomr documentation: https://docs.ropensci.org/phylocomr/reference/ph_comstruct.html

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(phylocomr)
})

set.seed(20260829)

project_dir <- "/Users/rosa/Desktop/ALLWork/Madison/Project/ryan-nitfix/intbio-nitfix-git/comparing_trees"

tree_root <- file.path(project_dir, "result", "gene_trees_full")
metadata_file <- file.path(project_dir, "result", "tables", "intbio_metadata_draft4.csv")
blast_rows_file <- file.path(
  project_dir,
  "result/comparative_tree_analysis/tables/11b_blast_assignment_rows_deduplicated_strict_iupac.tsv"
)
old_blast_counts_file <- file.path(
  project_dir,
  "result/comparative_tree_analysis/tables/12_blast_closest_genus_summary_grouped_min5.tsv"
)

out_dir <- file.path(project_dir, "result", "phylocom_clustering")
out_tables <- file.path(out_dir, "tables")
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)

genes <- c("nifH", "nifD", "nifK", "nifJ", "nodL", "nolG", "nolF", "noeA", "noeB", "nodX")
tree_sets <- c("strict", "iupac")
tree_set_labels <- c(strict = "Strict single-dominant", iupac = "Mixed-IUPAC")

metadata_labels <- c(
  "sample_type",
  "site",
  "state",
  "host_genus",
  "host_tribe",
  "native_status"
)

metadata_label_names <- c(
  sample_type = "Sample type",
  site = "Site",
  state = "State",
  host_genus = "Host genus",
  host_tribe = "Host tribe",
  native_status = "Native status",
  closest_blast_genus = "Closest BLAST genus"
)

min_group_size <- as.integer(Sys.getenv("MIN_GROUP_SIZE", "5"))
randomizations <- as.integer(Sys.getenv("PHYLOCOM_RANDOMIZATIONS", "999"))

clean_sample_id <- function(x) {
  x %>%
    str_replace("\\|.*$", "") %>%
    str_replace_all("_", "-") %>%
    str_trim()
}

normalize_col <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "") %>%
    str_to_lower()
}

first_existing_col <- function(df, choices) {
  hit <- choices[choices %in% names(df)][1]
  if (is.na(hit)) {
    rep(NA_character_, nrow(df))
  } else {
    as.character(df[[hit]])
  }
}

safe_name <- function(x) {
  x %>%
    str_replace_all("[^A-Za-z0-9]+", "_") %>%
    str_replace_all("^_+|_+$", "") %>%
    str_to_lower()
}

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(paste0(label, " does not exist:\n", path), call. = FALSE)
  }
}

get_tree_file <- function(gene, tree_set) {
  if (tree_set == "strict") {
    file.path(
      tree_root, gene, "03_tree_strict_nm5000",
      paste0(gene, "_consensus50_strict_single_dominant_nm5000.treefile")
    )
  } else {
    file.path(
      tree_root, gene, "04_tree_iupac_nm5000",
      paste0(gene, "_consensus50_iupac_all_pass_nm5000.treefile")
    )
  }
}

read_metadata <- function() {
  stop_if_missing(metadata_file, "Metadata file")

  metadata <- read_csv(metadata_file, show_col_types = FALSE)
  names(metadata) <- normalize_col(names(metadata))

  sample_col <- names(metadata)[names(metadata) %in% c("sample_id", "sample", "sampleid")][1]
  if (is.na(sample_col)) stop("Could not find sample ID column in metadata.", call. = FALSE)

  metadata %>%
    mutate(
      sample_type = first_existing_col(., c("sample_type", "type")),
      host_genus = first_existing_col(., c("host_genus", "genus", "plant_genus")),
      host_tribe = first_existing_col(., c("host_tribe", "tribe", "plant_tribe")),
      native_status = first_existing_col(., c("native_status", "native", "native_non_native")),
      site = first_existing_col(., c("site")),
      state = first_existing_col(., c("state"))
    ) %>%
    mutate(sample_root_id = clean_sample_id(.data[[sample_col]])) %>%
    distinct(sample_root_id, .keep_all = TRUE) %>%
    mutate(across(any_of(metadata_labels), ~replace_na(as.character(.x), "Unknown")))
}

read_blast_rows <- function() {
  stop_if_missing(blast_rows_file, "Deduplicated BLAST assignment file")

  blast_rows <- read_tsv(blast_rows_file, show_col_types = FALSE)
  required <- c("tree_set", "gene", "sample_id")
  missing_required <- setdiff(required, names(blast_rows))
  if (length(missing_required) > 0) {
    stop("BLAST assignment file is missing columns: ", paste(missing_required, collapse = ", "), call. = FALSE)
  }

  if (!"closest_genus_report" %in% names(blast_rows)) {
    if ("closest_genus" %in% names(blast_rows)) {
      blast_rows <- blast_rows %>% mutate(closest_genus_report = closest_genus)
    } else if ("best_reference_genus" %in% names(blast_rows)) {
      blast_rows <- blast_rows %>% mutate(closest_genus_report = best_reference_genus)
    } else {
      stop("Could not find closest genus column in BLAST assignment file.", call. = FALSE)
    }
  }

  blast_rows %>%
    filter(gene %in% genes, tree_set %in% tree_sets) %>%
    mutate(
      sample_id = as.character(sample_id),
      sample_root_id = clean_sample_id(sample_id),
      tree_set_label = recode(tree_set, !!!tree_set_labels),
      closest_genus_report = as.character(closest_genus_report),
      closest_blast_genus = closest_genus_report
    ) %>%
    filter(!is.na(closest_blast_genus), closest_blast_genus != "") %>%
    distinct(tree_set, gene, sample_id, .keep_all = TRUE)
}

match_blast_rows_to_tree <- function(tree, blast_rows, gene, tree_set) {
  gene_arg <- gene
  tree_set_arg <- tree_set
  tree_set_label_arg <- unname(tree_set_labels[[tree_set_arg]])

  tip_table <- tibble(
    tip_label = tree$tip.label,
    sample_root_id = clean_sample_id(tip_label)
  )

  blast_sub <- blast_rows %>%
    filter(.data$gene == !!gene_arg, .data$tree_set == !!tree_set_arg)

  blast_exact <- blast_sub %>%
    select(-sample_root_id) %>%
    rename(tip_label = sample_id)

  exact <- tip_table %>%
    left_join(blast_exact, by = "tip_label")

  needs_clean_match <- is.na(exact$closest_blast_genus)
  if (any(needs_clean_match)) {
    clean_blast <- blast_sub %>%
      group_by(sample_root_id) %>%
      slice(1) %>%
      ungroup() %>%
      select(
        sample_root_id,
        clean_sample_id_from_blast = sample_id,
        closest_blast_genus,
        closest_genus_report,
        any_of(c(
          "closest_taxon",
          "closest_genus",
          "best_reference_id",
          "best_reference_taxon",
          "best_reference_genus",
          "percent_identity",
          "alignment_length",
          "query_coverage_percent",
          "evalue",
          "bitscore"
        ))
      )

    clean_matched <- tip_table %>%
      filter(tip_label %in% exact$tip_label[needs_clean_match]) %>%
      left_join(clean_blast, by = "sample_root_id") %>%
      mutate(sample_id = tip_label)

    exact <- exact %>%
      filter(!needs_clean_match) %>%
      bind_rows(clean_matched)
  }

  exact %>%
    mutate(
      gene = gene_arg,
      tree_set = tree_set_arg,
      tree_set_label = tree_set_label_arg,
      sample_id = tip_label
    ) %>%
    filter(!is.na(closest_blast_genus), closest_blast_genus != "") %>%
    distinct(tree_set, gene, sample_id, .keep_all = TRUE)
}

standardize_phylocom_output <- function(result) {
  result <- as_tibble(result)
  names(result) <- normalize_col(names(result))

  if (!"plot" %in% names(result) && "sample" %in% names(result)) {
    result <- result %>% rename(plot = sample)
  }

  for (col in c("ntaxa", "mpd", "mpd_rnd", "mpd_sd", "nri", "mpd_rankhi",
                "mntd", "mntd_rnd", "mntd_sd", "nti", "mntd_rankhi", "runs")) {
    if (!col %in% names(result)) result[[col]] <- NA_real_
  }

  result
}

run_phylocom_one_label <- function(tree, tip_table, label_col, gene, tree_set) {
  gene_arg <- gene
  tree_set_arg <- tree_set
  tree_set_label_arg <- unname(tree_set_labels[[tree_set_arg]])
  label_type_report_arg <- unname(metadata_label_names[[label_col]] %||% label_col)

  dat <- tip_table %>%
    mutate(
      label_value = as.character(.data[[label_col]]),
      label_value = ifelse(is.na(label_value) | label_value == "", "Unknown", label_value)
    ) %>%
    filter(label_value != "Unknown") %>%
    add_count(label_value, name = "expected_ntaxa") %>%
    filter(expected_ntaxa >= min_group_size)

  if (nrow(dat) < 4 || n_distinct(dat$label_value) < 1) {
    return(tibble())
  }

  # IQ-TREE stores support values as internal node labels. The external
  # Phylocom binary can fail on those labels, and they are not used for MPD.
  tree$node.label <- NULL
  if (is.null(tree$edge.length)) tree$edge.length <- rep(1, nrow(tree$edge))
  tree$edge.length[!is.finite(tree$edge.length) | tree$edge.length <= 0] <- 1e-6

  keep_tips <- intersect(tree$tip.label, dat$tip_label)
  tr <- keep.tip(tree, keep_tips)
  tr$node.label <- NULL
  if (is.null(tr$edge.length)) tr$edge.length <- rep(1, nrow(tr$edge))
  tr$edge.length[!is.finite(tr$edge.length) | tr$edge.length <= 0] <- 1e-6

  tip_map <- tibble(
    tip_label = tr$tip.label,
    safe_tip = paste0("t", seq_along(tr$tip.label))
  )
  tr$tip.label <- tip_map$safe_tip

  sample_df <- dat %>%
    inner_join(tip_map, by = "tip_label") %>%
    transmute(
      group = as.character(safe_name(paste0(label_col, "__", label_value))),
      abundance = as.integer(1),
      taxon = as.character(safe_tip)
    ) %>%
    arrange(group, taxon) %>%
    as.data.frame(stringsAsFactors = FALSE)

  sample_df[[1]] <- as.character(sample_df[[1]])
  sample_df[[2]] <- as.integer(sample_df[[2]])
  sample_df[[3]] <- as.character(sample_df[[3]])

  group_map <- dat %>%
    distinct(label_value, expected_ntaxa) %>%
    mutate(plot = safe_name(paste0(label_col, "__", label_value))) %>%
    group_by(plot) %>%
    summarise(
      label_value = paste(sort(unique(label_value)), collapse = " / "),
      expected_ntaxa = sum(expected_ntaxa[!duplicated(label_value)]),
      .groups = "drop"
    )

  result <- tryCatch(
    ph_comstruct(
      sample = sample_df,
      phylo = tr,
      null_model = 0,
      randomizations = randomizations,
      abundance = FALSE
    ),
    error = function(e) {
      message("Phylocom failed for ", gene, " ", tree_set, " ", label_col, ": ", conditionMessage(e))
      return(NULL)
    }
  )

  if (is.null(result)) {
    return(tibble(
      gene = gene_arg,
      tree_set = tree_set_arg,
      tree_set_label = tree_set_label_arg,
      label_type = label_col,
      label_type_report = label_type_report_arg,
      label_value = NA_character_,
      expected_ntaxa = NA_integer_,
      phylocom_status = "failed_phylocom_binary_crash_or_error",
      ntaxa = NA_integer_,
      mpd = NA_real_,
      mpd_random = NA_real_,
      mpd_sd = NA_real_,
      nri = NA_real_,
      mpd_rankhi = NA_real_,
      mntd = NA_real_,
      mntd_random = NA_real_,
      mntd_sd = NA_real_,
      nti = NA_real_,
      mntd_rankhi = NA_real_,
      runs = NA_real_,
      p_mpd_cluster = NA_real_,
      p_mntd_cluster = NA_real_
    ))
  }

  standardize_phylocom_output(result) %>%
    left_join(group_map, by = "plot") %>%
    mutate(
      gene = gene_arg,
      tree_set = tree_set_arg,
      tree_set_label = tree_set_label_arg,
      label_type = label_col,
      label_type_report = label_type_report_arg,
      phylocom_status = "ok",
      p_mpd_cluster = (mpd_rankhi + 1) / (runs + 1),
      p_mntd_cluster = (mntd_rankhi + 1) / (runs + 1),
      ntaxa_matches_expected = !is.na(ntaxa) & ntaxa == expected_ntaxa,
      .before = 1
    )
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x)) y else x
}

metadata <- read_metadata()
blast_rows <- read_blast_rows()

run_one_tree <- function(gene, tree_set) {
  tree_file <- get_tree_file(gene, tree_set)
  if (!file.exists(tree_file)) {
    message("Missing tree: ", tree_file)
    return(list(results = tibble(), tip_table = tibble()))
  }

  message("Running Phylocom-style test: ", gene, " ", tree_set)
  tree <- read.tree(tree_file)

  tip_table <- match_blast_rows_to_tree(tree, blast_rows, gene, tree_set) %>%
    left_join(metadata, by = "sample_root_id")

  missing_tip_n <- length(setdiff(tree$tip.label, tip_table$tip_label))
  if (missing_tip_n > 0) {
    message("  Note: ", missing_tip_n, " tree tips did not have a BLAST genus row in 11b and were not tested for BLAST genus.")
  }

  write_tsv(
    tip_table,
    file.path(out_tables, paste0("10_tip_labels_matched_", gene, "_", tree_set, ".tsv"))
  )

  labels_to_test <- c(metadata_labels, "closest_blast_genus")
  labels_to_test <- labels_to_test[labels_to_test %in% names(tip_table)]

  results <- map_dfr(
    labels_to_test,
    ~run_phylocom_one_label(tree, tip_table, .x, gene, tree_set)
  )

  list(results = results, tip_table = tip_table)
}

all_runs <- expand_grid(gene = genes, tree_set = tree_sets) %>%
  pmap(~run_one_tree(..1, ..2))

all_results <- map_dfr(all_runs, "results") %>%
  group_by(tree_set, label_type) %>%
  mutate(
    fdr_mpd = p.adjust(p_mpd_cluster, method = "BH"),
    fdr_mntd = p.adjust(p_mntd_cluster, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    interpretation = case_when(
      fdr_mpd < 0.05 & nri > 0 ~ "significant clustering by MPD",
      fdr_mntd < 0.05 & nti > 0 ~ "significant clustering by MNTD",
      TRUE ~ "not significant / weak evidence"
    )
  )

all_tip_labels <- map_dfr(all_runs, "tip_table")

blast_counts_matched <- all_tip_labels %>%
  filter(!is.na(closest_genus_report), closest_genus_report != "") %>%
  count(tree_set, tree_set_label, gene, closest_genus_report, name = "n") %>%
  group_by(tree_set, gene) %>%
  mutate(
    gene_total = sum(n),
    percent = 100 * n / gene_total
  ) %>%
  ungroup() %>%
  filter(n >= min_group_size) %>%
  arrange(tree_set, gene, desc(n), closest_genus_report)

phylocom_genus_check <- all_results %>%
  filter(label_type == "closest_blast_genus") %>%
  transmute(
    tree_set,
    tree_set_label,
    gene,
    closest_genus_report = label_value,
    phylocom_ntaxa = as.integer(ntaxa),
    expected_ntaxa = as.integer(expected_ntaxa),
    ntaxa_matches_expected
  )

blast_phylocom_check <- blast_counts_matched %>%
  left_join(
    phylocom_genus_check,
    by = c("tree_set", "tree_set_label", "gene", "closest_genus_report")
  ) %>%
  mutate(
    n_matches_phylocom_ntaxa = !is.na(phylocom_ntaxa) & n == phylocom_ntaxa,
    n_matches_expected_ntaxa = !is.na(expected_ntaxa) & n == expected_ntaxa
  )

old_count_check <- tibble()
if (file.exists(old_blast_counts_file)) {
  old_count_check <- read_tsv(old_blast_counts_file, show_col_types = FALSE) %>%
    filter(gene %in% genes, tree_set %in% tree_sets) %>%
    select(tree_set, tree_set_label, gene, closest_genus_report, old_n = n, old_gene_total = gene_total, old_percent = percent) %>%
    right_join(
      blast_counts_matched,
      by = c("tree_set", "tree_set_label", "gene", "closest_genus_report")
    ) %>%
    mutate(
      old_n_matches_new_n = !is.na(old_n) & old_n == n,
      old_gene_total_matches_new_gene_total = !is.na(old_gene_total) & old_gene_total == gene_total
    )
}

if (!"mpd_random" %in% names(all_results)) all_results$mpd_random <- NA_real_
if (!"mntd_random" %in% names(all_results)) all_results$mntd_random <- NA_real_
if (!"mpd_rnd" %in% names(all_results)) all_results$mpd_rnd <- NA_real_
if (!"mntd_rnd" %in% names(all_results)) all_results$mntd_rnd <- NA_real_

report_table <- all_results %>%
  mutate(
    mpd_random_report = coalesce(as.numeric(mpd_random), as.numeric(mpd_rnd)),
    mntd_random_report = coalesce(as.numeric(mntd_random), as.numeric(mntd_rnd))
  ) %>%
  transmute(
    gene,
    tree_set,
    tree_set_label,
    label_type,
    label_type_report,
    label_value,
    expected_ntaxa,
    ntaxa = as.integer(ntaxa),
    ntaxa_matches_expected,
    mpd,
    mpd_random = mpd_random_report,
    mpd_sd,
    nri,
    p_mpd_cluster,
    fdr_mpd,
    mntd,
    mntd_random = mntd_random_report,
    mntd_sd,
    nti,
    p_mntd_cluster,
    fdr_mntd,
    randomizations,
    min_group_size,
    phylocom_status,
    interpretation
  ) %>%
  arrange(tree_set, gene, label_type, fdr_mpd)

write_tsv(all_results, file.path(out_tables, "10_phylocom_comstruct_all_results_matched.tsv"))
write_tsv(report_table, file.path(out_tables, "10_phylocom_report_table_matched.tsv"))
write_tsv(all_tip_labels, file.path(out_tables, "10_tip_labels_all_matched_to_blast_rows.tsv"))
write_tsv(blast_counts_matched, file.path(out_tables, "10_blast_closest_genus_summary_matched_min5.tsv"))
write_tsv(blast_phylocom_check, file.path(out_tables, "10_blast_count_vs_phylocom_ntaxa_check.tsv"))
if (nrow(old_count_check) > 0) {
  write_tsv(old_count_check, file.path(out_tables, "10_old_vs_new_blast_count_check.tsv"))
}

problem_rows <- blast_phylocom_check %>%
  filter(!n_matches_phylocom_ntaxa | !n_matches_expected_ntaxa)

if (nrow(problem_rows) > 0) {
  warning(
    "Some BLAST n values do not match Phylocom ntaxa. See: ",
    file.path(out_tables, "10_blast_count_vs_phylocom_ntaxa_check.tsv"),
    call. = FALSE
  )
} else {
  message("All matched BLAST genus counts agree with Phylocom ntaxa.")
}

message("Done.")
message("Phylocom report table: ", file.path(out_tables, "10_phylocom_report_table_matched.tsv"))
message("Matched BLAST count table: ", file.path(out_tables, "10_blast_closest_genus_summary_matched_min5.tsv"))
message("Count check table: ", file.path(out_tables, "10_blast_count_vs_phylocom_ntaxa_check.tsv"))
