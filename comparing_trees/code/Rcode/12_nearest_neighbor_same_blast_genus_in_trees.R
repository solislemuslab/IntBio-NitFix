#!/usr/bin/env Rscript

# Step 12: Sample-level nearest-neighbor analysis in functional-gene trees.
#
# Question answered:
#   For each sample/tip, is its nearest neighbor in the tree assigned to the
#   same closest BLAST genus?
#
# Why this is different from Phylocom:
#   Phylocom gives one p-value for a whole group, for example all 597
#   nifH + Mesorhizobium tips together.
#   This script gives one row per sample/tip, then summarizes the fraction of
#   samples whose nearest tree neighbor has the same closest BLAST genus.
#
# Inputs are matched to Step 10:
#   10_tip_labels_all_matched_to_blast_rows.tsv
#   10_blast_closest_genus_summary_matched_min5.tsv
#
# Outputs:
#   12_nearest_neighbor_by_sample.tsv
#   12_nearest_neighbor_summary_by_gene_genus.tsv
#   12_nearest_neighbor_summary_with_random_test.tsv
#   12_nearest_neighbor_same_genus_heatmap.png/pdf

suppressPackageStartupMessages({
  library(ape)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tibble)
  library(tidyr)
  library(ggplot2)
})

set.seed(20260830)

project_dir <- "/Users/rosa/Desktop/ALLWork/Madison/Project/ryan-nitfix/intbio-nitfix-git/comparing_trees"
tree_root <- file.path(project_dir, "result", "gene_trees_full")

in_tables <- file.path(project_dir, "result/phylocom_clustering/tables")
tip_labels_file <- file.path(in_tables, "10_tip_labels_all_matched_to_blast_rows.tsv")
matched_counts_file <- file.path(in_tables, "10_blast_closest_genus_summary_matched_min5.tsv")

out_dir <- file.path(project_dir, "result/tree_nearest_neighbor")
out_tables <- file.path(out_dir, "tables")
out_figs <- file.path(out_dir, "figures")
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figs, recursive = TRUE, showWarnings = FALSE)

genes <- c("nifH", "nifD", "nifK", "nifJ", "nodL", "nolG", "nolF", "noeA", "noeB", "nodX")
tree_sets <- c("strict", "iupac")
tree_set_labels <- c(strict = "Strict single-dominant", iupac = "Mixed-IUPAC")
tree_set_levels <- c("Strict single-dominant", "Mixed-IUPAC")

min_group_size <- as.integer(Sys.getenv("MIN_GROUP_SIZE", "5"))
randomizations <- as.integer(Sys.getenv("NEAREST_RANDOMIZATIONS", "999"))

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

same_label_for_nearest_ties <- function(labels, nearest_list) {
  map2_lgl(seq_along(labels), nearest_list, function(i, nn) {
    any(labels[nn] == labels[i], na.rm = TRUE)
  })
}

nearest_neighbor_for_tree <- function(gene, tree_set, tips_all) {
  tree_file <- get_tree_file(gene, tree_set)
  if (!file.exists(tree_file)) {
    message("Missing tree: ", tree_file)
    return(tibble())
  }

  message("Nearest-neighbor analysis: ", gene, " ", tree_set)

  tr <- read.tree(tree_file)
  tr$node.label <- NULL
  if (is.null(tr$edge.length)) tr$edge.length <- rep(1, nrow(tr$edge))
  tr$edge.length[!is.finite(tr$edge.length) | tr$edge.length <= 0] <- 1e-6

  tip_data <- tips_all %>%
    filter(.data$gene == !!gene, .data$tree_set == !!tree_set) %>%
    filter(!is.na(closest_genus_report), closest_genus_report != "") %>%
    distinct(tip_label, .keep_all = TRUE)

  keep_tips <- intersect(tr$tip.label, tip_data$tip_label)
  if (length(keep_tips) < 3) return(tibble())

  tr <- keep.tip(tr, keep_tips)
  tip_data <- tip_data %>%
    filter(tip_label %in% tr$tip.label) %>%
    mutate(tip_order = match(tip_label, tr$tip.label)) %>%
    arrange(tip_order)

  dist_mat <- cophenetic.phylo(tr)
  diag(dist_mat) <- Inf

  nearest_list <- apply(dist_mat, 1, function(x) {
    which(abs(x - min(x, na.rm = TRUE)) < .Machine$double.eps^0.5)
  })

  nearest_first <- map_int(nearest_list, 1)
  labels <- as.character(tip_data$closest_genus_report)
  same_genus <- same_label_for_nearest_ties(labels, nearest_list)

  nearest_tip_label <- rownames(dist_mat)[nearest_first]
  nearest_distance <- map2_dbl(seq_along(nearest_list), nearest_first, ~dist_mat[.x, .y])

  nearest_genus <- tip_data$closest_genus_report[match(nearest_tip_label, tip_data$tip_label)]

  tip_data %>%
    transmute(
      tree_set,
      tree_set_label = recode(tree_set, !!!tree_set_labels),
      gene,
      sample_id = sample_root_id,
      tip_label,
      closest_genus_report,
      nearest_neighbor_tip = nearest_tip_label,
      nearest_neighbor_sample_id = sample_root_id[match(nearest_tip_label, tip_label)],
      nearest_neighbor_closest_genus = nearest_genus,
      nearest_neighbor_distance = nearest_distance,
      nearest_neighbor_same_blast_genus = same_genus,
      nearest_neighbor_tie_count = lengths(nearest_list),
      sample_type,
      site,
      state,
      host_genus,
      host_tribe,
      native_status
    )
}

random_test_one_tree <- function(sample_rows) {
  if (nrow(sample_rows) < 3) return(tibble())

  labels <- as.character(sample_rows$closest_genus_report)
  groups <- sort(unique(labels))
  groups <- groups[tabulate(match(labels, groups)) >= min_group_size]
  if (length(groups) == 0) return(tibble())

  nearest_index <- match(sample_rows$nearest_neighbor_tip, sample_rows$tip_label)
  observed <- map_dfr(groups, function(g) {
    idx <- which(labels == g)
    tibble(
      closest_genus_report = g,
      n = length(idx),
      same_nearest_neighbor_n = sum(labels[nearest_index[idx]] == g, na.rm = TRUE),
      same_nearest_neighbor_percent = 100 * same_nearest_neighbor_n / n
    )
  })

  random_counts <- replicate(randomizations, {
    shuffled <- sample(labels, length(labels), replace = FALSE)
    map_int(groups, function(g) {
      idx <- which(shuffled == g)
      sum(shuffled[nearest_index[idx]] == g, na.rm = TRUE)
    })
  })

  if (is.null(dim(random_counts))) {
    random_counts <- matrix(random_counts, nrow = length(groups))
  }
  rownames(random_counts) <- groups

  observed %>%
    rowwise() %>%
    mutate(
      random_mean_same_n = mean(random_counts[closest_genus_report, ]),
      random_sd_same_n = sd(random_counts[closest_genus_report, ]),
      random_mean_same_percent = 100 * random_mean_same_n / n,
      nearest_neighbor_z = ifelse(
        random_sd_same_n > 0,
        (same_nearest_neighbor_n - random_mean_same_n) / random_sd_same_n,
        NA_real_
      ),
      p_more_same_than_random = (
        sum(random_counts[closest_genus_report, ] >= same_nearest_neighbor_n) + 1
      ) / (randomizations + 1),
      p_less_same_than_random = (
        sum(random_counts[closest_genus_report, ] <= same_nearest_neighbor_n) + 1
      ) / (randomizations + 1)
    ) %>%
    ungroup()
}

stop_if_missing(tip_labels_file, "Matched Step 10 tip-label table")
stop_if_missing(matched_counts_file, "Matched Step 10 BLAST genus count table")

tips_all <- read_tsv(tip_labels_file, show_col_types = FALSE)

sample_level <- expand_grid(gene = genes, tree_set = tree_sets) %>%
  pmap_dfr(~nearest_neighbor_for_tree(..1, ..2, tips_all))

summary_by_gene_genus <- sample_level %>%
  count(
    tree_set,
    tree_set_label,
    gene,
    closest_genus_report,
    nearest_neighbor_same_blast_genus,
    name = "sample_count"
  ) %>%
  group_by(tree_set, tree_set_label, gene, closest_genus_report) %>%
  mutate(
    n = sum(sample_count),
    same_nearest_neighbor_n = sum(sample_count[nearest_neighbor_same_blast_genus], na.rm = TRUE),
    same_nearest_neighbor_percent = 100 * same_nearest_neighbor_n / n
  ) %>%
  ungroup() %>%
  select(-nearest_neighbor_same_blast_genus, -sample_count) %>%
  distinct() %>%
  filter(n >= min_group_size) %>%
  arrange(tree_set, gene, desc(n), closest_genus_report)

random_summary <- sample_level %>%
  group_by(tree_set, tree_set_label, gene) %>%
  group_modify(~random_test_one_tree(.x)) %>%
  ungroup() %>%
  group_by(tree_set) %>%
  mutate(
    fdr_more_same_than_random = p.adjust(p_more_same_than_random, method = "BH"),
    fdr_less_same_than_random = p.adjust(p_less_same_than_random, method = "BH")
  ) %>%
  ungroup() %>%
  mutate(
    interpretation = case_when(
      fdr_more_same_than_random < 0.05 & same_nearest_neighbor_percent > random_mean_same_percent ~
        "more same-genus nearest neighbors than random",
      fdr_less_same_than_random < 0.05 & same_nearest_neighbor_percent < random_mean_same_percent ~
        "fewer same-genus nearest neighbors than random",
      TRUE ~ "similar to random / weak evidence"
    )
  )

matched_counts <- read_tsv(matched_counts_file, show_col_types = FALSE) %>%
  select(tree_set, tree_set_label, gene, closest_genus_report, blast_heatmap_n = n)

summary_with_random <- summary_by_gene_genus %>%
  left_join(
    random_summary,
    by = c("tree_set", "tree_set_label", "gene", "closest_genus_report", "n",
           "same_nearest_neighbor_n", "same_nearest_neighbor_percent")
  ) %>%
  left_join(
    matched_counts,
    by = c("tree_set", "tree_set_label", "gene", "closest_genus_report")
  ) %>%
  mutate(
    n_matches_blast_heatmap = !is.na(blast_heatmap_n) & n == blast_heatmap_n
  )

write_tsv(sample_level, file.path(out_tables, "12_nearest_neighbor_by_sample.tsv"))
write_tsv(summary_by_gene_genus, file.path(out_tables, "12_nearest_neighbor_summary_by_gene_genus.tsv"))
write_tsv(summary_with_random, file.path(out_tables, "12_nearest_neighbor_summary_with_random_test.tsv"))

check_table <- summary_with_random %>%
  count(n_matches_blast_heatmap, name = "rows")
write_tsv(check_table, file.path(out_tables, "12_nearest_neighbor_n_matches_blast_heatmap_check.tsv"))

bad_rows <- summary_with_random %>%
  filter(!n_matches_blast_heatmap)
if (nrow(bad_rows) > 0) {
  write_tsv(bad_rows, file.path(out_tables, "12_bad_rows_n_does_not_match_blast_heatmap.tsv"))
  stop(
    "Nearest-neighbor n does not match the matched BLAST heatmap n for some rows.\n",
    "Check: ", file.path(out_tables, "12_bad_rows_n_does_not_match_blast_heatmap.tsv"),
    call. = FALSE
  )
}

plot_data <- summary_with_random %>%
  mutate(
    gene = factor(gene, levels = genes),
    tree_set_label = factor(tree_set_label, levels = tree_set_levels),
    closest_genus_report = as.character(closest_genus_report),
    label_text = paste0(
      same_nearest_neighbor_n, "/", n,
      "\n", sprintf("%.0f%%", same_nearest_neighbor_percent),
      ifelse(!is.na(fdr_more_same_than_random) & fdr_more_same_than_random < 0.05, "*", "")
    )
  )

genus_order <- plot_data %>%
  group_by(closest_genus_report) %>%
  summarise(total_n = sum(n, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_n), closest_genus_report) %>%
  pull(closest_genus_report)

plot_data <- plot_data %>%
  mutate(closest_genus_report = factor(closest_genus_report, levels = rev(genus_order)))

p <- ggplot(plot_data, aes(x = gene, y = closest_genus_report, fill = same_nearest_neighbor_percent)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = label_text), size = 3.0, lineheight = 0.9) +
  facet_wrap(~ tree_set_label, ncol = 1) +
  scale_fill_gradient(low = "#f7fbff", high = "#238b45", limits = c(0, 100)) +
  labs(
    title = "Do samples have a nearest tree neighbor with the same closest BLAST genus?",
    subtitle = paste0(
      "Text = same-genus nearest neighbors / total samples and percent; * = more than random, FDR < 0.05; randomizations = ",
      randomizations
    ),
    x = "Functional gene",
    y = "Closest BLAST genus",
    fill = "% same-genus\nnearest neighbor"
  ) +
  theme_bw(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_line(color = "grey90", linewidth = 0.25),
    strip.background = element_rect(fill = "grey85", color = "grey35"),
    strip.text = element_text(face = "bold"),
    legend.position = "right"
  )

png_path <- file.path(out_figs, "12_nearest_neighbor_same_blast_genus_heatmap.png")
pdf_path <- file.path(out_figs, "12_nearest_neighbor_same_blast_genus_heatmap.pdf")

ggsave(png_path, p, width = 13, height = 10, dpi = 300)
ggsave(pdf_path, p, width = 13, height = 10)

message("Done.")
message("Sample-level table: ", file.path(out_tables, "12_nearest_neighbor_by_sample.tsv"))
message("Summary table: ", file.path(out_tables, "12_nearest_neighbor_summary_with_random_test.tsv"))
message("Check table: ", file.path(out_tables, "12_nearest_neighbor_n_matches_blast_heatmap_check.tsv"))
message("Figure PNG: ", png_path)
message("Figure PDF: ", pdf_path)
