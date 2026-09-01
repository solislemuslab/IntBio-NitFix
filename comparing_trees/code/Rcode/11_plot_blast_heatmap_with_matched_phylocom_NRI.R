#!/usr/bin/env Rscript

# Step 11: Plot the BLAST genus heatmap with matched Phylocom NRI values.
#
# This script does not run Phylocom. It only plots Step 10 results.
#
# It checks that:
#   n in the BLAST heatmap cell == ntaxa tested by phylocomr::ph_comstruct()
#
# If those values do not match, the script stops before making the figure.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(stringr)
})

project_dir <- "/Users/rosa/Desktop/ALLWork/Madison/Project/ryan-nitfix/intbio-nitfix-git/comparing_trees"

out_tables <- file.path(project_dir, "result/phylocom_clustering/tables")
out_figs <- file.path(project_dir, "result/phylocom_clustering/figures")
dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figs, recursive = TRUE, showWarnings = FALSE)

blast_counts_file <- file.path(out_tables, "10_blast_closest_genus_summary_matched_min5.tsv")
phylocom_file <- file.path(out_tables, "10_phylocom_report_table_matched.tsv")

stop_if_missing <- function(path, label) {
  if (!file.exists(path)) {
    stop(paste0(label, " does not exist:\n", path), call. = FALSE)
  }
}

stop_if_missing(blast_counts_file, "Matched BLAST genus count table from Step 10")
stop_if_missing(phylocom_file, "Matched Phylocom report table from Step 10")

genes <- c("nifH", "nifD", "nifK", "nifJ", "nodL", "nolG", "nolF", "noeA", "noeB", "nodX")
tree_set_levels <- c("Strict single-dominant", "Mixed-IUPAC")

blast_counts <- read_tsv(blast_counts_file, show_col_types = FALSE) %>%
  filter(gene %in% genes, tree_set %in% c("strict", "iupac")) %>%
  mutate(
    gene = factor(gene, levels = genes),
    tree_set_label = factor(tree_set_label, levels = tree_set_levels),
    closest_genus_report = as.character(closest_genus_report)
  )

phylocom_genus <- read_tsv(phylocom_file, show_col_types = FALSE) %>%
  filter(
    gene %in% genes,
    tree_set %in% c("strict", "iupac"),
    label_type == "closest_blast_genus",
    !is.na(label_value)
  ) %>%
  transmute(
    tree_set,
    gene,
    closest_genus_report = as.character(label_value),
    phylocom_ntaxa = as.integer(ntaxa),
    expected_ntaxa = as.integer(expected_ntaxa),
    nri = as.numeric(nri),
    mpd = as.numeric(mpd),
    mpd_random = as.numeric(mpd_random),
    p_mpd_cluster = as.numeric(p_mpd_cluster),
    fdr_mpd = as.numeric(fdr_mpd),
    nti = as.numeric(nti),
    p_mntd_cluster = as.numeric(p_mntd_cluster),
    fdr_mntd = as.numeric(fdr_mntd),
    randomizations = as.integer(randomizations),
    min_group_size = as.integer(min_group_size),
    interpretation
  )

plot_data <- blast_counts %>%
  left_join(phylocom_genus, by = c("tree_set", "gene", "closest_genus_report")) %>%
  mutate(
    n_matches_phylocom_ntaxa = !is.na(phylocom_ntaxa) & n == phylocom_ntaxa,
    n_matches_expected_ntaxa = !is.na(expected_ntaxa) & n == expected_ntaxa
  )

join_check <- plot_data %>%
  count(n_matches_phylocom_ntaxa, n_matches_expected_ntaxa, name = "rows")

write_tsv(join_check, file.path(out_tables, "11_matched_heatmap_join_check.tsv"))
write_tsv(plot_data, file.path(out_tables, "11_blast_genus_counts_with_matched_phylocom_NRI.tsv"))

bad_rows <- plot_data %>%
  filter(!n_matches_phylocom_ntaxa | !n_matches_expected_ntaxa)

if (nrow(bad_rows) > 0) {
  write_tsv(bad_rows, file.path(out_tables, "11_bad_rows_n_does_not_match_phylocom_ntaxa.tsv"))
  stop(
    "The figure was not made because some BLAST n values do not match Phylocom ntaxa.\n",
    "Check: ", file.path(out_tables, "11_bad_rows_n_does_not_match_phylocom_ntaxa.tsv"),
    call. = FALSE
  )
}

# Keep the same genus order as the original BLAST heatmap: most common genera
# across all genes/tree sets are easiest to compare visually.
genus_order <- plot_data %>%
  group_by(closest_genus_report) %>%
  summarise(total_n = sum(n, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(total_n), closest_genus_report) %>%
  pull(closest_genus_report)

plot_data <- plot_data %>%
  mutate(
    gene = factor(as.character(gene), levels = genes),
    tree_set_label = factor(as.character(tree_set_label), levels = tree_set_levels),
    closest_genus_report = factor(closest_genus_report, levels = rev(genus_order)),
    nri_label = paste0(
      n,
      "\nNRI=", ifelse(is.na(nri), "NA", sprintf("%.1f", nri)),
      ifelse(!is.na(fdr_mpd) & fdr_mpd < 0.05 & !is.na(nri) & nri > 0, "*", "")
    )
  )

p <- ggplot(plot_data, aes(x = gene, y = closest_genus_report, fill = n)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_text(aes(label = nri_label), size = 3.0, lineheight = 0.9) +
  facet_wrap(~ tree_set_label, ncol = 1) +
  scale_fill_gradient(low = "#f7fbff", high = "#08519c") +
  labs(
    title = "Gene-by-genus BLAST assignment summary with matched Phylocom NRI",
    subtitle = "Number = sample consensus sequences; NRI from phylocomr::ph_comstruct; * = FDR(MPD) < 0.05 and NRI > 0",
    x = "Functional gene",
    y = "Closest BLAST genus",
    fill = "Samples"
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

png_path <- file.path(out_figs, "11_blast_genus_by_gene_heatmap_min5_with_matched_phylocom_NRI.png")
pdf_path <- file.path(out_figs, "11_blast_genus_by_gene_heatmap_min5_with_matched_phylocom_NRI.pdf")

ggsave(png_path, p, width = 13, height = 10, dpi = 300)
ggsave(pdf_path, p, width = 13, height = 10)

message("Done.")
message("Figure PNG: ", png_path)
message("Figure PDF: ", pdf_path)
message("Joined table: ", file.path(out_tables, "11_blast_genus_counts_with_matched_phylocom_NRI.tsv"))
message("Join check: ", file.path(out_tables, "11_matched_heatmap_join_check.tsv"))
