# ============================================================
# 02_metadata_and_report_comparisons_consensus50.R
#
# Purpose:
#   Build report-ready comparisons for the consensus50 functional-gene trees.
#   This script turns tree/coverage/BLAST outputs into interpretable tables
#   and figures for Ryan/Ahmed's questions:
#     1. Which genes are recovered across sample types and host metadata?
#     2. Do nif, canonical nod, and accessory nodulation/symbiosis genes differ?
#     3. Which closest BLAST genera dominate each gene and sample type?
#     4. Which metadata categories are worth testing later with tree-based
#        phylogenetic clustering?
#
# Important:
#   BLAST labels are closest known reference hits, not confirmed species IDs.
# ============================================================

project_dir <- "/Users/rosa/Desktop/ALLWork/Madison/Project/ryan-nitfix/intbio-nitfix-git/symbiosis_sorted_all_sample_consensus_sequences_50percent"

result_dir <- file.path(project_dir, "result")
tables_dir <- file.path(result_dir, "tables")
comparative_dir <- file.path(result_dir, "comparative_tree_analysis")
comparative_tables <- file.path(comparative_dir, "tables")

out_dir <- file.path(result_dir, "metadata_tree_comparison")
out_tables <- file.path(out_dir, "tables")
out_figs <- file.path(out_dir, "figures")

dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figs, recursive = TRUE, showWarnings = FALSE)

genes <- c("nifH", "nifD", "nifK", "nifJ", "nodL", "nolG", "nolF", "noeA", "noeB", "nodX")
nif_genes <- c("nifH", "nifD", "nifK", "nifJ")
canonical_nod_genes <- c("nodL", "nodX")
accessory_nodulation_genes <- c("nolG", "nolF", "noeA", "noeB")

gene_group_levels <- c("nif", "canonical nod", "accessory nodulation/symbiosis-related")

gene_group_of <- function(x) {
  dplyr::case_when(
    x %in% nif_genes ~ "nif",
    x %in% canonical_nod_genes ~ "canonical nod",
    x %in% accessory_nodulation_genes ~ "accessory nodulation/symbiosis-related",
    TRUE ~ "other"
  )
}

# ----------------------------
# 1. Packages
# ----------------------------

need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Please install this R package first: ", pkg,
                "\nExample: install.packages('", pkg, "')"), call. = FALSE)
  }
}

for (pkg in c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "purrr", "forcats", "scales")) {
  need_pkg(pkg)
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)
library(forcats)
library(scales)

theme_set(theme_bw(base_size = 13))

# ----------------------------
# 2. Helpers
# ----------------------------

clean_sample_id <- function(x) {
  x %>%
    str_replace_all("_", "-") %>%
    str_replace("-R$", "-Rh") %>%
    str_replace("-N$", "-No")
}

sample_from_tip <- function(x) {
  str_replace(x, "\\|.*$", "")
}

sample_type_from_id <- function(x) {
  case_when(
    x %in% c("MC-1", "MC-2") ~ "MC",
    str_detect(x, "-No$") ~ "No",
    str_detect(x, "-Rh$") ~ "Rh",
    str_detect(x, "-Ro$") ~ "Ro",
    TRUE ~ "Unknown"
  )
}

safe_read_tsv <- function(path) {
  if (!file.exists(path)) stop(paste0("Missing required file:\n", path), call. = FALSE)
  read_tsv(path, show_col_types = FALSE)
}

metadata_file <- file.path(tables_dir, "intbio_metadata_draft4.csv")
if (!file.exists(metadata_file)) {
  metadata_file <- file.path(tables_dir, "intbio_metadata_draft3.csv")
}
if (!file.exists(metadata_file)) {
  stop("Could not find intbio_metadata_draft4.csv or intbio_metadata_draft3.csv in result/tables.", call. = FALSE)
}

coverage_file <- file.path(tables_dir, "consensus50_gene_coverage_all_samples.tsv")
tree_summary_file <- file.path(comparative_tables, "05_tree_summary_strict_and_mixed_iupac.tsv")
blast_file <- file.path(comparative_tables, "11b_blast_assignment_rows_deduplicated_strict_iupac.tsv")

coverage <- safe_read_tsv(coverage_file)
tree_summary <- safe_read_tsv(tree_summary_file)
blast <- safe_read_tsv(blast_file)
metadata <- read_csv(metadata_file, show_col_types = FALSE)

metadata_norm <- metadata %>%
  mutate(
    sample = clean_sample_id(`Sample ID`),
    metadata_type = Type,
    host_family = Family,
    host_genus = Genus,
    host_species = Species,
    host_tribe = Tribe,
    site = Site,
    county = County,
    state = State,
    native_status = native,
    latitude = Lat,
    longitude = Long,
    field_collection_date = `Date (field collection)`
  ) %>%
  select(sample, metadata_type, host_family, host_genus, host_species, host_tribe,
         site, county, state, native_status, latitude, longitude, field_collection_date) %>%
  distinct(sample, .keep_all = TRUE)

coverage_good <- coverage %>%
  filter(gene %in% genes) %>%
  mutate(
    good_coverage = percent_covered >= 80 & mean_depth >= 10,
    sample_type = sample_type_from_id(sample),
    gene_group = factor(gene_group_of(gene), levels = gene_group_levels)
  ) %>%
  left_join(metadata_norm, by = "sample")

tree_summary2 <- tree_summary %>%
  filter(gene %in% genes) %>%
  mutate(
    gene_group = factor(gene_group_of(gene), levels = gene_group_levels),
    tree_set_label = recode(tree_set,
                            strict = "Strict single-dominant",
                            iupac = "Mixed-IUPAC",
                            .default = tree_set)
  )

blast2 <- blast %>%
  filter(gene %in% genes, tree_set %in% c("strict", "iupac")) %>%
  mutate(
    sample = sample_from_tip(sample_id),
    sample_type = sample_type_from_id(sample),
    gene_group = factor(gene_group_of(gene), levels = gene_group_levels),
    tree_set_label = recode(tree_set,
                            strict = "Strict single-dominant",
                            iupac = "Mixed-IUPAC",
                            .default = tree_set),
    closest_genus = ifelse(is.na(closest_genus) | closest_genus == "", "Unassigned", closest_genus)
  ) %>%
  left_join(metadata_norm, by = "sample")

# ----------------------------
# 3. Report table: gene categories
# ----------------------------

gene_category_table <- tibble(
  gene = genes,
  gene_group = gene_group_of(genes),
  report_wording = case_when(
    gene %in% nif_genes ~ "nif nitrogen-fixation gene",
    gene %in% canonical_nod_genes ~ "canonical nod gene",
    gene %in% accessory_nodulation_genes ~ "accessory nodulation/symbiosis-related gene",
    TRUE ~ "other"
  )
)

write_tsv(gene_category_table, file.path(out_tables, "01_gene_categories_for_report.tsv"))

# ----------------------------
# 4. Recovery by metadata category
# ----------------------------

summarise_recovery <- function(data, metadata_col, min_group_n = 5) {
  metadata_sym <- rlang::sym(metadata_col)

  data %>%
    mutate(metadata_value = as.character(!!metadata_sym)) %>%
    mutate(metadata_value = ifelse(is.na(metadata_value) | metadata_value == "", "Unknown", metadata_value)) %>%
    group_by(gene, gene_group, metadata_category = metadata_col, metadata_value) %>%
    summarise(
      total_samples = n_distinct(sample),
      good_samples = n_distinct(sample[good_coverage]),
      percent_good = 100 * good_samples / total_samples,
      .groups = "drop"
    ) %>%
    filter(total_samples >= min_group_n)
}

metadata_categories <- c(
  "sample_type",
  "metadata_type",
  "host_family",
  "host_genus",
  "host_tribe",
  "site",
  "state",
  "native_status"
)

recovery_by_metadata <- map_dfr(metadata_categories, ~ summarise_recovery(coverage_good, .x, min_group_n = 5))

write_tsv(recovery_by_metadata, file.path(out_tables, "02_gene_recovery_by_metadata_category_pct80_depth10.tsv"))

sample_type_recovery <- recovery_by_metadata %>%
  filter(metadata_category == "sample_type") %>%
  mutate(gene = factor(gene, levels = genes))

p_sample_type <- sample_type_recovery %>%
  ggplot(aes(x = gene, y = percent_good, fill = metadata_value)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = paste0(round(percent_good, 1), "%")),
            position = position_dodge(width = 0.8), vjust = -0.25, size = 2.8) +
  scale_fill_manual(values = c(MC = "#E41A1C", No = "#377EB8", Rh = "#FF9900", Ro = "#7B52AB", Unknown = "gray70")) +
  labs(
    title = "Functional-gene recovery by sample type",
    subtitle = "Good coverage: percent covered >= 80% and mean depth >= 10",
    x = "Gene",
    y = "Samples with good coverage (%)",
    fill = "Sample type"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_figs, "01_gene_recovery_percent_by_sample_type.png"), p_sample_type, width = 12, height = 6, dpi = 300)

selected_metadata_for_heatmap <- c("host_tribe", "site", "native_status")

top_recovery_groups <- recovery_by_metadata %>%
  filter(metadata_category %in% selected_metadata_for_heatmap) %>%
  group_by(metadata_category, metadata_value) %>%
  summarise(max_total_samples = max(total_samples), .groups = "drop") %>%
  group_by(metadata_category) %>%
  slice_max(max_total_samples, n = 12, with_ties = FALSE) %>%
  ungroup()

recovery_heatmap_data <- recovery_by_metadata %>%
  semi_join(top_recovery_groups, by = c("metadata_category", "metadata_value")) %>%
  mutate(
    gene = factor(gene, levels = genes),
    metadata_label = paste0(metadata_category, ": ", metadata_value)
  )

p_recovery_meta <- recovery_heatmap_data %>%
  ggplot(aes(x = gene, y = fct_reorder(metadata_label, percent_good, .fun = max), fill = percent_good)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = paste0(round(percent_good), "%")), size = 2.7) +
  facet_wrap(~ metadata_category, scales = "free_y", ncol = 1) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", limits = c(0, 100)) +
  labs(
    title = "Gene recovery across host and geography metadata",
    subtitle = "Top metadata groups by sample count; groups with at least 5 samples are included",
    x = "Gene",
    y = "Metadata group",
    fill = "% good"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_figs, "02_gene_recovery_heatmap_host_geography.png"), p_recovery_meta, width = 12, height = 11, dpi = 300)

# ----------------------------
# 5. Strict vs mixed and tree reliability for report
# ----------------------------

tree_reliability <- tree_summary2 %>%
  select(gene, gene_group, tree_set, tree_set_label, input_alignment_sequences,
         alignment_length, iqtree_taxa_or_tips, best_model, tree_length,
         final_bootstrap_correlation, bootstrap_status) %>%
  arrange(gene_group, gene, tree_set)

write_tsv(tree_reliability, file.path(out_tables, "03_tree_reliability_summary_for_report.tsv"))

p_tree_reliability <- tree_reliability %>%
  mutate(
    gene = factor(gene, levels = genes),
    bootstrap_status = ifelse(is.na(bootstrap_status), "missing", bootstrap_status)
  ) %>%
  ggplot(aes(x = gene, y = input_alignment_sequences, fill = bootstrap_status)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = input_alignment_sequences), vjust = -0.25, size = 3) +
  facet_wrap(~ tree_set_label, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c(converged = "#1B9E77", not_converged = "#D95F02", missing = "gray70")) +
  labs(
    title = "Tree input size and bootstrap convergence",
    subtitle = "Strict trees are primary; mixed-IUPAC trees are sensitivity/exploratory when not converged",
    x = "Gene",
    y = "Input consensus sequences",
    fill = "Bootstrap status"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_figs, "03_tree_input_size_and_convergence.png"), p_tree_reliability, width = 11, height = 8, dpi = 300)

# ----------------------------
# 6. BLAST genera by sample type and gene
# ----------------------------

blast_genus_by_type <- blast2 %>%
  count(tree_set, tree_set_label, gene, gene_group, sample_type, closest_genus, name = "n") %>%
  group_by(tree_set, tree_set_label, gene, sample_type) %>%
  mutate(
    gene_sample_type_total = sum(n),
    percent = 100 * n / gene_sample_type_total
  ) %>%
  ungroup() %>%
  arrange(tree_set, gene, sample_type, desc(n))

write_tsv(blast_genus_by_type, file.path(out_tables, "04_blast_closest_genus_by_gene_sample_type.tsv"))

top_blast_genera <- blast2 %>%
  count(closest_genus, name = "total_n") %>%
  slice_max(total_n, n = 12, with_ties = FALSE) %>%
  pull(closest_genus)

blast_genus_plot_data <- blast_genus_by_type %>%
  mutate(
    closest_genus_report = ifelse(closest_genus %in% top_blast_genera, closest_genus, "Other"),
    gene = factor(gene, levels = genes)
  ) %>%
  group_by(tree_set, tree_set_label, gene, sample_type, closest_genus_report) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  group_by(tree_set, tree_set_label, gene, sample_type) %>%
  mutate(percent = 100 * n / sum(n)) %>%
  ungroup()

p_blast_type <- blast_genus_plot_data %>%
  filter(tree_set == "strict") %>%
  ggplot(aes(x = gene, y = percent, fill = closest_genus_report)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.15) +
  facet_wrap(~ sample_type, ncol = 2) +
  labs(
    title = "Closest BLAST genera by gene and sample type",
    subtitle = "Strict single-dominant sequences only; labels are closest known references, not confirmed species",
    x = "Gene",
    y = "Percent of strict consensus sequences",
    fill = "Closest BLAST genus"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_figs, "04_strict_blast_genus_composition_by_sample_type.png"), p_blast_type, width = 13, height = 9, dpi = 300)

p_blast_type_iupac <- blast_genus_plot_data %>%
  filter(tree_set == "iupac") %>%
  ggplot(aes(x = gene, y = percent, fill = closest_genus_report)) +
  geom_col(width = 0.8, color = "white", linewidth = 0.15) +
  facet_wrap(~ sample_type, ncol = 2) +
  labs(
    title = "Closest BLAST genera by gene and sample type",
    subtitle = "Mixed-IUPAC sequence set; use as sensitivity/exploratory",
    x = "Gene",
    y = "Percent of mixed-IUPAC consensus sequences",
    fill = "Closest BLAST genus"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_figs, "05_iupac_blast_genus_composition_by_sample_type.png"), p_blast_type_iupac, width = 13, height = 9, dpi = 300)

# ----------------------------
# 7. Metadata association screening
# ----------------------------

screen_association <- function(data, metadata_col, min_group_n = 10) {
  metadata_sym <- rlang::sym(metadata_col)

  d <- data %>%
    mutate(metadata_value = as.character(!!metadata_sym)) %>%
    mutate(metadata_value = ifelse(is.na(metadata_value) | metadata_value == "", "Unknown", metadata_value)) %>%
    filter(metadata_value != "Unknown") %>%
    group_by(metadata_value) %>%
    filter(n_distinct(sample) >= min_group_n) %>%
    ungroup()

  map_dfr(genes, function(g) {
    dg <- d %>% filter(gene == g)
    if (nrow(dg) == 0 || n_distinct(dg$metadata_value) < 2 || n_distinct(dg$good_coverage) < 2) {
      return(tibble(gene = g, metadata_category = metadata_col, test = NA_character_,
                    p_value = NA_real_, note = "not_enough_data"))
    }

    tab <- table(dg$metadata_value, dg$good_coverage)
    p <- tryCatch({
      if (any(tab < 5) || nrow(tab) > 2) {
        suppressWarnings(chisq.test(tab, simulate.p.value = TRUE, B = 10000)$p.value)
      } else {
        chisq.test(tab)$p.value
      }
    }, error = function(e) NA_real_)

    tibble(
      gene = g,
      metadata_category = metadata_col,
      test = "chi_square_or_simulated_chi_square",
      p_value = p,
      note = "screening_only_not_phylogenetic_clustering"
    )
  })
}

association_screen <- map_dfr(
  c("sample_type", "host_tribe", "host_genus", "site", "state", "native_status"),
  ~ screen_association(coverage_good, .x, min_group_n = 10)
) %>%
  mutate(
    p_adjust_BH_within_metadata = ave(p_value, metadata_category, FUN = function(x) p.adjust(x, method = "BH"))
  ) %>%
  arrange(metadata_category, p_value)

write_tsv(association_screen, file.path(out_tables, "05_metadata_gene_recovery_association_screen.tsv"))

p_assoc <- association_screen %>%
  filter(!is.na(p_value)) %>%
  mutate(
    gene = factor(gene, levels = genes),
    metadata_category = factor(metadata_category, levels = c("sample_type", "host_tribe", "host_genus", "site", "state", "native_status")),
    minus_log10_p = -log10(p_value)
  ) %>%
  ggplot(aes(x = gene, y = metadata_category, fill = minus_log10_p)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = ifelse(p_adjust_BH_within_metadata < 0.05, "*", "")), size = 5) +
  scale_fill_gradient(low = "#F7FBFF", high = "#CB181D", na.value = "gray90") +
  labs(
    title = "Screening test: gene recovery vs metadata",
    subtitle = "Asterisks mark BH-adjusted p < 0.05 within metadata category; this is not a phylogenetic clustering test",
    x = "Gene",
    y = "Metadata category",
    fill = "-log10(p)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_figs, "06_gene_recovery_metadata_association_screen.png"), p_assoc, width = 11, height = 6, dpi = 300)

# ----------------------------
# 8. Report-ready interpretation notes
# ----------------------------

report_notes <- tibble::tribble(
  ~section, ~key_point, ~report_text,
  "Gene categories", "Correct nod wording",
  "The selected genes include four nif genes, two canonical nod genes (nodL and nodX), and four accessory nodulation/symbiosis-related genes (nolG, nolF, noeA, and noeB). Therefore noe/nol genes should not be described simply as canonical nod genes.",
  "Tree interpretation", "Strict trees primary",
  "Strict single-dominant trees are the primary phylogenetic summaries because each included tip is a single dominant consensus sequence with <=20% N. Mixed-IUPAC trees include both strict and mixed possible multitemplate sequences and should be interpreted as sensitivity analyses, especially when bootstrap convergence is poor.",
  "BLAST interpretation", "Closest hit only",
  "BLAST taxon and genus labels represent the closest known reference sequence in the available functional-gene alignment. They are useful for biological orientation, but they are not confirmed species identifications.",
  "Ryan suggestion", "Functional gene vs 16S bridge",
  "The circle plots and associated tables can be used as a functional-gene counterpart to 16S summaries: for each sample or sample type, compare which functional genes are recovered and which closest BLAST genera dominate.",
  "Ahmed suggestion", "Cross-gene host/geography patterns",
  "The metadata recovery and BLAST summaries provide a first check of whether host or geographic patterns are consistent across nif, canonical nod, and accessory nodulation/symbiosis-related genes. Formal tree-based clustering tests can be added after these descriptive summaries."
)

write_tsv(report_notes, file.path(out_tables, "06_report_interpretation_notes.tsv"))

# ----------------------------
# 9. Output index
# ----------------------------

output_index <- tibble(
  output_type = c(rep("table", 6), rep("figure", 6)),
  file = c(
    file.path(out_tables, "01_gene_categories_for_report.tsv"),
    file.path(out_tables, "02_gene_recovery_by_metadata_category_pct80_depth10.tsv"),
    file.path(out_tables, "03_tree_reliability_summary_for_report.tsv"),
    file.path(out_tables, "04_blast_closest_genus_by_gene_sample_type.tsv"),
    file.path(out_tables, "05_metadata_gene_recovery_association_screen.tsv"),
    file.path(out_tables, "06_report_interpretation_notes.tsv"),
    file.path(out_figs, "01_gene_recovery_percent_by_sample_type.png"),
    file.path(out_figs, "02_gene_recovery_heatmap_host_geography.png"),
    file.path(out_figs, "03_tree_input_size_and_convergence.png"),
    file.path(out_figs, "04_strict_blast_genus_composition_by_sample_type.png"),
    file.path(out_figs, "05_iupac_blast_genus_composition_by_sample_type.png"),
    file.path(out_figs, "06_gene_recovery_metadata_association_screen.png")
  )
) %>%
  mutate(exists = file.exists(file))

write_tsv(output_index, file.path(out_tables, "00_metadata_tree_comparison_output_index.tsv"))

message("Done. Metadata/report comparison outputs saved here:")
message(out_dir)
print(output_index)

