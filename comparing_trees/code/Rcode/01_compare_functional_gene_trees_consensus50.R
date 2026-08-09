# ============================================================
# 01_compare_functional_gene_trees_consensus50.R
#
# Purpose:
#   Compare strict single-dominant and mixed-IUPAC functional-gene trees
#   for the consensus_sequences_50percent analysis.
#
# Main questions addressed:
#   1. Which genes have strong sample coverage?
#   2. How many strict and mixed-IUPAC consensus sequences were used per gene?
#   3. Which trees converged in IQ-TREE?
#   4. How do nif and nod-related genes compare?
#   5. Which samples have good coverage across multiple genes?
#
# This script is local-only. It does not change tree files.
# It writes new summary tables and figures into:
#   result/comparative_tree_analysis/
# ============================================================

# ----------------------------
# 0. User paths
# ----------------------------

project_dir <- "/Users/rosa/Desktop/ALLWork/Madison/Project/ryan-nitfix/intbio-nitfix-git/symbiosis_sorted_all_sample_consensus_sequences_50percent"

result_dir <- file.path(project_dir, "result")
tables_dir <- file.path(result_dir, "tables")
gene_trees_dir <- file.path(result_dir, "gene_trees_full")

out_dir <- file.path(result_dir, "comparative_tree_analysis")
out_tables <- file.path(out_dir, "tables")
out_figs <- file.path(out_dir, "figures")

dir.create(out_tables, recursive = TRUE, showWarnings = FALSE)
dir.create(out_figs, recursive = TRUE, showWarnings = FALSE)

genes <- c("nifH", "nifD", "nifK", "nifJ", "nodL", "nolG", "nolF", "noeA", "noeB", "nodX")

nif_genes <- c("nifH", "nifD", "nifK", "nifJ")
nod_related_genes <- c("nodL", "nolG", "nolF", "noeA", "noeB", "nodX")

# ----------------------------
# 1. Packages
# ----------------------------

need_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste0("Please install this R package first: ", pkg,
                "\nExample: install.packages('", pkg, "')"), call. = FALSE)
  }
}

for (pkg in c("ggplot2", "dplyr", "tidyr", "readr", "stringr", "purrr")) {
  need_pkg(pkg)
}

library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(purrr)

theme_set(theme_bw(base_size = 13))

# ----------------------------
# 2. Helper functions
# ----------------------------

safe_file <- function(...) file.path(...)

first_existing <- function(paths) {
  hit <- paths[file.exists(paths)]
  if (length(hit) == 0) return(NA_character_)
  hit[1]
}

count_fasta <- function(path) {
  if (is.na(path) || !file.exists(path)) return(NA_integer_)
  sum(grepl("^>", readLines(path, warn = FALSE)))
}

fasta_lengths <- function(path) {
  if (is.na(path) || !file.exists(path)) return(integer())
  x <- readLines(path, warn = FALSE)
  starts <- which(grepl("^>", x))
  if (length(starts) == 0) return(integer())
  lens <- integer(length(starts))
  for (i in seq_along(starts)) {
    s <- starts[i] + 1
    e <- if (i < length(starts)) starts[i + 1] - 1 else length(x)
    if (s > e) {
      lens[i] <- 0
    } else {
      lens[i] <- nchar(paste0(x[s:e], collapse = ""))
    }
  }
  lens
}

first_fasta_length <- function(path) {
  lens <- fasta_lengths(path)
  if (length(lens) == 0) return(NA_integer_)
  lens[1]
}

sequence_composition <- function(path) {
  if (is.na(path) || !file.exists(path)) {
    return(tibble(
      ACGT_count = NA_real_, IUPAC_count_excluding_N = NA_real_, N_count = NA_real_, gap_count = NA_real_,
      total_characters = NA_real_, ACGT_percent = NA_real_, IUPAC_percent_excluding_N = NA_real_,
      N_percent = NA_real_, gap_percent = NA_real_
    ))
  }
  x <- readLines(path, warn = FALSE)
  seq <- toupper(paste0(x[!grepl("^>", x)], collapse = ""))
  chars <- strsplit(seq, "")[[1]]
  total <- length(chars)
  if (total == 0) total <- NA_real_

  acgt <- sum(chars %in% c("A", "C", "G", "T"), na.rm = TRUE)
  iupac <- sum(chars %in% c("R", "Y", "S", "W", "K", "M", "B", "D", "H", "V"), na.rm = TRUE)
  n <- sum(chars == "N", na.rm = TRUE)
  gap <- sum(chars == "-", na.rm = TRUE)

  tibble(
    ACGT_count = acgt,
    IUPAC_count_excluding_N = iupac,
    N_count = n,
    gap_count = gap,
    total_characters = total,
    ACGT_percent = 100 * acgt / total,
    IUPAC_percent_excluding_N = 100 * iupac / total,
    N_percent = 100 * n / total,
    gap_percent = 100 * gap / total
  )
}

read_log_lines <- function(log_file) {
  if (is.na(log_file) || !file.exists(log_file)) return(character())
  readLines(log_file, warn = FALSE)
}

get_model <- function(log_file) {
  lines <- read_log_lines(log_file)
  hit <- lines[grepl("Best-fit model", lines, ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_character_)
  line <- tail(hit, 1)
  line <- sub(".*Best-fit model: ", "", line)
  line <- sub(" chosen.*", "", line)
  line <- sub(" according.*", "", line)
  line
}

get_tree_length <- function(log_file) {
  lines <- read_log_lines(log_file)
  hit <- lines[grepl("Total tree length", lines, ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_real_)
  as.numeric(tail(strsplit(tail(hit, 1), "\\s+")[[1]], 1))
}

get_final_boot_corr <- function(log_file) {
  lines <- read_log_lines(log_file)
  hit <- lines[grepl("Bootstrap correlation coefficient", lines, ignore.case = TRUE)]
  if (length(hit) == 0) return(NA_real_)
  as.numeric(tail(strsplit(tail(hit, 1), "\\s+")[[1]], 1))
}

get_boot_status <- function(log_file) {
  lines <- read_log_lines(log_file)
  if (length(lines) == 0) return(NA_character_)
  if (any(grepl("WARNING: bootstrap analysis did not converge", lines, ignore.case = TRUE))) {
    "not_converged"
  } else {
    "converged"
  }
}

get_iqtree_taxa <- function(log_file) {
  lines <- read_log_lines(log_file)
  hit <- lines[grepl("^[0-9]+ taxa", lines)]
  if (length(hit) == 0) return(NA_integer_)
  as.integer(strsplit(tail(hit, 1), "\\s+")[[1]][1])
}

get_tree_paths <- function(gene, tree_set) {
  if (tree_set == "strict") {
    aln_candidates <- c(
      safe_file(gene_trees_dir, gene, "02_alignment", paste0(gene, "_consensus50_strict_single_dominant.mafft.fasta"))
    )
    prefix_candidates <- c(
      safe_file(gene_trees_dir, gene, "03_tree_strict_nm5000", paste0(gene, "_consensus50_strict_single_dominant_nm5000"))
    )
  } else {
    aln_candidates <- c(
      safe_file(gene_trees_dir, gene, "02_alignment", paste0(gene, "_consensus50_iupac_all_pass.mafft.fasta"))
    )
    prefix_candidates <- c(
      safe_file(gene_trees_dir, gene, "04_tree_iupac_nm5000", paste0(gene, "_consensus50_iupac_all_pass_nm5000")),
      safe_file(gene_trees_dir, paste0(gene, "_iupac_nm5000"), paste0(gene, "_consensus50_iupac_all_pass_nm5000"))
    )
  }

  prefix <- first_existing(paste0(prefix_candidates, ".treefile"))
  if (!is.na(prefix)) prefix <- sub("\\.treefile$", "", prefix)

  tibble(
    gene = gene,
    tree_set = tree_set,
    alignment_file = first_existing(aln_candidates),
    prefix = prefix,
    treefile = ifelse(is.na(prefix), NA_character_, paste0(prefix, ".treefile")),
    contree = ifelse(is.na(prefix), NA_character_, paste0(prefix, ".contree")),
    iqtree_report = ifelse(is.na(prefix), NA_character_, paste0(prefix, ".iqtree")),
    log = ifelse(is.na(prefix), NA_character_, paste0(prefix, ".log"))
  )
}

# ----------------------------
# 3. Coverage summaries
# ----------------------------

coverage_file <- file.path(tables_dir, "consensus50_gene_coverage_all_samples.tsv")
if (!file.exists(coverage_file)) {
  stop(paste0("Coverage table not found:\n", coverage_file), call. = FALSE)
}

coverage <- read_tsv(coverage_file, show_col_types = FALSE)

coverage_good <- coverage %>%
  filter(gene %in% genes) %>%
  mutate(
    pass80_depth10 = percent_covered >= 80 & mean_depth >= 10,
    sample_type = case_when(
      sample %in% c("MC-1", "MC-2") ~ "MC",
      str_detect(sample, "-No$") ~ "No",
      str_detect(sample, "-Rh$") ~ "Rh",
      str_detect(sample, "-Ro$") ~ "Ro",
      TRUE ~ "Unknown"
    ),
    gene_group = case_when(
      gene %in% nif_genes ~ "nif",
      gene %in% nod_related_genes ~ "nod/nol/noe",
      TRUE ~ "other"
    )
  )

coverage_gene_summary <- coverage_good %>%
  group_by(gene, gene_group) %>%
  summarise(
    total_samples = n_distinct(sample),
    good_samples_pct80_depth10 = n_distinct(sample[pass80_depth10]),
    good_sample_gene_rows = sum(pass80_depth10),
    percent_samples_good = 100 * good_samples_pct80_depth10 / total_samples,
    mean_percent_covered = mean(percent_covered, na.rm = TRUE),
    median_percent_covered = median(percent_covered, na.rm = TRUE),
    mean_depth = mean(mean_depth, na.rm = TRUE),
    median_depth = median(mean_depth, na.rm = TRUE),
    .groups = "drop"
  )

write_tsv(coverage_gene_summary, file.path(out_tables, "01_coverage_summary_by_gene_pct80_depth10.tsv"))

coverage_by_type <- coverage_good %>%
  group_by(gene, gene_group, sample_type) %>%
  summarise(
    total_samples = n_distinct(sample),
    good_samples = n_distinct(sample[pass80_depth10]),
    percent_good = 100 * good_samples / total_samples,
    .groups = "drop"
  )

write_tsv(coverage_by_type, file.path(out_tables, "02_coverage_by_gene_and_sample_type_pct80_depth10.tsv"))

# ----------------------------
# 4. Consensus QC summaries
# ----------------------------

qc_summary_long <- map_dfr(genes, function(g) {
  qc_file <- first_existing(c(
    file.path(gene_trees_dir, g, "01_consensus", paste0(g, "_consensus_qc.tsv")),
    file.path(project_dir, "gene_trees", g, "01_consensus", paste0(g, "_consensus_qc.tsv"))
  ))

  if (is.na(qc_file)) {
    return(tibble(gene = g, status = "missing_qc_file", count = NA_integer_))
  }

  read_tsv(qc_file, show_col_types = FALSE) %>%
    count(status, name = "count") %>%
    mutate(gene = g, .before = 1)
})

qc_summary_wide <- qc_summary_long %>%
  pivot_wider(names_from = status, values_from = count, values_fill = 0)

write_tsv(qc_summary_long, file.path(out_tables, "03_consensus_qc_status_summary_long.tsv"))
write_tsv(qc_summary_wide, file.path(out_tables, "04_consensus_qc_status_summary_wide.tsv"))

# ----------------------------
# 5. Tree summaries
# ----------------------------

tree_paths <- expand_grid(gene = genes, tree_set = c("strict", "iupac")) %>%
  pmap_dfr(function(gene, tree_set) get_tree_paths(gene, tree_set))

tree_summary <- tree_paths %>%
  rowwise() %>%
  mutate(
    input_alignment_sequences = count_fasta(alignment_file),
    alignment_length = first_fasta_length(alignment_file),
    iqtree_taxa_or_tips = get_iqtree_taxa(log),
    best_model = get_model(log),
    tree_length = get_tree_length(log),
    final_bootstrap_correlation = get_final_boot_corr(log),
    bootstrap_status = get_boot_status(log),
    treefile_exists = !is.na(treefile) && file.exists(treefile),
    contree_exists = !is.na(contree) && file.exists(contree),
    iqtree_report_exists = !is.na(iqtree_report) && file.exists(iqtree_report),
    log_exists = !is.na(log) && file.exists(log)
  ) %>%
  ungroup() %>%
  mutate(
    tree_set_label = recode(tree_set,
                            strict = "Strict single-dominant",
                            iupac = "Mixed-IUPAC"),
    gene_group = case_when(
      gene %in% nif_genes ~ "nif",
      gene %in% nod_related_genes ~ "nod/nol/noe",
      TRUE ~ "other"
    )
  )

write_tsv(tree_summary, file.path(out_tables, "05_tree_summary_strict_and_mixed_iupac.tsv"))
write_tsv(filter(tree_summary, tree_set == "strict"), file.path(out_tables, "06_strict_tree_summary.tsv"))
write_tsv(filter(tree_summary, tree_set == "iupac"), file.path(out_tables, "07_mixed_iupac_tree_summary.tsv"))

# ----------------------------
# 6. Alignment composition summaries
# ----------------------------

composition_summary <- tree_paths %>%
  rowwise() %>%
  mutate(comp = list(sequence_composition(alignment_file))) %>%
  unnest(comp) %>%
  ungroup() %>%
  mutate(
    tree_set_label = recode(tree_set,
                            strict = "Strict single-dominant",
                            iupac = "Mixed-IUPAC")
  )

write_tsv(composition_summary, file.path(out_tables, "08_alignment_sequence_composition.tsv"))

# ----------------------------
# 7. Sample overlap across genes
# ----------------------------

sample_gene_matrix <- coverage_good %>%
  group_by(sample, sample_type, gene) %>%
  summarise(pass80_depth10 = any(pass80_depth10), .groups = "drop") %>%
  pivot_wider(names_from = gene, values_from = pass80_depth10, values_fill = FALSE) %>%
  mutate(
    n_good_genes = rowSums(across(all_of(genes))),
    n_good_nif_genes = rowSums(across(all_of(nif_genes))),
    n_good_nod_related_genes = rowSums(across(all_of(nod_related_genes)))
  ) %>%
  select(sample, sample_type, n_good_genes, n_good_nif_genes, n_good_nod_related_genes, all_of(genes))

write_tsv(sample_gene_matrix, file.path(out_tables, "09_sample_overlap_good_genes_pct80_depth10.tsv"))

sample_overlap_summary <- sample_gene_matrix %>%
  count(n_good_genes, sample_type, name = "n_samples") %>%
  arrange(sample_type, n_good_genes)

write_tsv(sample_overlap_summary, file.path(out_tables, "10_sample_overlap_summary_by_sample_type.tsv"))

sample_gene_combinations <- sample_gene_matrix %>%
  rowwise() %>%
  mutate(
    recovered_genes = {
      vals <- c_across(all_of(genes))
      hits <- genes[as.logical(vals)]
      if (length(hits) == 0) "None" else paste(hits, collapse = "+")
    }
  ) %>%
  ungroup()

write_tsv(sample_gene_combinations, file.path(out_tables, "10a_sample_gene_recovery_combinations_pct80_depth10.tsv"))

gene_combination_summary <- sample_gene_combinations %>%
  count(recovered_genes, n_good_genes, name = "total_samples") %>%
  arrange(desc(total_samples), desc(n_good_genes), recovered_genes)

write_tsv(gene_combination_summary, file.path(out_tables, "10b_top_gene_recovery_combinations_pct80_depth10.tsv"))

gene_combination_by_type <- sample_gene_combinations %>%
  count(recovered_genes, n_good_genes, sample_type, name = "n_samples") %>%
  left_join(gene_combination_summary %>% select(recovered_genes, total_samples), by = "recovered_genes") %>%
  arrange(desc(total_samples), recovered_genes, sample_type)

write_tsv(gene_combination_by_type, file.path(out_tables, "10c_gene_recovery_combinations_by_sample_type_pct80_depth10.tsv"))

# ----------------------------
# 8. Figures
# ----------------------------

p_coverage <- coverage_gene_summary %>%
  ggplot(aes(x = reorder(gene, -good_samples_pct80_depth10), y = good_samples_pct80_depth10, fill = gene_group)) +
  geom_col(width = 0.75) +
  geom_text(aes(label = good_samples_pct80_depth10), vjust = -0.3, size = 3.2) +
  scale_fill_manual(values = c(nif = "#2C7FB8", `nod/nol/noe` = "#7B52AB", other = "gray60")) +
  labs(
    title = "Samples with good coverage by functional gene",
    subtitle = "Good coverage: percent covered >= 80% and mean depth >= 10",
    x = "Gene",
    y = "Number of samples",
    fill = "Gene group"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_figs, "01_good_coverage_samples_by_gene.png"), p_coverage, width = 10, height = 6, dpi = 300)

p_type_heatmap <- coverage_by_type %>%
  ggplot(aes(x = gene, y = sample_type, fill = percent_good)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.1f%%", percent_good)), size = 3) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C", na.value = "gray90") +
  labs(
    title = "Good coverage by gene and sample type",
    subtitle = "Good coverage: percent covered >= 80% and mean depth >= 10",
    x = "Gene",
    y = "Sample type",
    fill = "% good"
  )

ggsave(file.path(out_figs, "02_coverage_heatmap_gene_by_sample_type.png"), p_type_heatmap, width = 10, height = 5, dpi = 300)

p_tree_counts <- tree_summary %>%
  ggplot(aes(x = gene, y = input_alignment_sequences, fill = tree_set_label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  geom_text(aes(label = input_alignment_sequences),
            position = position_dodge(width = 0.8), vjust = -0.25, size = 3) +
  scale_fill_manual(values = c("Strict single-dominant" = "#2C7FB8", "Mixed-IUPAC" = "#F28E2B")) +
  labs(
    title = "Strict vs mixed-IUPAC tree input sequences",
    subtitle = "Strict = single-dominant only; mixed-IUPAC = strict single-dominant + mixed possible multitemplate",
    x = "Gene",
    y = "Number of sample-derived consensus sequences / tree tips",
    fill = "Tree sequence set"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_figs, "03_strict_vs_mixed_tree_input_sequences.png"), p_tree_counts, width = 11, height = 6, dpi = 300)

p_boot <- tree_summary %>%
  ggplot(aes(x = gene, y = final_bootstrap_correlation, color = bootstrap_status, shape = tree_set_label)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0.99, linetype = "dashed", color = "gray40") +
  scale_color_manual(values = c(converged = "#1B9E77", not_converged = "#D95F02")) +
  labs(
    title = "Bootstrap convergence by gene and tree set",
    subtitle = "Dashed line marks the approximate 0.99 convergence target",
    x = "Gene",
    y = "Final bootstrap correlation",
    color = "Bootstrap status",
    shape = "Tree set"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_figs, "04_bootstrap_convergence_by_gene_tree_set.png"), p_boot, width = 11, height = 6, dpi = 300)

top_n_combinations <- 20

top_combos <- gene_combination_summary %>%
  slice_max(total_samples, n = top_n_combinations, with_ties = FALSE) %>%
  pull(recovered_genes)

combo_plot_data <- gene_combination_by_type %>%
  filter(recovered_genes %in% top_combos) %>%
  mutate(
    combo_label = paste0(recovered_genes, " (n=", total_samples, ")"),
    combo_label = factor(combo_label, levels = rev(unique(combo_label[match(top_combos, recovered_genes)])))
  )

combo_label_levels <- gene_combination_summary %>%
  filter(recovered_genes %in% top_combos) %>%
  mutate(combo_label = paste0(recovered_genes, " (n=", total_samples, ")")) %>%
  pull(combo_label)

combo_plot_data <- combo_plot_data %>%
  mutate(combo_label = factor(combo_label, levels = rev(combo_label_levels)))

p_gene_combinations <- combo_plot_data %>%
  ggplot(aes(x = combo_label, y = n_samples, fill = sample_type)) +
  geom_col(color = "white", width = 0.75) +
  coord_flip() +
  scale_fill_manual(values = c(MC = "#E41A1C", No = "#377EB8", Rh = "#FF9900", Ro = "#7B52AB", Unknown = "gray70")) +
  labs(
    title = "Most common recovered functional-gene combinations",
    subtitle = "Genes counted when percent covered >= 80% and mean depth >= 10; labels show total samples",
    x = "Recovered gene combination",
    y = "Number of samples",
    fill = "Sample type"
  )

ggsave(file.path(out_figs, "05_top_gene_recovery_combinations.png"), p_gene_combinations, width = 12, height = 8, dpi = 300)

combo_heatmap_data <- expand_grid(
  recovered_genes = top_combos,
  gene = genes
) %>%
  mutate(
    recovered = recovered_genes != "None" & str_detect(recovered_genes, paste0("(^|\\+)", gene, "(\\+|$)"))
  ) %>%
  left_join(
    gene_combination_summary %>%
      filter(recovered_genes %in% top_combos) %>%
      mutate(combo_label = paste0(recovered_genes, " (n=", total_samples, ")")) %>%
      select(recovered_genes, combo_label, total_samples, n_good_genes),
    by = "recovered_genes"
  ) %>%
  mutate(
    combo_label = factor(combo_label, levels = rev(combo_label_levels)),
    gene = factor(gene, levels = genes)
  )

p_combo_heatmap <- combo_heatmap_data %>%
  ggplot(aes(x = gene, y = combo_label, fill = recovered)) +
  geom_tile(color = "white", linewidth = 0.35) +
  scale_fill_manual(values = c(`TRUE` = "#1B9E77", `FALSE` = "gray92"), labels = c(`TRUE` = "Recovered", `FALSE` = "Not recovered")) +
  labs(
    title = "Which genes are recovered together?",
    subtitle = "Top recovered gene combinations across samples at percent covered >= 80% and mean depth >= 10",
    x = "Gene",
    y = "Recovered gene combination",
    fill = "Gene status"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_figs, "05b_gene_recovery_combination_heatmap.png"), p_combo_heatmap, width = 13, height = 8, dpi = 300)

composition_long <- composition_summary %>%
  select(gene, tree_set_label, ACGT_percent, IUPAC_percent_excluding_N, N_percent, gap_percent) %>%
  pivot_longer(cols = c(ACGT_percent, IUPAC_percent_excluding_N, N_percent, gap_percent),
               names_to = "category", values_to = "percent") %>%
  mutate(category = recode(category,
                           ACGT_percent = "A/C/G/T",
                           IUPAC_percent_excluding_N = "IUPAC excluding N",
                           N_percent = "N",
                           gap_percent = "Gap"))

p_comp <- composition_long %>%
  ggplot(aes(x = gene, y = percent, fill = category)) +
  geom_col(width = 0.75) +
  facet_wrap(~ tree_set_label, ncol = 1) +
  scale_fill_manual(values = c("A/C/G/T" = "#1B9E77", "IUPAC excluding N" = "#7570B3", "N" = "#D95F02", "Gap" = "gray70")) +
  labs(
    title = "Alignment composition by gene and tree set",
    x = "Gene",
    y = "Percent of aligned characters",
    fill = "Character class"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(out_figs, "06_alignment_composition_by_gene_tree_set.png"), p_comp, width = 11, height = 8, dpi = 300)

# ----------------------------
# 9. BLAST summaries if available locally
# ----------------------------

# The BLAST/taxon files may live either in result/tables or inside the
# per-gene tree folders. This section collects every available file, labels each
# row by tree set, removes duplicate files for the same gene/tree/sample, and
# plots strict and mixed-IUPAC assignments separately. This avoids adding strict,
# mixed, and relaxed nifH annotations together in one bar.
blast_candidates <- unique(c(
  list.files(
    tables_dir,
    pattern = "best_reference_hit_with_taxon\\.tsv$",
    full.names = TRUE,
    recursive = FALSE
  ),
  list.files(
    gene_trees_dir,
    pattern = "best_reference_hit_with_taxon\\.tsv$",
    full.names = TRUE,
    recursive = TRUE
  )
))

blast_min_samples <- 5

if (length(blast_candidates) > 0) {
  blast_long_all <- map_dfr(blast_candidates, function(f) {
    x <- read_tsv(f, show_col_types = FALSE)
    if (!all(c("sample_id", "best_reference_taxon") %in% names(x))) return(tibble())

    file_name <- basename(f)
    guessed_gene <- str_extract(file_name, paste(genes, collapse = "|"))
    if (is.na(guessed_gene)) {
      guessed_gene <- str_extract(f, paste(genes, collapse = "|"))
    }

    tree_set <- case_when(
      str_detect(file_name, "strict_pct60") ~ "relaxed_strict_pct60",
      str_detect(file_name, "iupac") ~ "iupac",
      str_detect(file_name, "strict|sample_") ~ "strict",
      TRUE ~ "unknown"
    )

    x %>%
      mutate(
        source_file = f,
        source_file_name = file_name,
        tree_set = tree_set,
        gene = ifelse("gene" %in% names(.), gene, guessed_gene),
        closest_taxon = best_reference_taxon,
        closest_genus = str_extract(best_reference_taxon, "^[A-Za-z]+"),
        closest_genus = ifelse(is.na(closest_genus) | closest_genus == "", "Unassigned", closest_genus)
      ) %>%
      filter(gene %in% genes) %>%
      select(source_file, source_file_name, tree_set, gene, sample_id, closest_taxon, closest_genus, everything())
  })

  if (nrow(blast_long_all) > 0) {
    blast_source_counts <- blast_long_all %>%
      count(gene, tree_set, source_file_name, name = "rows") %>%
      arrange(gene, tree_set, source_file_name)

    write_tsv(blast_source_counts, file.path(out_tables, "11a_blast_source_file_row_counts.tsv"))

    blast_long <- blast_long_all %>%
      filter(tree_set %in% c("strict", "iupac")) %>%
      arrange(gene, tree_set, sample_id, desc(bitscore), desc(percent_identity), desc(query_coverage_percent)) %>%
      distinct(gene, tree_set, sample_id, .keep_all = TRUE) %>%
      mutate(
        tree_set_label = recode(tree_set,
                                strict = "Strict single-dominant",
                                iupac = "Mixed-IUPAC")
      )

    write_tsv(blast_long, file.path(out_tables, "11b_blast_assignment_rows_deduplicated_strict_iupac.tsv"))

    blast_summary_raw <- blast_long %>%
      count(tree_set, tree_set_label, gene, closest_genus, name = "n") %>%
      arrange(tree_set, gene, desc(n), closest_genus)

    write_tsv(blast_summary_raw, file.path(out_tables, "11_blast_closest_genus_summary_available_files.tsv"))

    blast_summary_grouped <- blast_summary_raw %>%
      mutate(
        closest_genus_report = ifelse(n >= blast_min_samples, closest_genus, paste0("Other (<", blast_min_samples, ")"))
      ) %>%
      group_by(tree_set, tree_set_label, gene, closest_genus_report) %>%
      summarise(n = sum(n), .groups = "drop") %>%
      group_by(tree_set, tree_set_label, gene) %>%
      mutate(
        gene_total = sum(n),
        percent = 100 * n / gene_total
      ) %>%
      ungroup() %>%
      arrange(tree_set, gene, desc(n), closest_genus_report)

    write_tsv(blast_summary_grouped, file.path(out_tables, "12_blast_closest_genus_summary_grouped_min5.tsv"))

    p_blast <- blast_summary_grouped %>%
      mutate(
        gene = factor(gene, levels = genes),
        tree_set_label = factor(tree_set_label, levels = c("Strict single-dominant", "Mixed-IUPAC"))
      ) %>%
      ggplot(aes(x = gene, y = n, fill = closest_genus_report)) +
      geom_col(width = 0.8, color = "white", linewidth = 0.15) +
      facet_wrap(~ tree_set_label, ncol = 1, scales = "free_y") +
      labs(
        title = "Closest BLAST genera by functional gene and tree set",
        subtitle = paste0("Duplicate BLAST files are removed by gene/tree/sample; genera with fewer than ",
                          blast_min_samples, " assignments per gene are grouped as Other"),
        x = "Functional gene",
        y = "Number of sample consensus sequences",
        fill = "Closest BLAST genus"
      ) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right"
      )

    ggsave(file.path(out_figs, "07_top_closest_blast_genera_by_gene.png"), p_blast, width = 13, height = 10, dpi = 300)

    p_blast_heatmap <- blast_summary_grouped %>%
      filter(closest_genus_report != paste0("Other (<", blast_min_samples, ")")) %>%
      mutate(
        gene = factor(gene, levels = genes),
        tree_set_label = factor(tree_set_label, levels = c("Strict single-dominant", "Mixed-IUPAC"))
      ) %>%
      ggplot(aes(x = gene, y = closest_genus_report, fill = n)) +
      geom_tile(color = "white", linewidth = 0.3) +
      geom_text(aes(label = n), size = 3) +
      facet_wrap(~ tree_set_label, ncol = 1) +
      scale_fill_gradient(low = "#f7fbff", high = "#08519c") +
      labs(
        title = "Gene-by-genus BLAST assignment summary",
        subtitle = paste0("Only gene-genus combinations with at least ", blast_min_samples,
                          " sample consensus sequences; strict and mixed-IUPAC are separated"),
        x = "Functional gene",
        y = "Closest BLAST genus",
        fill = "Samples"
      ) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(out_figs, "08_blast_genus_by_gene_heatmap_min5.png"), p_blast_heatmap, width = 12, height = 10, dpi = 300)
  }
}

# ----------------------------
# 10. Output index
# ----------------------------

output_index <- tibble(
  output_type = c(rep("table", 17), rep("figure", 9)),
  file = c(
    file.path(out_tables, "01_coverage_summary_by_gene_pct80_depth10.tsv"),
    file.path(out_tables, "02_coverage_by_gene_and_sample_type_pct80_depth10.tsv"),
    file.path(out_tables, "03_consensus_qc_status_summary_long.tsv"),
    file.path(out_tables, "04_consensus_qc_status_summary_wide.tsv"),
    file.path(out_tables, "05_tree_summary_strict_and_mixed_iupac.tsv"),
    file.path(out_tables, "06_strict_tree_summary.tsv"),
    file.path(out_tables, "07_mixed_iupac_tree_summary.tsv"),
    file.path(out_tables, "08_alignment_sequence_composition.tsv"),
    file.path(out_tables, "09_sample_overlap_good_genes_pct80_depth10.tsv"),
    file.path(out_tables, "10_sample_overlap_summary_by_sample_type.tsv"),
    file.path(out_tables, "10a_sample_gene_recovery_combinations_pct80_depth10.tsv"),
    file.path(out_tables, "10b_top_gene_recovery_combinations_pct80_depth10.tsv"),
    file.path(out_tables, "10c_gene_recovery_combinations_by_sample_type_pct80_depth10.tsv"),
    file.path(out_tables, "11_blast_closest_genus_summary_available_files.tsv"),
    file.path(out_tables, "11a_blast_source_file_row_counts.tsv"),
    file.path(out_tables, "11b_blast_assignment_rows_deduplicated_strict_iupac.tsv"),
    file.path(out_tables, "12_blast_closest_genus_summary_grouped_min5.tsv"),
    file.path(out_figs, "01_good_coverage_samples_by_gene.png"),
    file.path(out_figs, "02_coverage_heatmap_gene_by_sample_type.png"),
    file.path(out_figs, "03_strict_vs_mixed_tree_input_sequences.png"),
    file.path(out_figs, "04_bootstrap_convergence_by_gene_tree_set.png"),
    file.path(out_figs, "05_top_gene_recovery_combinations.png"),
    file.path(out_figs, "05b_gene_recovery_combination_heatmap.png"),
    file.path(out_figs, "06_alignment_composition_by_gene_tree_set.png"),
    file.path(out_figs, "07_top_closest_blast_genera_by_gene.png"),
    file.path(out_figs, "08_blast_genus_by_gene_heatmap_min5.png")
  )
) %>%
  mutate(exists = file.exists(file))

write_tsv(output_index, file.path(out_tables, "00_comparative_analysis_output_index.tsv"))

message("Done. Comparative analysis outputs saved here:")
message(out_dir)
print(output_index)
