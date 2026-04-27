# NicheNet ligand activity analysis for the internship assignment
# Comparison 4: matched donors, Myotubes vs SkMVECs 2D confluent
#
# This script reads the provided raw gene count files, performs paired
# differential expression with edgeR, and runs NicheNet on the genes
# upregulated in the receiver cell type.

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(nichenetr)
  library(cowplot)
  library(ggpubr)
})

### 0. Settings ---------------------------------------------------------------

# Main biological question:
# Which ligands expressed by SkMVECs could explain the transcriptional program
# observed in Myotubes?
#
# To run the opposite direction, change these two values:
sender_celltype <- "SkMVECs"
receiver_celltype <- "Myotubes"

organism <- "human"
top_n_ligands <- 30
min_cpm <- 1
min_samples_expressed <- 2
fdr_cutoff <- 0.05
logfc_cutoff <- 1

get_script_dir <- function() {
  file_arg <- "--file="
  cmd_args <- commandArgs(trailingOnly = FALSE)
  script_path <- sub(file_arg, "", cmd_args[grepl(file_arg, cmd_args)])

  if (length(script_path) > 0) {
    return(dirname(normalizePath(script_path)))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    active_path <- rstudioapi::getActiveDocumentContext()$path
    if (!is.null(active_path) && nzchar(active_path)) {
      return(dirname(normalizePath(active_path)))
    }
  }

  getwd()
}

script_dir <- get_script_dir()
project_dir <- if (basename(script_dir) == "scripts") {
  normalizePath(file.path(script_dir, ".."), mustWork = FALSE)
} else {
  script_dir
}

count_dir <- file.path(project_dir, "data", "counts")
if (!dir.exists(count_dir)) {
  count_dir <- project_dir
}

output_dir <- file.path(project_dir, "results", "nichenet_assignment_results")
network_dir <- file.path(project_dir, "nichenet_networks")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(network_dir, showWarnings = FALSE, recursive = TRUE)

sample_info <- tribble(
  ~sample_id, ~donor, ~cell_type, ~file_name,
  "GC122806", "788", "Myotubes", "GC122806_KG.count",
  "GC122807", "790", "Myotubes", "GC122807_KG.count",
  "GC122808", "792", "Myotubes", "GC122808_KG.count",
  "GC122809", "794", "Myotubes", "GC122809_KG.count",
  "GC122810", "796", "Myotubes", "GC122810_KG.count",
  "GC122786", "788", "SkMVECs", "GC122786_KG.count",
  "GC122787", "790", "SkMVECs", "GC122787_KG.count",
  "GC122788", "792", "SkMVECs", "GC122788_KG.count",
  "GC124881", "794", "SkMVECs", "GC124881_KG.count",
  "GC122790", "796", "SkMVECs", "GC122790_KG.count"
) %>%
  mutate(
    file_path = file.path(count_dir, file_name),
    donor = factor(donor),
    cell_type = relevel(factor(cell_type), ref = sender_celltype)
  )

if (!all(file.exists(sample_info$file_path))) {
  stop(
    "Missing count file(s): ",
    paste(sample_info$file_path[!file.exists(sample_info$file_path)], collapse = ", "),
    call. = FALSE
  )
}

### 1. Read count files -------------------------------------------------------

read_count_file <- function(path, sample_id) {
  readr::read_tsv(
    path,
    col_names = c("gene", sample_id),
    col_types = cols(gene = col_character(), .default = col_double())
  ) %>%
    filter(!is.na(gene), gene != "") %>%
    group_by(gene) %>%
    summarise("{sample_id}" := sum(.data[[sample_id]], na.rm = TRUE), .groups = "drop")
}

count_tables <- map2(sample_info$file_path, sample_info$sample_id, read_count_file)

counts_tbl <- reduce(count_tables, full_join, by = "gene") %>%
  mutate(across(-gene, ~replace_na(.x, 0))) %>%
  arrange(gene)

counts <- counts_tbl %>%
  column_to_rownames("gene") %>%
  select(all_of(sample_info$sample_id)) %>%
  as.matrix()

storage.mode(counts) <- "integer"

write_csv(counts_tbl, file.path(output_dir, "combined_counts_matrix.csv"))
write_csv(sample_info, file.path(output_dir, "sample_metadata.csv"))

### 2. Paired differential expression ----------------------------------------

y <- DGEList(counts = counts, samples = sample_info)
keep <- filterByExpr(y, group = sample_info$cell_type)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

design <- model.matrix(~ donor + cell_type, data = sample_info)
colnames(design) <- make.names(colnames(design))

receiver_coef <- make.names(paste0("cell_type", receiver_celltype))
if (!receiver_coef %in% colnames(design)) {
  stop(
    "Could not find receiver coefficient in design matrix: ",
    receiver_coef,
    "\nDesign columns are: ",
    paste(colnames(design), collapse = ", "),
    call. = FALSE
  )
}

y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef = receiver_coef)

de_results <- topTags(qlf, n = Inf)$table %>%
  rownames_to_column("gene") %>%
  as_tibble() %>%
  arrange(FDR)

write_csv(de_results, file.path(output_dir, "DE_receiver_vs_sender_edgeR.csv"))

geneset_oi_raw <- de_results %>%
  filter(FDR <= fdr_cutoff, logFC >= logfc_cutoff) %>%
  pull(gene)

if (length(geneset_oi_raw) < 20) {
  warning(
    "The strict gene set has fewer than 20 genes. ",
    "For exploratory analysis, consider fdr_cutoff <- 0.10 or logfc_cutoff <- 0.5."
  )
}

### 3. Load NicheNet prior model ---------------------------------------------

download_if_missing <- function(url, destfile) {
  if (!file.exists(destfile)) {
    message("Downloading: ", basename(destfile))
    options(timeout = 1200)
    download.file(url, destfile = destfile, mode = "wb")
  }
  destfile
}

if (organism == "human") {
  lr_network_file <- download_if_missing(
    "https://zenodo.org/record/7074291/files/lr_network_human_21122021.rds",
    file.path(network_dir, "lr_network_human_21122021.rds")
  )
  ligand_target_file <- download_if_missing(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final.rds",
    file.path(network_dir, "ligand_target_matrix_nsga2r_final.rds")
  )
  weighted_networks_file <- download_if_missing(
    "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final.rds",
    file.path(network_dir, "weighted_networks_nsga2r_final.rds")
  )
} else if (organism == "mouse") {
  lr_network_file <- download_if_missing(
    "https://zenodo.org/record/7074291/files/lr_network_mouse_21122021.rds",
    file.path(network_dir, "lr_network_mouse_21122021.rds")
  )
  ligand_target_file <- download_if_missing(
    "https://zenodo.org/record/7074291/files/ligand_target_matrix_nsga2r_final_mouse.rds",
    file.path(network_dir, "ligand_target_matrix_nsga2r_final_mouse.rds")
  )
  weighted_networks_file <- download_if_missing(
    "https://zenodo.org/record/7074291/files/weighted_networks_nsga2r_final_mouse.rds",
    file.path(network_dir, "weighted_networks_nsga2r_final_mouse.rds")
  )
} else {
  stop("organism must be 'human' or 'mouse'.", call. = FALSE)
}

lr_network <- readRDS(lr_network_file) %>% distinct(from, to)
ligand_target_matrix <- readRDS(ligand_target_file)
weighted_networks <- readRDS(weighted_networks_file)

### 4. Define sender ligands and receiver receptors --------------------------

plain_cpm <- cpm(y)

get_expressed_genes <- function(cell_type) {
  sample_ids <- sample_info %>%
    filter(.data$cell_type == .env$cell_type) %>%
    pull(sample_id)

  rownames(plain_cpm)[
    rowSums(plain_cpm[, sample_ids, drop = FALSE] >= min_cpm) >= min_samples_expressed
  ]
}

expressed_genes_sender <- get_expressed_genes(sender_celltype)
expressed_genes_receiver <- get_expressed_genes(receiver_celltype)

ligands <- lr_network %>% pull(from) %>% unique()
receptors <- lr_network %>% pull(to) %>% unique()

expressed_ligands <- intersect(ligands, expressed_genes_sender)
expressed_receptors <- intersect(receptors, expressed_genes_receiver)

potential_ligands <- lr_network %>%
  filter(from %in% expressed_ligands, to %in% expressed_receptors) %>%
  pull(from) %>%
  unique() %>%
  intersect(colnames(ligand_target_matrix))

geneset_oi <- intersect(geneset_oi_raw, rownames(ligand_target_matrix))
background_expressed_genes <- intersect(expressed_genes_receiver, rownames(ligand_target_matrix))

write_csv(tibble(gene = geneset_oi), file.path(output_dir, "geneset_of_interest_receiver_upregulated.csv"))
write_csv(tibble(gene = background_expressed_genes), file.path(output_dir, "background_receiver_expressed_genes.csv"))
write_csv(tibble(ligand = potential_ligands), file.path(output_dir, "potential_ligands.csv"))

message("Sender cell type: ", sender_celltype)
message("Receiver cell type: ", receiver_celltype)
message("Expressed sender ligands: ", length(expressed_ligands))
message("Expressed receiver receptors: ", length(expressed_receptors))
message("Potential ligands: ", length(potential_ligands))
message("Geneset of interest: ", length(geneset_oi))
message("Background genes: ", length(background_expressed_genes))

if (length(geneset_oi) < 10) {
  stop(
    "Too few receiver-upregulated genes overlap the NicheNet model. ",
    "Try relaxing fdr_cutoff/logfc_cutoff or check gene symbols.",
    call. = FALSE
  )
}

### 5. NicheNet ligand activity analysis -------------------------------------

ligand_activities <- predict_ligand_activities(
  geneset = geneset_oi,
  background_expressed_genes = background_expressed_genes,
  ligand_target_matrix = ligand_target_matrix,
  potential_ligands = potential_ligands
) %>%
  arrange(desc(aupr_corrected)) %>%
  mutate(rank = row_number())

write_csv(ligand_activities, file.path(output_dir, "ligand_activities.csv"))

best_upstream_ligands <- ligand_activities %>%
  slice_head(n = top_n_ligands) %>%
  pull(test_ligand)

write_lines(best_upstream_ligands, file.path(output_dir, "best_upstream_ligands.txt"))

### 6. Infer ligand-target and ligand-receptor links -------------------------

active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(
    get_weighted_ligand_target_links,
    geneset = geneset_oi,
    ligand_target_matrix = ligand_target_matrix,
    n = 200
  ) %>%
  bind_rows()

write_csv(active_ligand_target_links_df, file.path(output_dir, "active_ligand_target_links.csv"))

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25
)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>%
  unique() %>%
  intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets, order_ligands])

p_ligand_target <- make_heatmap_ggplot(
  vis_ligand_target,
  y_name = paste("Prioritized", sender_celltype, "ligands"),
  x_name = paste(receiver_celltype, "upregulated genes"),
  color = "purple",
  legend_title = "Regulatory potential"
) +
  scale_fill_gradient2(low = "whitesmoke", high = "purple")

ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands,
  expressed_receptors,
  lr_network,
  weighted_networks$lr_sig
)

write_csv(ligand_receptor_links_df, file.path(output_dir, "active_ligand_receptor_links.csv"))

vis_ligand_receptor <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both"
)

p_ligand_receptor <- make_heatmap_ggplot(
  t(vis_ligand_receptor),
  y_name = paste("Prioritized", sender_celltype, "ligands"),
  x_name = paste("Receptors expressed by", receiver_celltype),
  color = "mediumvioletred",
  legend_title = "Prior interaction potential"
)

### 7. Summary figures --------------------------------------------------------

vis_ligand_aupr <- ligand_activities %>%
  filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>%
  select(aupr_corrected) %>%
  arrange(aupr_corrected) %>%
  as.matrix(ncol = 1)

p_ligand_aupr <- make_heatmap_ggplot(
  vis_ligand_aupr,
  y_name = paste("Prioritized", sender_celltype, "ligands"),
  x_name = "Ligand activity",
  color = "darkorange",
  legend_title = "AUPR"
) +
  theme(axis.text.x.top = element_blank())

ggsave(
  file.path(output_dir, "01_ligand_activity_aupr.png"),
  p_ligand_aupr,
  width = 5,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(output_dir, "02_ligand_target_heatmap.png"),
  p_ligand_target,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(output_dir, "03_ligand_receptor_heatmap.png"),
  p_ligand_receptor,
  width = 10,
  height = 7,
  dpi = 300
)

combined_without_legends <- plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_ligand_target + theme(
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.title.y = element_blank()
  ),
  p_ligand_receptor + theme(legend.position = "none"),
  ncol = 1,
  align = "v",
  rel_heights = c(1, 2, 2)
)

final_figure <- plot_grid(
  combined_without_legends,
  get_legend(p_ligand_aupr),
  ncol = 2,
  rel_widths = c(0.9, 0.1)
)

ggsave(
  file.path(output_dir, "04_nichenet_summary_figure.png"),
  final_figure,
  width = 12,
  height = 13,
  dpi = 300
)

message("Analysis finished. Results are in: ", normalizePath(output_dir))
