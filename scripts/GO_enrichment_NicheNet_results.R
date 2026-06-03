# Gene Ontology enrichment analysis for the NicheNet results
#
# Run this script after running both NicheNet scripts:
# 1. NicheNet_assignment_Myotubes_SkMVECs.R
# 2. NicheNet_reverse_Myotubes_sender_SkMVECs_receiver.R

# The script performs GO Biological Process enrichment for:
# - receiver-upregulated genes
# - NicheNet-predicted target genes

# Install once if needed:
#install.packages("BiocManager")
#BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

### 1. Settings ---------------------------------------------------------------

project_dir <- "C:/Users/royde/OneDrive/Documenten/BIT11 Internship"

forward_results_dir <- file.path(
  project_dir,
  "repo_upload",
  "results",
  "nichenet_assignment_results"
)

reverse_results_dir <- file.path(
  project_dir,
  "repo_upload",
  "results",
  "nichenet_assignment_results_reverse"
)

go_output_dir <- file.path(
  project_dir,
  "results",
  "gene_ontology_enrichment"
)

dir.create(go_output_dir, recursive = TRUE, showWarnings = FALSE)

### 2. Read gene lists --------------------------------------------------------

# Forward analysis: SkMVECs --> Myotubes
forward_receiver_genes <- read_csv(
  file.path(forward_results_dir, "geneset_of_interest_receiver_upregulated.csv"),
  show_col_types = FALSE
) %>%
  pull(gene) %>%
  unique()

forward_target_genes <- read_csv(
  file.path(forward_results_dir, "active_ligand_target_links.csv"),
  show_col_types = FALSE
) %>%
  pull(target) %>%
  unique()

# Reverse analysis: Myotubes --> SkMVECs 
reverse_receiver_genes <- read_csv(
  file.path(reverse_results_dir, "geneset_of_interest_receiver_upregulated.csv"),
  show_col_types = FALSE
) %>%
  pull(gene) %>%
  unique()

reverse_target_genes <- read_csv(
  file.path(reverse_results_dir, "active_ligand_target_links.csv"),
  show_col_types = FALSE
) %>%
  pull(target) %>%
  unique()

### Read background gene lists -----------------------------------------------

# These background genes are all genes expressed in the receiver cell type.
# They make the GO enrichment more fair, because GO will comapre the gene list
# against genes that were actually detectable in the RNA-seq data.

forward_background_genes <- read_csv(
  file.path(forward_results_dir, "background_receiver_expressed_genes.csv"),
  show_col_types = FALSE
) %>%
  pull(gene) %>%
  unique()

reverse_background_genes <- read_csv(
  file.path(reverse_results_dir, "background_receiver_expressed_genes.csv"),
  show_col_types = FALSE
) %>%
  pull(gene) %>%
  unique()

### 3. Convert gene symbols to Entrez IDs -------------------------------------

# clusterProfiler works best with Entrez gene IDs. 
# Therefore, gene symbols are converted first.
forward_receiver_ids <- bitr(
  forward_receiver_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

forward_target_ids <- bitr(
  forward_target_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

reverse_receiver_ids <- bitr(
  reverse_receiver_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

reverse_target_ids <- bitr(
  reverse_target_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

# Convert background genes to Entrez IDs
forward_background_ids <- bitr(
  forward_background_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

reverse_background_ids <- bitr(
  reverse_background_genes,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

### 4. Run GO enrichment ------------------------------------------------------

# GO biological process enrichment for forward receiver-upregulated genes
go_forward_receiver <- enrichGO(
  gene = forward_receiver_ids$ENTREZID,
  universe = forward_background_ids$ENTREZID,
  OrgDb = org.Hs.eg.db, # database for Homo Sapiens
  keyType = "ENTREZID", # 
  ont = "BP", # Biological Process
  pAdjustMethod = "BH", # Benjamini-Hochberg: controls the false discovery rate: correct the p-values so the enrichment results are less likely to include false positives.
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
# GO Biological Process enrichment for forward NicheNet target genes.
go_forward_targets <- enrichGO(
  gene = forward_target_ids$ENTREZID,
  universe = forward_background_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
# GO Biological Process enrichment for reverse receiver-upregulated genes
go_reverse_receiver <- enrichGO(
  gene = reverse_receiver_ids$ENTREZID,
  universe = reverse_background_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)
# GO Biological Process enrichment for reverse NicheNet target genes
go_reverse_targets <- enrichGO(
  gene = reverse_target_ids$ENTREZID,
  universe = reverse_background_ids$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

### 5. Save GO result tables --------------------------------------------------

write_csv(
  as.data.frame(go_forward_receiver),
  file.path(go_output_dir, "SkMVECs_to_Myotubes_receiver_upregulated_GO_BP_results.csv")
)

write_csv(
  as.data.frame(go_forward_targets),
  file.path(go_output_dir, "SkMVECs_to_Myotubes_nichenet_targets_GO_BP_results.csv")
)

write_csv(
  as.data.frame(go_reverse_receiver),
  file.path(go_output_dir, "Myotubes_to_SkMVECs_receiver_upregulated_GO_BP_results.csv")
)

write_csv(
  as.data.frame(go_reverse_targets),
  file.path(go_output_dir, "Myotubes_to_SkMVECs_nichenet_targets_GO_BP_results.csv")
)

### 6. Make dotplots ----------------------------------------------------------
#make a dotplot from the GO enrichment result and show the top 15 enriched GO terms.
p1 <- dotplot(go_forward_receiver, showCategory = 15) +
  ggtitle("SkMVECs to Myotubes - receiver upregulated genes")

p2 <- dotplot(go_forward_targets, showCategory = 15) +
  ggtitle("SkMVECs to Myotubes - NicheNet target genes")

p3 <- dotplot(go_reverse_receiver, showCategory = 15) +
  ggtitle("Myotubes to SkMVECs - receiver upregulated genes")

p4 <- dotplot(go_reverse_targets, showCategory = 15) +
  ggtitle("Myotubes to SkMVECs - NicheNet target genes")

### 7. Save dotplots ----------------------------------------------------------

ggsave(
  file.path(go_output_dir, "SkMVECs_to_Myotubes_receiver_upregulated_GO_BP_dotplot.png"),
  p1,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(go_output_dir, "SkMVECs_to_Myotubes_nichenet_targets_GO_BP_dotplot.png"),
  p2,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(go_output_dir, "Myotubes_to_SkMVECs_receiver_upregulated_GO_BP_dotplot.png"),
  p3,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  file.path(go_output_dir, "Myotubes_to_SkMVECs_nichenet_targets_GO_BP_dotplot.png"),
  p4,
  width = 10,
  height = 7,
  dpi = 300
)

###############################################################################

run_go_enrichment <- function(gene_symbols, analysis_name, gene_set_name) {

  gene_symbols <- unique(gene_symbols)
  gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]

  if (length(gene_symbols) < 10) {
    warning("Too few genes for GO enrichment: ", analysis_name, " - ", gene_set_name)
    return(NULL)
  }

  # Convert gene symbols to Entrez IDs, because clusterProfiler uses Entrez IDs.
  gene_ids <- bitr(
    gene_symbols,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  if (nrow(gene_ids) < 10) {
    warning("Too few genes could be mapped to Entrez IDs: ", analysis_name, " - ", gene_set_name)
    return(NULL)
  }

  go_result <- enrichGO(
    gene = gene_ids$ENTREZID,
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = "BP",
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.05,
    readable = TRUE
  )

  result_table <- as.data.frame(go_result)

  output_prefix <- paste(analysis_name, gene_set_name, "GO_BP", sep = "_")

  write_csv(
    result_table,
    file.path(go_output_dir, paste0(output_prefix, "_results.csv"))
  )

  if (nrow(result_table) > 0) {
    dot_plot <- dotplot(go_result, showCategory = 15) +
      ggtitle(paste(analysis_name, "-", gene_set_name))

    ggsave(
      file.path(go_output_dir, paste0(output_prefix, "_dotplot.png")),
      dot_plot,
      width = 10,
      height = 7,
      dpi = 300
    )

    bar_plot <- barplot(go_result, showCategory = 15) +
      ggtitle(paste(analysis_name, "-", gene_set_name))

    ggsave(
      file.path(go_output_dir, paste0(output_prefix, "_barplot.png")),
      bar_plot,
      width = 10,
      height = 7,
      dpi = 300
    )
  }

  return(go_result)
}

### 3. Function to read NicheNet output and analyze it ------------------------

analyze_nichenet_folder <- function(results_dir, analysis_name) {

  geneset_file <- file.path(results_dir, "geneset_of_interest_receiver_upregulated.csv")
  target_file <- file.path(results_dir, "active_ligand_target_links.csv")

  if (!file.exists(geneset_file)) {
    warning("Missing file: ", geneset_file)
  } else {
    receiver_upregulated_genes <- read_csv(geneset_file, show_col_types = FALSE) %>%
      pull(gene)

    run_go_enrichment(
      receiver_upregulated_genes,
      analysis_name,
      "receiver_upregulated_genes"
    )
  }

  if (!file.exists(target_file)) {
    warning("Missing file: ", target_file)
  } else {
    nichenet_targets <- read_csv(target_file, show_col_types = FALSE) %>%
      pull(target)

    run_go_enrichment(
      nichenet_targets,
      analysis_name,
      "nichenet_predicted_targets"
    )
  }
}

### 4. Run GO enrichment for both NicheNet directions -------------------------

analyze_nichenet_folder(
  forward_results_dir,
  "SkMVECs_to_Myotubes"
)

analyze_nichenet_folder(
  reverse_results_dir,
  "Myotubes_to_SkMVECs"
)



