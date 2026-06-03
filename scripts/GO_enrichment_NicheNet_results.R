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

### Save gene ID mapping tables -----------------------------------

# These files show which gene symbols were successfully converted to Entrez IDs.
# This is useful because genes that cannot be mapped are not included in GO enrichment.

write_csv(
  forward_receiver_ids,
  file.path(go_output_dir, "mapped_forward_receiver_genes.csv")
)

write_csv(
  forward_target_ids,
  file.path(go_output_dir, "mapped_forward_target_genes.csv")
)

write_csv(
  reverse_receiver_ids,
  file.path(go_output_dir, "mapped_reverse_receiver_genes.csv")
)

write_csv(
  reverse_target_ids,
  file.path(go_output_dir, "mapped_reverse_target_genes.csv")
)

write_csv(
  forward_background_ids,
  file.path(go_output_dir, "mapped_forward_background_genes.csv")
)

write_csv(
  reverse_background_ids,
  file.path(go_output_dir, "mapped_reverse_background_genes.csv")
)

### 6. Make and save dotplots/barplots ---------------------------------------

# Convert GeneRatio from text such as "12/150" into a number.
calculate_gene_ratio <- function(gene_ratio) {
  map_dbl(
    gene_ratio,
    function(x) {
      parts <- str_split(x, "/", simplify = TRUE)
      as.numeric(parts[1]) / as.numeric(parts[2])
    }
  )
}

make_go_dotplot_and_barplot <- function(go_result, plot_title, output_prefix) {

  # Select one shared top 15 GO-term table.
  # Both the dotplot and barplot are made from this exact same table,
  # so the GO term list and p.adjust values are identical.
  plot_table <- as_tibble(go_result@result) %>%
    arrange(p.adjust) %>%
    slice_head(n = 15) %>%
    mutate(
      GeneRatio_number = calculate_gene_ratio(GeneRatio),
      Description = factor(Description, levels = rev(Description))
    )

  dot_plot <- ggplot(
    plot_table,
    aes(
      x = GeneRatio_number,
      y = Description,
      size = Count,
      color = p.adjust
    )
  ) +
    geom_point() +
    scale_color_gradient(low = "red", high = "blue") +
    labs(
      title = plot_title,
      x = "GeneRatio",
      y = NULL,
      color = "p.adjust",
      size = "Count"
    ) +
    theme_bw()

  bar_plot <- ggplot(
    plot_table,
    aes(
      x = Count,
      y = Description,
      fill = p.adjust
    )
  ) +
    geom_col() +
    scale_fill_gradient(low = "red", high = "blue") +
    labs(
      title = plot_title,
      x = "Count",
      y = NULL,
      fill = "p.adjust"
    ) +
    theme_bw()

  ggsave(
    file.path(go_output_dir, paste0(output_prefix, "_dotplot.png")),
    dot_plot,
    width = 10,
    height = 7,
    dpi = 300
  )

  ggsave(
    file.path(go_output_dir, paste0(output_prefix, "_barplot.png")),
    bar_plot,
    width = 10,
    height = 7,
    dpi = 300
  )
}

make_go_dotplot_and_barplot(
  go_forward_receiver,
  "SkMVECs to Myotubes - receiver upregulated genes",
  "SkMVECs_to_Myotubes_receiver_upregulated_GO_BP"
)

make_go_dotplot_and_barplot(
  go_forward_targets,
  "SkMVECs to Myotubes - NicheNet target genes",
  "SkMVECs_to_Myotubes_nichenet_targets_GO_BP"
)

make_go_dotplot_and_barplot(
  go_reverse_receiver,
  "Myotubes to SkMVECs - receiver upregulated genes",
  "Myotubes_to_SkMVECs_receiver_upregulated_GO_BP"
)

make_go_dotplot_and_barplot(
  go_reverse_targets,
  "Myotubes to SkMVECs - NicheNet target genes",
  "Myotubes_to_SkMVECs_nichenet_targets_GO_BP"
)




