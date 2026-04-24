# NicheNet's ligand activity analysis on a gene set of interest: predict active 
# ligand and their target genes.
# Following vignette: https://github.com/saeyslab/nichenetr/blob/master/vignettes/ligand_activity_geneset.md

'''
install.packages("devtools")
library(devtools)
install.packages("gitcreds")
install.packages("usethis")
usethis::create_github_token()
gitcreds::gitcreds_set()
system("git --version")
gitcreds::gitcreds_set()
devtools::install_github("saeyslab/nichenetr")
library(nichenetr)
vignette("ligand_activity_geneset", package="nichenetr")
'''
# Load packages
library(nichenetr)
library(tidyverse)

## Read in NicheNet networks
# Read in ligand-target prior model
lr_network <- readRDS("C:/Users/royde/OneDrive/Documenten/BIT11 Internship/lr_network_human_21122021.rds")
#url("https://zenodo.org/records/7074291/files/lr_network_human_21122021.rds")
ligand_target_matrix <- readRDS("C:/Users/royde/OneDrive/Documenten/BIT11 Internship/ligand_target_matrix_nsga2r_final.rds")
weighted_networks <- readRDS("C:/Users/royde/OneDrive/Documenten/BIT11 Internship/weighted_networks_nsga2r_final.rds")

# Keep only unique combinations of the columns ´from` and `to` 
lr_network <- lr_network %>% distinct(from,to)
#distinct(lr_network, from, to)
head(lr_network)

# Target genes in rows, ligands in columns
ligand_target_matrix[1:5, 1:5]

# Interactions and their weights in the ligand-receptor + signaling network
head(weighted_networks$lr_sig)

# Interactions and their weights in the gene regulatory network
head(weighted_networks$gr)

## Read in the expression data of interacting cells
hnscc_expression <- readRDS("C:/Users/royde/OneDrive/Documenten/BIT11 Internship/hnscc_expression.rds")
expression <- hnscc_expression$expression
sample_info <- hnscc_expression$sample_info # contains meta-information about the cells

# Convert to official gene symbols
colnames(expression) <- convert_alias_to_symbols(colnames(expression), "human", verbose = FALSE)

### 1. Define a set of potential ligands
# Sender-focused approach: CAF = sender and malignant cells = receivers
tumors_remove <- c("HN10","HN","HN12", "HN13", "HN24", "HN7", "HN8","HN23")

CAF_ids <- sample_info %>%
  filter(`Lymph node` == 0 & !(tumor %in% tumors_remove) &
           `non-cancer cell type` == "CAF") %>% pull(cell)
malignant_ids <- sample_info %>% filter(`Lymph node` == 0 &
                                          !(tumor %in% tumors_remove) &
                                          `classified  as cancer cell` == 1) %>% pull(cell)

expressed_genes_sender <- expression[CAF_ids,] %>%
  apply(2,function(x){10*(2**x - 1)}) %>%
  apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>% 
  names()

expressed_genes_receiver <- expression[malignant_ids,] %>%
  apply(2,function(x){10*(2**x - 1)}) %>%
  apply(2,function(x){log2(mean(x) + 1)}) %>% .[. >= 4] %>%
  names()

length(expressed_genes_sender)

length(expressed_genes_receiver)

# Filter the expressed ligands and receptors to only those that bind together
ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands <- intersect(ligands,expressed_genes_sender)

receptors <- lr_network %>% pull(to) %>% unique()
expressed_receptors <- intersect(receptors,expressed_genes_receiver)

potential_ligands <-  lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% 
  pull(from) %>% unique()

head(potential_ligands)

### 2. Define the gene set of interest and a background of genes
# The gene set of interest consists of genes for which the expression is possibly affected due to communication with other cells.
# The definition of this gene set depends on your research question! (crucial step)

# p-EMT gene set --> gene set of interest to investigate how CAFs can induce p-EMT in malignant cells.
geneset_oi <- readr::read_tsv(url("https://zenodo.org/record/3260758/files/pemt_signature.txt"),
                              col_names = "gene") %>%
  pull(gene) %>% .[. %in% rownames(ligand_target_matrix)] 

length(geneset_oi)

### 3. Define background genes
# malignant cells --> background set
background_expressed_genes <- expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]

length(background_expressed_genes)

### 4. Perform NicheNet ligand activity analysis
# assess how well each CAF-ligand can predict the p-EMT gene set compared to the background of expressed genes.
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                               background_expressed_genes = background_expressed_genes,
                                               ligand_target_matrix = ligand_target_matrix,
                                               potential_ligands = potential_ligands)
# ligands are ranked based on the area under the precision-recall curve (AUPR) between a ligand's target predictions and the observed transcriptional response.
# AUPR is the most informative measure to define ligand activity
(ligand_activities <- ligand_activities %>% arrange(-aupr_corrected) %>%
    mutate(rank = rank(desc(aupr_corrected))))

best_upstream_ligands <- ligand_activities %>% top_n(30, aupr_corrected) %>%
  arrange(-aupr_corrected) %>% pull(test_ligand)

best_upstream_ligands

### 5. Infer target genes and receptors of top-ranked ligands
# Active targets genes = highest regulatory potential
active_ligand_target_links_df <- best_upstream_ligands %>%
  lapply(get_weighted_ligand_target_links,
         geneset = geneset_oi,
         ligand_target_matrix = ligand_target_matrix,
         n = 200) %>% bind_rows()

active_ligand_target_links <- prepare_ligand_target_visualization(
  ligand_target_df = active_ligand_target_links_df,
  ligand_target_matrix = ligand_target_matrix,
  cutoff = 0.25)

order_ligands <- intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev()
order_targets <- active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links))

vis_ligand_target <- t(active_ligand_target_links[order_targets,order_ligands])

p_ligand_target_network <- make_heatmap_ggplot(vis_ligand_target, "Prioritized CAF-ligands", "p-EMT genes in malignant cells",
                                               color = "purple", legend_title = "Regulatory potential") +
  scale_fill_gradient2(low = "whitesmoke",  high = "purple")

p_ligand_target_network

# look at which receptors of the receiver cell population (malignant cells) can potentially bind to the prioritized ligands from the sender cell population (CAFs).
ligand_receptor_links_df <- get_weighted_ligand_receptor_links(
  best_upstream_ligands, expressed_receptors,
  lr_network, weighted_networks$lr_sig) 

vis_ligand_receptor_network <- prepare_ligand_receptor_visualization(
  ligand_receptor_links_df,
  best_upstream_ligands,
  order_hclust = "both") 

(make_heatmap_ggplot(t(vis_ligand_receptor_network), 
                     y_name = "Prioritized CAF-ligands", x_name = "Receptors expressed by malignant cells",  
                     color = "mediumvioletred", legend_title = "Prior interaction potential"))

### 6. Summary visualizations
library(RColorBrewer)
library(cowplot)
library(ggpubr)

# Ligand activity matrix
vis_ligand_aupr <- ligand_activities %>% filter(test_ligand %in% best_upstream_ligands) %>%
  column_to_rownames("test_ligand") %>% select(aupr_corrected) %>% arrange(aupr_corrected) %>% as.matrix(ncol = 1)

p_ligand_aupr <- make_heatmap_ggplot(vis_ligand_aupr,
                                     "Prioritized CAF-ligands", "Ligand activity",
                                     color = "darkorange", legend_title = "AUPR") + 
  theme(axis.text.x.top = element_blank())
p_ligand_aupr

# Combine plots


#############
p_ligand_receptor <- make_heatmap_ggplot(
  t(vis_ligand_receptor_network), 
  y_name = "Prioritized CAF-ligands", 
  x_name = "Receptors expressed by malignant cells",  
  color = "mediumvioletred", 
  legend_title = "Prior interaction potential"
)

figures_without_legend <- plot_grid(
  p_ligand_aupr + theme(legend.position = "none"),
  p_ligand_target_network + theme(
    legend.position = "none",
    axis.ticks = element_blank(),
    axis.title.y = element_blank()
  ),
  p_ligand_receptor + theme(legend.position = "none"),
  ncol = 1,
  align = "v",
  rel_heights = c(1, 2, 2)
)

figures_without_legend

legend <- get_legend(p_ligand_aupr)

final_figure <- plot_grid(
  figures_without_legend,
  legend,
  ncol = 2,
  rel_widths = c(0.9, 0.1)
)

final_figure



