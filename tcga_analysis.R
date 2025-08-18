
# tcga_analysis.R
# Transcriptomic analysis of TCGA-THCA for the subtype "Papillary adenocarcinoma, NOS"
# Hugo Tovar, National Inmstitute of Genomic Medicine, Mexico hatovar@inmegen.gob.mx

#############################################
# 1. Load required libraries
#############################################
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(DESeq2)
  library(dplyr)
  library(org.Hs.eg.db)
  library(tibble)
  library(corto)
  library(viper)
  library(clusterProfiler)
  library(forcats)
  library(ggplot2)
  library(readr)
  library(limma)
  library(edgeR)
  library(EnhancedVolcano)
  library(pheatmap)
  library(RColorBrewer)
  library(gplots)
})
# Load helper functions
source("scripts/helpers.R")

# Create output directories
outputsFolder <- paste0(getwd(), "/results/")
plotsFolder <- paste0(getwd(), "/plots/")
rdsFolder <- paste0(getwd(), "/rds/")
dir.create(rdsFolder, showWarnings = FALSE)
dir.create(outputsFolder, showWarnings = FALSE)
dir.create(plotsFolder, showWarnings = FALSE)

#############################################
# 2. Download and prepare raw count data
#############################################
query <- GDCquery(
  project = "TCGA-THCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

download_if_needed(query)
raw_counts <- GDCprepare(query, summarizedExperiment = FALSE) %>% as.data.frame()

# Keep only gene_name, gene_id, gene_type and unstranded counts
count_cols <- grep("^unstranded_", colnames(raw_counts), value = TRUE)
raw_counts_unstranded <- raw_counts[, c("gene_id", "gene_name", "gene_type", count_cols)]

# Remove non-gene rows like N_noFeature, N_unmapped, etc.
raw_counts_unstranded <- raw_counts_unstranded[!grepl("^N_", raw_counts$gene_id), ]

# Keep only protein_coding genes
raw_counts_clean <- raw_counts_unstranded[raw_counts_unstranded$gene_type == "protein_coding", ]

# Drop gene_type column (no longer needed)
raw_counts_clean <- raw_counts_clean[, -which(colnames(raw_counts_clean) == "gene_type")]

#############################################
# 3. Filter samples by subtype "Papillary adenocarcinoma, NOS"
#############################################
clinical <- GDCquery_clinic(project = "TCGA-THCA", type = "clinical")
nos_ids <- clinical %>%
  filter(primary_diagnosis == "Papillary adenocarcinoma, NOS") %>%
  pull(submitter_id)

clean_colnames <- gsub("^.*TCGA", "TCGA", colnames(raw_counts_clean))
col_short_ids <- sapply(strsplit(clean_colnames, "-"), function(x) paste(x[1:3], collapse = "-"))
keep_cols <- col_short_ids %in% nos_ids
keep_cols[2] <- TRUE  # keep gene_name column

filtered_counts <- raw_counts_clean[, keep_cols]

#############################################
# 4. Collapse duplicated gene names using statistical mode
#############################################
collapsed <- filtered_counts %>% 
  group_by(gene_name) %>%
  summarise(across(everything(), ~getmode(.)), .groups = "drop")

#############################################
# 5. Filter low-expressed genes and normalize (TMM + CPM)
#############################################
keep_genes <- rowSums(collapsed[,-1] > 10) > 0.8 * ncol(collapsed[,-1])
exp_matrix <- as.data.frame(collapsed[keep_genes, ])
rownames(exp_matrix) <- exp_matrix$gene_name
exp_matrix <- exp_matrix[, -1]

sample_ids <- colnames(exp_matrix)
tumors <- sample_ids[grepl("-01A", sample_ids)]
healty <- sample_ids[grepl("-11A", sample_ids)]

dge <- DGEList(counts = exp_matrix)
dge <- calcNormFactors(dge, method = "TMM")
cpm_matrix <- cpm(dge, normalized.lib.sizes = TRUE)

# Save CPM matrix (tumor samples only) for ARACNe-AP
inmat_NOS <- as.data.frame(cpm_matrix[, tumors]) %>%
  rownames_to_column(var = "gene_name")
write_tsv(inmat_NOS, "inmat_TCGA.txt")

#############################################
# 6. Differential Expression Analysis (limma + voom)
#############################################
group <- factor(ifelse(colnames(cpm_matrix) %in% tumors, "tumor", "normal"))
design <- model.matrix(~group)

# voom transformation and linear model fitting
v <- voom(dge, design, plot = FALSE)
fit <- lmFit(v, design)
fit <- eBayes(fit)

# extract all DEGs table
deg_results <- topTable(fit, coef = "grouptumor", number = Inf, sort.by = "P")
deg_results <- rownames_to_column(deg_results, var = "gene_name")
write_tsv(deg_results, paste0(outputsFolder, "DEGs_limma_TCGA.tsv"))

# ----------------------------
# Volcano plot
# ----------------------------
pdf(paste0(plotsFolder, "volcano_TCGA.pdf"), width = 8, height = 8)
EnhancedVolcano(deg_results,
                lab = deg_results$gene_name,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'TCGA Papillary NOS',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.5,
                labSize = 3)
dev.off()

# ----------------------------
# Heatmap of top 100 DEGs
# ----------------------------

# Extract top 100 DEGs by adjusted p-value
top_degs <- deg_results %>%
  arrange(adj.P.Val) %>%
  slice(1:100) %>%
  pull(gene_name)

# Extract raw matrix and compute log2(CPM + 1)
log_cpm <- log2(cpm_matrix[top_degs, ] + 1)

# Generate color vector by condition
group_colors <- ifelse(colnames(log_cpm) %in% tumors, "#F8766D", "#00BFC4")

# Add color bar above heatmap
col_side_colors <- matrix(group_colors, nrow = 1)

# Generate heatmap with Z-score per gene
pdf(paste0(plotsFolder, "heatmap_TCGA.pdf"), width = 8, height = 8)
heatmap.2(
  as.matrix(log_cpm),
  scale = "row",                  # z-transformation per row
  trace = "none",
  col = colorRampPalette(rev(brewer.pal(11, "RdYlBu")))(100),
  ColSideColors = group_colors,
  dendrogram = "both",
  key = TRUE,
  key.title = "Z-score",
  key.xlab = "Expression",
  cexRow = 0.6,
  cexCol = 0.5,
  margins = c(5, 9),
  main = "Top 100 DEGs TCGA Papillary NOS"
)
legend("topright", legend = c("Tumor", "Normal"),
       fill = c("#F8766D", "#00BFC4"), cex = 0.8, border = NA, bty = "n")
dev.off()

#############################################
# 7. GSEA: GO (BP and MF) and KEGG using clusterProfiler
#############################################

# Run GSEA for GO
gsea_go_out <- run_gsea_go(df = deg_results)

# Save GO results
write_tsv(gsea_go_out$gsea_result@result,
          paste0(outputsFolder, "tcga_gsea_go_bp_mf.tsv"))

saveRDS(gsea_go_out, paste0(rdsFolder, "tcga_gsea_go_out.rds"))

# Ridgeplot for GO
pdf(paste0(plotsFolder, "tcga_gsea_go.pdf"), width = 12, height = 15)
plot_ridge_panels(gsea_go_out$gsea_result, 20, title = "TCGA GSEA GO Ridgeplot")
dev.off()

# GSEA for KEGG
kegg_gsea <- gseKEGG(geneList = gsea_go_out$gene_list,
                     organism = "hsa",
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

# Save KEGG results
write_tsv(kegg_gsea@result, paste0(outputsFolder, "tcga_gsea_kegg.tsv"))

saveRDS(kegg_gsea, paste0(rdsFolder, "tcga_gsea_kegg_out.rds"))

# Ridgeplot for KEGG
pdf(paste0(plotsFolder, "tcga_gsea_kegg.pdf"), width = 12, height = 15)
ridgeplot(kegg_gsea, showCategory = 20, fill = "p.adjust") +
  ggplot2::ggtitle("TCGA GSEA KEGG Ridgeplot")
dev.off()

#############################################
# 8. Master Regulator Analysis (VIPER)
#############################################

regulon <- aracne2regulon("tcga_tumor_network.txt", cpm_matrix[, tumors])

signature <- rowTtest(cpm_matrix[, tumors], cpm_matrix[, healty])
signature <- (qnorm(signature$p.value / 2, lower.tail = FALSE) *
              sign(signature$statistic))[, 1]

nullmodel <- ttestNull(cpm_matrix[, tumors], 
                       cpm_matrix[, healty], 
                       per = 1000, repos = TRUE)

mrs <- msviper(signature, regulon, nullmodel)
mrs_all <- viper_mrsTopTable(mrs, p_threshold = 1)

# Save results
saveRDS(mrs, paste0(rdsFolder, "tcga_mrs.rds"))
write_tsv(as_tibble(mrs_all), paste0(outputsFolder, "all_mrs_table_TCGA.tsv"))
saveRDS(mrs_all, paste0(rdsFolder, "tcga_mrs_table.rds"))

# Plots
pdf(paste0(plotsFolder, "mrs_plot_TCGA.pdf"), width = 7, height = 8)
plot(mrs, density = 100, color = c("#0571b0", "#ca0020"),
     smooth = 0, cex = 0, bins = 800, sep = 1, hybrid = TRUE, gama = 3)
dev.off()

png(paste0(plotsFolder, "mrs_plot_TCGA.png"), width = 960, height = 1200, res = 150)
plot(mrs, density = 100, color = c("#0571b0", "#ca0020"),
     smooth = 0, cex = 0, bins = 800, sep = 1, hybrid = TRUE, gama = 3)
dev.off()

#############################################
# 9. ORA: Hallmark enrichment with MSigDB gene sets
#############################################

# Universe of expressed genes used in ORA
gene_universe <- rownames(cpm_matrix)

# Load Hallmark data (GSEA Broad GMT)
hallmark_gmt <- read.gmt("extra_data/h.all.v2025.1.Hs.symbols.gmt")

# Extract regulons for the top 10 MRTFs
top_mrs <- mrs_all  %>%  slice(1:20)
regulon_df <- getregulon(regulon, tf_list = top_mrs$TF)

# Run ORA mode 1: per TF
ora_tf <- ora_hallmarks(regulon_df = regulon_df,
                        hallmark_gmt = hallmark_gmt,
                        universe = gene_universe,
                        mode = "by_tf")
write_tsv(ora_tf, paste0(outputsFolder, "TCGA_ORA_hallmark_byTF.tsv"))
saveRDS(ora_tf, paste0(rdsFolder, "tcga_ora_hallmark_byTF.rds"))

# Run ORA mode 2: all targets combined
ora_all <- ora_hallmarks(regulon_df = regulon_df,
                         hallmark_gmt = hallmark_gmt,
                         universe = gene_universe,
                         mode = "all_targets")
write_tsv(ora_all, paste0(outputsFolder, "TCGA_ORA_hallmark_allTargets.tsv"))
saveRDS(ora_all, paste0(rdsFolder, "tcga_ora_hallmark_alltargets.rds"))

#############################################
# 10. ORA: GO + KEGG enrichment per TF and all targets
#############################################

# ORA per TF with GO (ALL categories)
ora_go_by_tf <- ora_go(regulon_df,
                       universe = gene_universe,
                       mode = "by_tf",
                       organism = org.Hs.eg.db)

# ORA with GO for all targets combined
ora_go_all_targets <- ora_go(regulon_df,
                             universe = gene_universe,
                             mode = "all_targets",
                             organism = org.Hs.eg.db)

# ORA per TF with KEGG
ora_kegg_by_tf <- ora_kegg(regulon_df,
                           universe = gene_universe,
                           mode = "by_tf")

# ORA with KEGG for all targets combined
ora_kegg_all_targets <- ora_kegg(regulon_df,
                                 universe = gene_universe,
                                 mode = "all_targets")

# Save combined results to files
write_tsv(ora_go_by_tf, paste0(outputsFolder, "tcga_ora_go_by_tf.tsv"))
write_tsv(ora_go_all_targets, paste0(outputsFolder, "tcga_ora_go_all_targets.tsv"))
write_tsv(ora_kegg_by_tf, paste0(outputsFolder, "tcga_ora_kegg_by_tf.tsv"))
write_tsv(ora_kegg_all_targets, paste0(outputsFolder, "tcga_ora_kegg_all_targets.tsv"))

saveRDS(ora_go_by_tf, paste0(rdsFolder, "tcga_ora_go_by_tf.rds"))
saveRDS(ora_go_all_targets, paste0(rdsFolder, "tcga_ora_go_all_targets.rds"))
saveRDS(ora_kegg_by_tf, paste0(rdsFolder, "tcga_ora_kegg_by_tf.rds"))
saveRDS(ora_kegg_all_targets, paste0(rdsFolder, "tcga_ora_kegg_all_targets.rds"))

# save.image("TCGA_session.RData")