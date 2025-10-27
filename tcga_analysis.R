
# tcga_analysis.R
# Transcriptomic analysis of TCGA-THCA for paper:
# "Two Cohorts, One Network: Consensus Master Regulators 
# Orchestrating Papillary Thyroid Carcinoma"
# Hugo Tovar, National Institute of Genomic Medicine, Mexico 
# hatovar@inmegen.gob.mx

#############################################
# 1. Load required libraries
#############################################
suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(DESeq2)
  library(dplyr)
  library(org.Hs.eg.db)
  library(tibble)
  library(viper)
  library(clusterProfiler)
  library(forcats)
  library(ggplot2)
  library(readr)
  library(limma)
  library(edgeR)
})
# Load helper functions
source("helpers.R")

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
raw_counts <- GDCprepare(query, summarizedExperiment = FALSE) %>%
  as.data.frame()

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
# 4. Collapse duplicated gene names using the statistical mode
#############################################
collapsed <- filtered_counts %>% 
  group_by(gene_name) %>%
  summarise(across(everything(), 
                    ~getmode(.)), 
            .groups = "drop")

#############################################
# 5. Filter low-expressed genes and normalize (TMM + CPM)
#############################################
keep_genes <- rowSums(collapsed[,-1] > 10) > 0.8 * ncol(collapsed[,-1])
exp_matrix <- as.data.frame(collapsed[keep_genes, ])
rownames(exp_matrix) <- exp_matrix$gene_name
exp_matrix <- exp_matrix[, -1]

sample_ids <- colnames(exp_matrix)
tumors <- sample_ids[grepl("-01A", sample_ids)]
healthy <- sample_ids[grepl("-11A", sample_ids)]

dge <- DGEList(counts = exp_matrix)
dge <- calcNormFactors(dge, method = "TMM")
cpm_matrix <- cpm(dge, normalized.lib.sizes = TRUE)
save(tumors, healthy, cpm_matrix, file = "tcga_sample_groups_n_matrix.RData")
saveRDS(cpm_matrix[, tumors], paste0(rdsFolder, "inmat_TCGA.rds"))

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



#############################################
# 7. Master Regulator Analysis (implemented in VIPER)
#############################################


set.seed(123)
regulon <- aracne2regulon("../tiroides_networks/tcga_tumor_network_p1E-8.txt", cpm_matrix[, tumors])

signature <- rowTtest(cpm_matrix[, tumors], cpm_matrix[, healthy])
signature <- (qnorm(signature$p.value / 2, lower.tail = FALSE) *
              sign(signature$statistic))[, 1]

nullmodel <- ttestNull(cpm_matrix[, tumors], 
                       cpm_matrix[, healthy], 
                       per = 1000, repos = TRUE, cores = 30)

mrs <- msviper(signature, regulon, nullmodel, cores = 30)
mrs_all <- viper_mrsTopTable(mrs, p_threshold = 1)

# Save results
saveRDS(mrs, paste0(rdsFolder, "tcga_mrs.rds"))
write_tsv(as_tibble(mrs_all), paste0(outputsFolder, "all_mrs_table_TCGA.tsv"))
saveRDS(mrs_all, paste0(rdsFolder, "tcga_mrs_table.rds"))

# save.image("tcga_session.RData")
