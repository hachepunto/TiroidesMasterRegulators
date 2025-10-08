
# geo_analysis.R
# Transcriptomic analysis of GEO GSE33630 for the subtype "Papillary adenocarcinoma, NOS"
# Hugo Tovar, National Inmstitute of Genomic Medicine, Mexico hatovar@inmegen.gob.mx

#############################################
# 1. Load required libraries and setup folders (GEO analysis)
#############################################

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(affy)
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
  library(RColorBrewer)
  library(gplots)
  library(GEOquery)
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

# Download GSE33630 manually if not already present ~850 MB
geodataFolder <- paste0(getwd(), "/GEOdata/GSE33630/")
dir.create(geodataFolder, showWarnings = FALSE)
tar_file <- paste0(geodataFolder, "GSE33630_RAW.tar")
if (!file.exists("GEOdata/GSE33630/GSE33630_RAW.tar")) {
  download.file(
    url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE33nnn/GSE33630/suppl/GSE33630_RAW.tar",
    destfile = tar_file,
    mode = "wb"
  )
}

cel_directory <- paste0(geodataFolder, "CEL_files")

# Decompression of .tar archive containing CEL files
untar(tar_file, exdir = cel_directory)

# If you encounter errors during download or get a corrupted TAR file,
# you can manually download it from the following GEO page:
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE33630

#############################################
# 2. Normalize raw CEL files using RMA
#############################################

# Read all .CEL.gz files from the directory
cel_files <- list.files(cel_directory, pattern = ".CEL.gz$", full.names = TRUE)

# Read raw microarray data
raw_affy <- ReadAffy(filenames = cel_files)

# Normalize data using Robust Multi-array Average (RMA)
eset <- rma(raw_affy)

#############################################
# 3. Download metadata and annotate samples
#############################################

# Download GEO metadata for GSE33630
gse <- getGEO("GSE33630", GSEMatrix = TRUE, getGPL = FALSE)
metadata <- pData(phenoData(gse[[1]]))

# Ensure sample names in metadata match those in eset
sampleNames(eset) <- substr(sampleNames(eset), 1, 9)

# Match and filter metadata to match expression matrix
metadata <- metadata[match(sampleNames(eset), metadata$geo_accession), ]

# Assign group label based on 'pathological' column
metadata$Group <- ifelse(metadata$pathological == "papillary thyroid carcinoma (PTC)", "PTC",
                         ifelse(metadata$pathological == "patient-matched non-tumor control", "Normal", NA))

# Keep only samples with defined group
eset <- eset[, !is.na(metadata$Group)]
metadata <- metadata[!is.na(metadata$Group), ]

# Attach metadata to expression set
pData(eset) <- metadata

#############################################
# 4. Differential expression analysis with limma
#############################################

# Define the design matrix for comparing groups
design <- model.matrix(~ 0 + metadata$Group)
colnames(design) <- c("Normal", "PTC")

# Create contrast matrix to compare PTC vs Normal
contrast.matrix <- makeContrasts(PTC - Normal, levels = design)

# Fit linear model
fit <- lmFit(eset, design)

# Apply contrast and empirical Bayes moderation
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Extract all differential expression results
results <- topTable(fit2, adjust = "fdr", number = Inf)

# Convert Affymetrix probe IDs to gene symbols using gprofiler2
affy_ids <- rownames(results)
gconvert <- gprofiler2::gconvert(query = affy_ids,
                                 target = "HGNC",
                                 mthreshold = 1,
                                 filter_na = FALSE)

# Sort and align converted gene symbols with results
class(gconvert[, 1]) <- "integer"
gconvert <- gconvert[sort.list(gconvert[, 1]), ]
gene_symbol <- gconvert[, 5]

# ---- Filter out probes with missing gene symbols ----
valid_affy <- !is.na(gene_symbol)
results <- results[valid_affy, ]
gene_symbol <- gene_symbol[valid_affy]
results$gene_name <- gene_symbol

# Collapse probes by keeping the one with the highest B-statistic per gene symbol
results_with_row <- results %>% 
  mutate(rownum = row_number())

Bmax <- results_with_row %>%
  group_by(gene_name) %>%
  filter(B == max(B)) %>%
  ungroup() %>%
  pull(rownum)

results_collapsed <- results[Bmax, ]

# Save collapsed results
write_tsv(results_collapsed,
          file = paste0(outputsFolder,"DEGs_limma_GEO.tsv"))


# ----------------------------
# Volcano plot
# ----------------------------
pdf(paste0(plotsFolder, "volcano_GEO.pdf"), width = 8, height = 8)
EnhancedVolcano(results_collapsed,
                lab = results_collapsed$gene_name,
                x = 'logFC',
                y = 'adj.P.Val',
                title = 'GSE33630: PTC vs Normal',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 2.5,
                labSize = 3)
dev.off()


# ----------------------------
# Heatmap of top 100 DEGs
# ----------------------------

# Extract top 100 DEGs by adjusted p-value
top_degs <- results_collapsed %>%
  arrange(adj.P.Val) %>%
  slice(1:100) %>%
  pull(gene_name)

# Create matrix of expression data (log2-transformed)
expr_matrix <- exprs(eset)
expr_matrix <- expr_matrix[match(rownames(results_collapsed), rownames(expr_matrix)), ]
rownames(expr_matrix) <- results_collapsed$gene_name
log_expr <- log2(expr_matrix[top_degs, ] + 1)

# Assign colors to sample groups
group_colors <- ifelse(metadata$Group == "PTC", "#F8766D", "#00BFC4")

# Generate heatmap with z-score per row
pdf(paste0(plotsFolder, "heatmap_GEO.pdf"), width = 8, height = 8)
heatmap.2(
  as.matrix(log_expr),
  scale = "row",                  # z-score per gene
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
  main = "Top 100 DEGs - GSE33630"
)
legend("topright", legend = c("PTC", "Normal"),
       fill = c("#F8766D", "#00BFC4"), cex = 0.8, border = NA, bty = "n")
dev.off()

#############################################
# 5. Prepare expression matrix for ARACNe
#############################################

# Extract expression matrix from normalized ExpressionSet
geo_matrix <- exprs(eset)

# Annotate Affymetrix probe IDs with gene symbols using g:Profiler
gconvert <- gprofiler2::gconvert(query = rownames(geo_matrix),
                                 target = "HGNC",
                                 mthreshold = 1,
                                 filter_na = FALSE)

# Ensure correct order for alignment
class(gconvert[, 1]) <- "integer"
gconvert <- gconvert[sort.list(gconvert[, 1]), ]
# Extract HGNC symbols, fall back to probe IDs if unavailable
genesymbol <- gconvert[, 5]

# ---- Filter out probes with missing gene symbols ----
valid_idx <- !is.na(genesymbol)
geo_matrix <- geo_matrix[valid_idx, ]
genesymbol <- genesymbol[valid_idx]


# Create design matrix for differential expression: PTC (case) vs Normal (control)
design <- matrix(rep(0, 2 * nrow(metadata)), nrow = nrow(metadata))
colnames(design) <- c('case', 'control')
rownames(design) <- rownames(metadata)
design[metadata$Group == "PTC", "case"] <- 1
design[metadata$Group == "Normal", "control"] <- 1


# Build contrast matrix and fit linear model
cont.matrix <- makeContrasts('case - control', levels = design)
fit <- lmFit(geo_matrix, design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# Retrieve full statistics (unsorted)
statistics <- topTable(fit2, coef = 1, adjust = "fdr", number = nrow(geo_matrix), sort.by = "none")


# Collapse probes by selecting the one with the highest B-statistic per gene
Bmax <- statistics %>%
  mutate(gene_name = genesymbol) %>%
  mutate(rownum = row_number()) %>%
  group_by(gene_name) %>%
  filter(B == max(B)) %>%
  ungroup() %>%
  pull(rownum)

# Subset expression matrix to selected probes and assign gene symbols as rownames
collapsed <- geo_matrix[Bmax, ]
rownames(collapsed) <- genesymbol[Bmax]

# ---- Filter out lowly expressed genes (mean log2 expression <= 4) ----
mean_expr <- rowMeans(collapsed)
collapsed <- collapsed[mean_expr > 4, ]

# ---- Filter out low-variance genes (bottom 25% by variance) ----
gene_var <- apply(collapsed, 1, var)
var_threshold <- quantile(gene_var, probs = 0.25)
collapsed <- collapsed[gene_var > var_threshold, ]


tumors <- which(metadata$Group == "PTC")
healty <- which(metadata$Group == "Normal")

# Extract only PTC samples
ptc_matrix <- collapsed[, tumors]

# Convert to data.frame with gene names as a proper column
ptc_df <- ptc_matrix %>%
  as.data.frame() %>%
  rownames_to_column(var = "gene_name")

# Save as tab-delimited file
write.table(ptc_df,
            file = "inmat_GEO.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)

#############################################
# 6. GSEA: GO (BP and MF) and KEGG using clusterProfiler
#############################################

# Run GSEA for GO terms (BP and MF)
gsea_go_out <- run_gsea_go(df = results_collapsed)

# Save GO results
write_tsv(gsea_go_out$gsea_result@result,
          paste0(outputsFolder, "geo_gsea_go_bp_mf.tsv"))

# Ridgeplot for GO
pdf(paste0(plotsFolder, "geo_gsea_go.pdf"), width = 12, height = 15)
plot_ridge_panels(gsea_go_out$gsea_result, 20, title = "GEO GSEA GO Ridgeplot")
dev.off()

# Save GO RDS
saveRDS(gsea_go_out, paste0(rdsFolder, "geo_gsea_go.rds"))

# Run GSEA for KEGG
kegg_gsea <- gseKEGG(geneList = gsea_go_out$gene_list,
                     organism = "hsa",
                     pvalueCutoff = 0.05,
                     verbose = FALSE)

# Save KEGG results
write_tsv(kegg_gsea@result, paste0(outputsFolder, "geo_gsea_kegg.tsv"))

# Ridgeplot for KEGG
pdf(paste0(plotsFolder, "geo_gsea_kegg.pdf"), width = 12, height = 15)
ridgeplot(kegg_gsea, showCategory = 20, fill = "p.adjust") +
  ggtitle("GEO GSEA KEGG Ridgeplot")
dev.off()

# Save KEGG RDS
saveRDS(kegg_gsea, paste0(rdsFolder, "geo_gsea_kegg.rds"))


# Run GSEA: Hallmarks (MSigDB) using clusterProfiler
# Load Hallmark data (GSEA Broad GMT)
hallmark_gmt <- read.gmt("extra_data/h.all.v2025.1.Hs.symbols.gmt")
hallmark_names <- read_tsv("extra_data/hallmarks_names.txt")

rank_sym_geo <- gsea_go_out$gene_mapping %>%
  inner_join(
    AnnotationDbi::select(org.Hs.eg.db,
                          keys   = unique(.$ENTREZID),
                          keytype= "ENTREZID",
                          columns= "SYMBOL") %>% as_tibble(),
    by = "ENTREZID"
  ) %>%
  filter(!is.na(SYMBOL)) %>%
  mutate(SYMBOL = toupper(SYMBOL)) %>%
  group_by(SYMBOL) %>%
  summarize(logFC = mean(logFC), .groups = "drop") %>%
  arrange(desc(logFC))

gene_list_sym_geo <- rank_sym_geo$logFC
names(gene_list_sym_geo) <- rank_sym_geo$SYMBOL
gene_list_sym_geo <- sort(gene_list_sym_geo, decreasing = TRUE)

# 2) Asegurar TERM2GENE/TERM2NAME en formato limpio
hallmark_t2g <- hallmark_gmt %>%
  as_tibble() %>%
  transmute(term = as.character(term), gene = toupper(gene)) %>%
  distinct(term, gene)

hallmark_t2n <- hallmark_names %>%
  transmute(term = ID, name = label_ready)

# 3) Correr GSEA con símbolos
gsea_hallmark_geo <- clusterProfiler::GSEA(
  geneList     = gene_list_sym_geo,   # <-- ahora SYMBOL
  TERM2GENE    = hallmark_t2g,
  TERM2NAME    = hallmark_t2n,
  pvalueCutoff = 1,
  verbose      = FALSE
)

readr::write_tsv(as.data.frame(gsea_hallmark_geo@result),
                   file = paste0(outputsFolder, "geo_gsea_hallmarks.tsv"))

pdf(paste0(plotsFolder, "geo_gsea_hallmarks.pdf"), width = 10, height = 8)
print(enrichplot::ridgeplot(gsea_hallmark_geo, showCategory = 20, fill = "p.adjust") +
ggplot2::ggtitle("GEO GSEA Hallmarks (Top 20)"))
dev.off()

ggsave(filename = paste0(plotsFolder, "geo_gsea_hallmarks.png"),
         plot     = enrichplot::ridgeplot(gsea_hallmark_geo, showCategory = 20, fill = "p.adjust") +
                    ggplot2::ggtitle("GEO GSEA Hallmarks (Top 20)"),
         width = 10, height = 8, dpi = 300)

saveRDS(gsea_hallmark_geo, paste0(rdsFolder, "geo_gsea_hallmarks.rds"))

#############################################
# 7. Master Regulator Analysis (VIPER)
#############################################

regulon <- aracne2regulon("geo_tumor_network.txt", collapsed[, tumors])

signature <- rowTtest(collapsed[, tumors], collapsed[, healty])
signature <- (qnorm(signature$p.value / 2, lower.tail = FALSE) *
              sign(signature$statistic))[, 1]

nullmodel <- ttestNull(collapsed[, tumors], 
                       collapsed[, healty], 
                       per = 1000, 
                       repos = TRUE)

mrs <- msviper(signature, regulon, nullmodel)
mrs_all <- viper_mrsTopTable(mrs, p_threshold = 1)

# Save top results
saveRDS(mrs, paste0(rdsFolder, "geo_mrs.rds"))
write_tsv(as_tibble(mrs_all), paste0(outputsFolder, "all_mrs_table_GEO.tsv"))
saveRDS(mrs_all, paste0(rdsFolder, "geo_mrs_table.rds"))

# Plots
pdf(paste0(plotsFolder, "mrs_plot_GEO.pdf"), width = 7, height = 8)
plot(mrs, 20, density = 100, color = c("#0571b0", "#ca0020"),
     smooth = 0, cex = 0, bins = 800, sep = 1, hybrid = TRUE, gama = 3)
dev.off()

png(paste0(plotsFolder, "mrs_plot_GEO.png"), width = 960, height = 1200, res = 150)
plot(mrs, 20, density = 100, color = c("#0571b0", "#ca0020"),
     smooth = 0, cex = 0, bins = 800, sep = 1, hybrid = TRUE, gama = 3)
dev.off()


# save.image("geo_session.RData")