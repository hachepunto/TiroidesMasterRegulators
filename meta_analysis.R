# meta_analysis.R
# Meta-analysis of master regulator results from TCGA and GEO
# Hugo Tovar, National Institute of Genomic Medicine, Mexico

#############################################
# 1. Load required libraries and set up folders
#############################################

suppressPackageStartupMessages({
 library(tidyverse)
 library(tibble)
 library(magrittr)
 library(metap)
 library(EnhancedVolcano)
 library(viper)
 library(igraph)
 library(clusterProfiler)
 library(ggtext)
 library(ComplexUpset)
 library(ggplot2)
 library(scales)
 library(cowplot)
 library(patchwork)
 library(dplyr)
 library(tidyr)
 library(grid)
})

# Load helper functions
source("scripts/helpers.R")

# Create output directories
outputsFolder <- paste0(getwd(), "/meta_results/")
plotsFolder <- paste0(getwd(), "/meta_results/plots/")
dir.create(outputsFolder, showWarnings = FALSE)
dir.create(plotsFolder, showWarnings = FALSE)

#############################################
# 2. Load GSEA results and perform meta-analysis
#############################################

tcga_gsea_go <- readRDS("rds/tcga_gsea_go.rds")
tcga_gsea_kegg <- readRDS("rds/tcga_gsea_kegg.rds")
tcga_gsea_hallmarks <- readRDS("rds/tcga_gsea_hallmarks.rds")
geo_gsea_go <- readRDS("rds/geo_gsea_go.rds")
geo_gsea_kegg <- readRDS("rds/geo_gsea_kegg.rds")
geo_gsea_hallmarks  <- readRDS("rds/geo_gsea_hallmarks.rds")

# ------------ GO ------------
tcga_go_df <- get_gsea_df(tcga_gsea_go)  %>% select_gsea_cols(is_go = TRUE,  suffix = "tcga")
geo_go_df  <- get_gsea_df(geo_gsea_go)   %>% select_gsea_cols(is_go = TRUE,  suffix = "geo")

go_meta <- combine_two_gsea(
  df_geo  = geo_go_df,
  df_tcga = tcga_go_df,
  keys    = c("ONTOLOGY","ID","Description")
)

filt_go_meta <- go_meta %>%
  filter(!is.na(meta_pval),
         meta_padj < 0.05,
         concordant_sign) %>%
  arrange(desc(abs(mean_NES)), Description) %>%
  mutate(
    full_label = paste(ONTOLOGY, Description, sep = ": "),
    full_label = factor(full_label, levels = unique(full_label))
  )

write_tsv(filt_go_meta, file.path(outputsFolder, "meta_gsea_GO.tsv"))

go_top <- filt_go_meta %>%
  filter(is.finite(meta_padj), is.finite(mean_NES)) %>%
  arrange(desc(abs(mean_NES))) %>%
  slice_head(n = 30) %>%
  mutate(term = full_label) %>%
  mutate(term = factor(term, levels = rev(unique(term))))

p_go <- ggplot(go_top, aes(x = mean_NES, y = term)) +
  geom_segment(aes(x = 0, xend = mean_NES, y = term, yend = term),
               linewidth = 1.0, color = "grey60") +
  geom_point(aes(color = meta_padj), size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  scale_color_gradient(low = "#2c7fb8", high = "#253494", trans = "reverse",
                       name = "Meta-adjusted p-value") +
  labs(title = "Meta-GSEA GO (TCGA + GEO)", x = "Mean NES", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.y        = element_text(size = 8),
        axis.text.x        = element_text(size = 8),
        legend.position    = "right",
        plot.title         = element_text(face = "bold", size = 15))

ggsave(file.path(plotsFolder, "meta_gsea_GO_lollipop.pdf"), p_go, width = 16, height = 4, units = "in")
ggsave(file.path(plotsFolder, "meta_gsea_GO_lollipop.png"), p_go, width = 16, height = 4, units = "in", dpi = 300)


# ------------ KEGG ------------
tcga_kegg_df <- get_gsea_df(tcga_gsea_kegg) %>% select_gsea_cols(is_go = FALSE, suffix = "tcga")
geo_kegg_df  <- get_gsea_df(geo_gsea_kegg)  %>% select_gsea_cols(is_go = FALSE, suffix = "geo")

kegg_meta <- combine_two_gsea(
  df_geo  = geo_kegg_df,
  df_tcga = tcga_kegg_df,
  keys    = c("ID","Description")
)

filt_kegg_meta <- kegg_meta %>%
  filter(!is.na(meta_pval),
         meta_padj < 0.05,
         concordant_sign) %>%
  arrange(desc(abs(mean_NES)), Description) %>%
  mutate(
    full_label = factor(Description, levels = unique(Description))
  )

write_tsv(filt_kegg_meta, file.path(outputsFolder, "meta_gsea_KEGG.tsv"))

kegg_top <- filt_kegg_meta %>%
  filter(is.finite(meta_padj), is.finite(mean_NES)) %>%
  arrange(desc(abs(mean_NES))) %>%
  slice_head(n = 30) %>%
  mutate(term = full_label) %>% 
  mutate(term = factor(term, levels = rev(unique(term))))

p_kegg <- ggplot(kegg_top, aes(x = mean_NES, y = term)) +
  geom_segment(aes(x = 0, xend = mean_NES, y = term, yend = term),
               linewidth = 1.0, color = "grey60") +
  geom_point(aes(color = meta_padj), size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  scale_color_gradient(low = "#2c7fb8", high = "#253494", trans = "reverse",
                       name = "Meta-adjusted p-value") +
  labs(title = "Meta-GSEA KEGG (TCGA + GEO)", x = "Mean NES", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.y        = element_text(size = 8),
        axis.text.x        = element_text(size = 8),
        legend.position    = "right",
        plot.title         = element_text(face = "bold", size = 15))

ggsave(file.path(plotsFolder, "meta_gsea_KEGG_lollipop.pdf"), p_kegg, width = 16, height = 4, units = "in")
ggsave(file.path(plotsFolder, "meta_gsea_KEGG_lollipop.png"), p_kegg, width = 16, height = 4, units = "in", dpi = 300)

# ------------ HALLMARKS ------------
# Load Hallmark data (GSEA Broad GMT)
hallmark_gmt <- read.gmt("extra_data/h.all.v2025.1.Hs.symbols.gmt")
hallmark_names <- read_tsv("extra_data/hallmarks_names.txt")

tcga_hm_df <- get_gsea_df(tcga_gsea_hallmarks) %>%
  select_gsea_cols(is_go = FALSE, suffix = "tcga")

geo_hm_df  <- get_gsea_df(geo_gsea_hallmarks) %>%
  select_gsea_cols(is_go = FALSE, suffix = "geo")

hallmarks_meta <- combine_two_gsea(
  df_geo  = geo_hm_df,
  df_tcga = tcga_hm_df,
  keys    = c("ID","Description")   # no ONTOLOGY en Hallmarks
)

# Nombres "bonitos" si los tienes en hallmark_names (ID -> label_ready)
filt_hallmarks_meta <- hallmarks_meta %>%
  dplyr::filter(!is.na(meta_pval),
                meta_padj < 0.05,
                concordant_sign) %>%
  dplyr::left_join(hallmark_names, by = "ID") %>%
  dplyr::mutate(label = dplyr::coalesce(label_ready, Description),
                full_label = factor(label, levels = unique(label))) %>%
  dplyr::arrange(desc(abs(mean_NES)), label)

readr::write_tsv(filt_hallmarks_meta, file.path(outputsFolder, "meta_gsea_HALLMARKS.tsv"))


hallm_top <- filt_hallmarks_meta %>%
  filter(is.finite(meta_padj), is.finite(mean_NES)) %>%
  arrange(desc(abs(mean_NES))) %>%
  slice_head(n = 30) %>%
  mutate(term = full_label) %>%  
  mutate(term = factor(term, levels = rev(unique(term))))

p_hallm <- ggplot(hallm_top, aes(x = mean_NES, y = term)) +
  geom_segment(aes(x = 0, xend = mean_NES, y = term, yend = term),
               linewidth = 1.0, color = "grey60") +
  geom_point(aes(color = meta_padj), size = 3) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.4, color = "grey50") +
  scale_color_gradient(low = "#2c7fb8", high = "#253494", trans = "reverse",
                       name = "Meta-adjusted p-value") +
  labs(title = "Meta-GSEA Hallmarks (TCGA + GEO)", x = "Mean NES", y = NULL) +
  theme_minimal(base_size = 13) +
  theme(panel.grid.major.y = element_blank(),
        panel.grid.minor   = element_blank(),
        axis.text.y        = element_text(size = 8),
        axis.text.x        = element_text(size = 8),
        legend.position    = "right",
        plot.title         = element_text(face = "bold", size = 15))

ggsave(file.path(plotsFolder, "meta_gsea_HALLMARKS_lollipop.pdf"), p_hallm, width = 16, height = 4, units = "in")
ggsave(file.path(plotsFolder, "meta_gsea_HALLMARKS_lollipop.png"), p_hallm, width = 16, height = 4, units = "in", dpi = 300)


#############################################
# 2. Load master regulator results and perform meta-analysis
#############################################

# Load VIPER MRA results for TCGA and GEO
mrs_table_tcga <- readRDS("rds/tcga_mrs_table.rds")
mrs_table_geo <- readRDS("rds/geo_mrs_table.rds")

# Ensure both have the same naming convention for TF column
colnames(mrs_table_tcga)[colnames(mrs_table_tcga) == "TF"] <- "tf"
colnames(mrs_table_geo)[colnames(mrs_table_geo) == "TF"] <- "tf"

# Inner join by TF
merged_mrs_tbl <- inner_join(mrs_table_tcga, mrs_table_geo, by = "tf", suffix = c("_tcga", "_geo"))

# Perform meta-analysis using sumlog on p-values
meta_results <- merged_mrs_tbl %>%
 rowwise() %>%
 mutate(
  meta_pval = sumlog(c(p.value_tcga, p.value_geo))$p,
  mean_nes = mean(c(nes_tcga, nes_geo)),
  concordant_sign = sign(nes_tcga) == sign(nes_geo)
 ) %>%
 ungroup() %>%
 mutate(meta_padj = p.adjust(meta_pval, method = "fdr")) %>%
 arrange(meta_padj)


sig_tmrs <- meta_results %>% filter(meta_padj < 0.05)

# Save meta-analysis results
write_tsv(meta_results, paste0(outputsFolder,"meta_mrs_results.tsv"))

# Volcano plot
pdf(paste0(plotsFolder, "meta_volcano_TMRs.pdf"), width = 8, height = 12)
EnhancedVolcano(meta_results,
  lab = meta_results$tf,
  x = 'mean_nes',
  y = 'meta_padj',
  pCutoff = 0.05,
  FCcutoff = 1,
  xlim = c(-4, 4),
  ylim = c(0, 3),
  xlab = bquote("Mean NES"),
  ylab = bquote(-Log[10]~italic(P)),
  title = "Meta–analysis Volcano Plot: TMRs",
  subtitle = "Fisher's combined p–values from TCGA and GEO",
  caption = "Only common TFs across both datasets",
  legendLabels = c("NS", "Mean NES", "p–value", "p–value and Mean NES"),
  legendPosition = "top",
  pointSize = 2.5,
  labSize = 3.5,
  colAlpha = 0.7,
  drawConnectors = TRUE,
  widthConnectors = 0.5,
  max.overlaps = 50
)
dev.off()

#############################################
# 3. Shadow significative TMRs
#############################################

# load VIPER objects
tcga_mrs <- readRDS("rds/tcga_mrs.rds")
geo_mrs <- readRDS("rds/geo_mrs.rds")

# Regulones por dataset
regulon_tcga <- tcga_mrs$regulon
regulon_geo <- geo_mrs$regulon

# TMRs validados por meta-análisis
tmrs_vector <- meta_results %>%
 filter(meta_padj < 0.05) %>%
 pull(tf)

tmr_shadow_df <- purrr::map_dfr(tmrs_vector, function(tf) {
 tryCatch({
  targets <- find_consensus_targets(tf, tmrs_vector, regulon_tcga, regulon_geo)
  if (length(targets) > 0) {
   tibble(regulator = tf, target_tmr = targets)
  } else {
   NULL
  }
 }, error = function(e) NULL)
})


g <- graph_from_data_frame(tmr_shadow_df, directed = TRUE)

pdf(file = paste0(plotsFolder, "autoregulated_TMR_network.pdf"), width = 8, height = 8)
plot(
 g,
 vertex.size = 10,       # tamaño de los nodos
 vertex.label.cex = .5,    # tamaño de la etiqueta (más grande)
 vertex.label.color = "blue", # color del texto
 vertex.color = "grey",   # color del nodo
 edge.arrow.size = 0.4,    # tamaño de la flecha (más pequeño)
 edge.color = "gray40",    # color del borde
 layout = layout_with_fr    # puedes cambiar el layout si quieres
)
dev.off()

#############################################
# 4. Meta-Regulon
#############################################

# Built the combined regulon only for the validated TMRs.
meta_regulon <- purrr::compact(purrr::map(setNames(tmrs_vector, tmrs_vector), 
                     ~build_meta_regulon(.x, regulon_tcga, regulon_geo)))

met_reg_tbl <- regulon2sif_tbl(meta_regulon)

write_tsv(met_reg_tbl, file = paste0(outputsFolder, "meta_regulon.txt"),)

# Shared targets by TMR
meta_reg_stats <- met_reg_tbl %>% 
 count(source, sort = TRUE)

write_tsv(meta_reg_stats, paste0(outputsFolder, "shared_tarjets_by_tmr.tsv"))

#############################################
# 9. ORA: Hallmark enrichment with MSigDB gene sets
#############################################

# ----------------------------
# TCGA data
# ----------------------------
# Universe of expressed genes used in ORA
tcga_gene_universe <- read_tsv("inmat_TCGA.txt", col_names = TRUE)[[1]]

tcga_reg_df <- getregulon(regulon_tcga, tf_list = tmrs_vector)

# Run ORA mode 1: per TF
tcga_ora_hallmarks_byTF <- ora_hallmarks(regulon_df = tcga_reg_df,
            hallmark_gmt = hallmark_gmt,
            universe = tcga_gene_universe,
            mode = "by_tf")


# Run ORA mode 2: all targets combined
tcga_ora_hallmarks_all <- ora_hallmarks(regulon_df = tcga_reg_df,
             hallmark_gmt = hallmark_gmt,
             universe = tcga_gene_universe,
             mode = "all_targets")

# ----------------------------
# GEO data
# ----------------------------
# Universe of expressed genes used in ORA
geo_gene_universe <- read_tsv("inmat_GEO.txt", col_names = TRUE)[[1]]

geo_reg_df <- getregulon(regulon_geo, tf_list = tmrs_vector)

# Run ORA mode 1: per TF
geo_ora_hallmarks_byTF <- ora_hallmarks(regulon_df = geo_reg_df,
            hallmark_gmt = hallmark_gmt,
            universe = geo_gene_universe,
            mode = "by_tf")


# Run ORA mode 2: all targets combined
geo_ora_hallmarks_all <- ora_hallmarks(regulon_df = geo_reg_df,
             hallmark_gmt = hallmark_gmt,
             universe = geo_gene_universe,
             mode = "all_targets")

#############################################
# 10. Combine and perform the meta‑analysis by TF and Hallmark.
#############################################

# ----------------------------
# by TF
# ----------------------------

hallm_meta_byTF <- full_join(geo_ora_hallmarks_byTF, tcga_ora_hallmarks_byTF,
               by = c("set", "ID"),
               suffix = c("_geo", "_tcga")) %>%
 rowwise() %>%
 mutate(
  meta_pval = if (is.finite(pvalue_geo) & is.finite(pvalue_tcga)) {
   sumlog(c(pvalue_geo, pvalue_tcga))$p
  } else {
   NA_real_
  }
 ) %>%
 ungroup() %>% 
 mutate(
  meta_padj = p.adjust(meta_pval, method = "fdr"),
  log10_meta_pval = -log10(meta_pval)
 ) %>%
 left_join(hallmark_names, by = "ID") %>%
 rename(label = label_ready) %>%
 left_join(sig_tmrs %>% select(tf, mean_nes), by = c("set" = "tf")) %>%
 arrange(desc(abs(mean_nes)))


filt_hallm_meta_byTF <- hallm_meta_byTF %>%
 filter(meta_padj < 0.05) %>%
 mutate(
  mean_enrichment = (FoldEnrichment_geo + FoldEnrichment_tcga) / 2
 ) %>%
 arrange(abs(mean_nes), desc(set), mean_enrichment) %>%
 mutate(
  full_label = paste(set, label, sep = ": "),
  full_label = factor(full_label, levels = full_label) 
 )

write_tsv(filt_hallm_meta_byTF, paste0(outputsFolder, "meta_hallm_byTF.tsv"))

######################################################################################


meta_hallm_byTF_barplot <- ggplot(filt_hallm_meta_byTF, aes(x = full_label, y = mean_enrichment, fill = meta_padj)) +
 geom_col(width = 0.7, color = "grey20") +
 coord_flip() +
 scale_fill_gradient(
  low = "#2c7fb8", high = "#253494",
  name = "Meta-adjusted\np-value",
  trans = "reverse"
 ) +
 labs(
  title = "Hallmark pathway enrichment across TMR regulons",
  subtitle = "Pathways with meta-adjusted p-value < 0.05",
  x = NULL,
  y = "Average fold enrichment",
  caption = "Data from GEO and TCGA"
 ) +
 theme_minimal(base_size = 14) +
 theme(
  axis.text.y = element_markdown(size = 12),
  axis.text.x = element_text(size = 12),
  legend.position = "right",
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_blank(),
  plot.title = element_text(face = "bold", size = 18),
  plot.subtitle = element_text(size = 13, margin = margin(b = 10)),
  plot.caption = element_text(size = 10, face = "italic")
 )

# LOLLIPOP para meta hallmarks por TF
meta_hallm_byTF_lollipop <- ggplot(
  filt_hallm_meta_byTF,
  aes(x = full_label, y = mean_enrichment)
) +
  geom_segment(aes(xend = full_label, y = 0, yend = mean_enrichment),
               linewidth = 0.9, color = "grey60") +
  geom_point(aes(color = meta_padj),
             size = 3.4, stroke = 0.2) +
  coord_flip() +
  scale_color_gradient(
    low = "#2c7fb8", high = "#253494",
    name = "Meta-adjusted\np-value",
    trans = "reverse"
  ) +
  labs(
    title = "Hallmark pathway enrichment across TMR regulons",
    subtitle = "Pathways with meta-adjusted p-value < 0.05",
    x = NULL,
    y = "Average fold enrichment",
    caption = "Data from GEO and TCGA"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_markdown(size = 8),
    axis.text.x = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 13, margin = margin(b = 10)),
    plot.caption = element_text(size = 10, face = "italic")
  )


######################################################################################

ggsave(paste0(plotsFolder, "meta_hallm_byTF_barplot.pdf"), plot = meta_hallm_byTF_barplot, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_hallm_byTF_barplot.png"), plot = meta_hallm_byTF_barplot, 
    width = 12, height = 10, units = "in", dpi = 300)

ggsave(paste0(plotsFolder, "meta_hallm_byTF_lollipop.pdf"), plot = meta_hallm_byTF_lollipop, 
  width = 12, height = 5, dpi = 300)
ggsave(paste0(plotsFolder, "meta_hallm_byTF_lollipop.png"), plot = meta_hallm_byTF_lollipop, 
  width = 12, height = 5, dpi = 300)


# ----------------------------
### All
# ----------------------------

hallm_meta_all <- full_join(geo_ora_hallmarks_all, tcga_ora_hallmarks_all,
           by = c("set", "ID"),
           suffix = c("_geo", "_tcga")) %>%
 rowwise() %>%
 mutate(
  meta_pval = if (is.finite(pvalue_geo) & is.finite(pvalue_tcga)) {
   sumlog(c(pvalue_geo, pvalue_tcga))$p
  } else {
   NA_real_
  }
 ) %>%
 ungroup() %>% 
 mutate(
  meta_padj = p.adjust(meta_pval, method = "fdr"),
  log10_meta_pval = -log10(meta_pval)
 ) %>%
 left_join(hallmark_names, by = "ID") %>%
 rename(label = label_ready) %>%
 left_join(sig_tmrs %>% select(tf, mean_nes), by = c("set" = "tf")) %>%
 arrange(desc(abs(mean_nes)))


filt_hallm_meta_all <- hallm_meta_all %>%
 filter(meta_padj < 0.05) %>%
 mutate(
  mean_enrichment = (FoldEnrichment_geo + FoldEnrichment_tcga) / 2
 ) %>%
 arrange(abs(mean_nes), desc(set), mean_enrichment) %>%
 mutate(
  full_label = paste(set, label, sep = ": "),
  full_label = factor(full_label, levels = full_label) 
 )

write_tsv(filt_hallm_meta_all, paste0(outputsFolder, "meta_hallm_all.tsv"))
######################################################################################


meta_hallm_all_barplot <- ggplot(filt_hallm_meta_all, aes(x = full_label, y = mean_enrichment, fill = meta_padj)) +
 geom_col(width = 0.7, color = "grey20") +
 coord_flip() +
 scale_fill_gradient(
  low = "#2c7fb8", high = "#253494",
  name = "Meta-adjusted\np-value",
  trans = "reverse"
 ) +
 labs(
  title = "Hallmark pathway enrichment across TMR regulons",
  subtitle = "Pathways with meta-adjusted p-value < 0.05",
  x = NULL,
  y = "Average fold enrichment",
  caption = "Data from GEO and TCGA"
 ) +
 theme_minimal(base_size = 14) +
 theme(
  axis.text.y = element_markdown(size = 12),
  axis.text.x = element_text(size = 12),
  legend.position = "right",
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_blank(),
  plot.title = element_text(face = "bold", size = 18),
  plot.subtitle = element_text(size = 13, margin = margin(b = 10)),
  plot.caption = element_text(size = 10, face = "italic")
 )

# LOLLIPOP para meta hallmarks (ALL)
meta_hallm_all_lollipop <- ggplot(
  filt_hallm_meta_all,
  aes(x = full_label, y = mean_enrichment)
) +
  # tallo
  geom_segment(aes(xend = full_label, y = 0, yend = mean_enrichment),
               linewidth = 0.9, color = "grey60") +
  geom_point(aes(color = meta_padj),
             size = 3.4, stroke = 0.2) +
  coord_flip() +
  scale_color_gradient(
    low = "#2c7fb8", high = "#253494",
    name = "Meta-adjusted\np-value",
    trans = "reverse"
  ) +
  labs(
    title = "Hallmark pathway enrichment across TMR regulons",
    subtitle = "Pathways with meta-adjusted p-value < 0.05",
    x = NULL,
    y = "Average fold enrichment",
    caption = "Data from GEO and TCGA"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_markdown(size = 8),
    axis.text.x = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 13, margin = margin(b = 10)),
    plot.caption = element_text(size = 10, face = "italic")
  )

######################################################################################

ggsave(paste0(plotsFolder, "meta_hallm_all_barplot.pdf"), plot = meta_hallm_all_barplot, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_hallm_all_barplot.png"), plot = meta_hallm_all_barplot, 
    width = 12, height = 10, units = "in", dpi = 300)

ggsave(paste0(plotsFolder, "meta_hallm_all_lollipop.pdf"),
       plot = meta_hallm_all_lollipop, width = 12, height = 4, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_hallm_all_lollipop.png"),
       plot = meta_hallm_all_lollipop, width = 12, height = 4, units = "in", dpi = 300)

# ----------------------------
### heatmap plot
# ----------------------------

# — 1. Prepare data —
tf_hallmark_long <- filt_hallm_meta_byTF %>% 
 distinct(set, label) %>% 
 mutate(present = 1) %>% 
 pivot_wider(names_from = label, values_from = present, values_fill = 0) %>% 
 pivot_longer(-set, names_to = "hallmark", values_to = "present") %>% 
 filter(present == 1)

hallmark_counts <- tf_hallmark_long %>% 
 count(hallmark) %>% 
 arrange(desc(n))

tf_counts <- tf_hallmark_long %>% 
 count(set) %>% 
 arrange(n)

hallmark_levels <- hallmark_counts$hallmark
tf_levels    <- tf_counts$set

tf_hallmark_long <- tf_hallmark_long %>% 
 mutate(
  hallmark = factor(hallmark, levels = hallmark_levels),
  set    = factor(set,    levels = tf_levels)
 )
hallmark_counts <- hallmark_counts %>% 
 mutate(hallmark = factor(hallmark, levels = hallmark_levels))
tf_counts    <- tf_counts    %>% 
 mutate(set    = factor(set, levels = tf_levels))

# — 2. Three panels —

# 2a) Dot matrix (heatmap‐like)
main_plot <- ggplot(tf_hallmark_long, aes(x = hallmark, y = set)) +
 geom_point(size = 3) +
 scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
 scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
 theme_minimal() +
 theme(
  axis.title = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1),
  panel.grid = element_blank()
 )

# 2b) Bars above with mixed labels
top_bar <- ggplot(hallmark_counts, aes(x = hallmark, y = n)) +
 geom_col(fill = "grey40", width = 0.4) +
 geom_text(
  data  = subset(hallmark_counts, n >= 4),
  aes(label = n),
  vjust  = 2, 
  color  = "white",
  size  = 3
 ) +
 geom_text(
  data  = subset(hallmark_counts, n <= 3),
  aes(label = n),
  vjust  = -0.3, 
  color  = "black",
  size  = 3
 ) +
 scale_x_discrete(limits = hallmark_levels, expand = c(0, 0)) +
 scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
 theme_minimal() +
 theme(
  axis.title = element_blank(),
  axis.text  = element_blank(),
  axis.ticks = element_blank(),
  panel.grid = element_blank()
 )

# 2c) Bars on the right with conditional labels
right_bar <- ggplot(tf_counts, aes(x = n, y = set)) +
 geom_col(fill = "grey40", width = 0.6) +
 geom_text(
  data = filter(tf_counts, n >= 2),
  aes(label = n),
  hjust = 2,
  color = "white",
  size = 3
 ) +
 geom_text(
  data = filter(tf_counts, n == 1),
  aes(label = n),
  hjust = -0.2,
  color = "black",
  size = 3
 ) +
 scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
 scale_y_discrete(limits = tf_levels, expand = c(0, 0)) +
 theme_minimal() +
 theme(
  axis.title  = element_blank(),
  axis.text  = element_blank(),
  axis.ticks  = element_blank(),
  panel.grid  = element_blank()
 )

# 3. Assemble it
upper <- cowplot::plot_grid(
 NULL, top_bar, NULL,
 ncol    = 3,
 rel_widths = c(0.3, 3.9, 1.15),
 align   = "h"
)

lower <- cowplot::plot_grid(
 main_plot, right_bar,
 ncol    = 2,
 rel_widths = c(4, 1),
 align   = "h"
)

co_ocurrence <- cowplot::plot_grid(
 upper, lower,
 nrow    = 2,
 rel_heights = c(0.5, 4),
 align    = "v"
)

co_ocurrence

ggsave(paste0(plotsFolder,"meta_tmrs-hallamrks_co-ocurrence-plot.pdf"), plot = co_ocurrence, 
    width = 15, height = 9, units = "in", dpi = 300)
ggsave(paste0(plotsFolder,"meta_tmrs-hallamrks_co-ocurrence-plot.png"), plot = co_ocurrence, 
    width = 15, height = 9, units = "in", dpi = 300)

#############################################
# 11. ORA: GO + KEGG enrichment per TF and all targets
#############################################

# ----------------------------
# TCGA data
# ----------------------------

# ORA per TF with GO (ALL categories)
tcga_ora_go_by_tf <- ora_go(regulon_df = tcga_reg_df,
            universe = tcga_gene_universe,
            mode = "by_tf",
            organism = org.Hs.eg.db)

# ORA with GO for all targets combined
tcga_ora_go_all_targets <- ora_go(regulon_df = tcga_reg_df,
               universe = tcga_gene_universe,
               mode = "all_targets",
               organism = org.Hs.eg.db)

# ORA per TF with KEGG
tcga_ora_kegg_by_tf <- ora_kegg(regulon_df = tcga_reg_df,
              universe = tcga_gene_universe,
              mode = "by_tf")

# ORA with KEGG for all targets combined
tcga_ora_kegg_all_targets <- ora_kegg(regulon_df = tcga_reg_df,
                 universe = tcga_gene_universe,
                 mode = "all_targets")


# ----------------------------
# GEO data
# ----------------------------

# ORA per TF with GO (ALL categories)
geo_ora_go_by_tf <- ora_go(regulon_df = geo_reg_df,
            universe = geo_gene_universe,
            mode = "by_tf",
            organism = org.Hs.eg.db)

# ORA with GO for all targets combined
geo_ora_go_all_targets <- ora_go(regulon_df = geo_reg_df,
               universe = geo_gene_universe,
               mode = "all_targets",
               organism = org.Hs.eg.db)

# ORA per TF with KEGG
geo_ora_kegg_by_tf <- ora_kegg(regulon_df = geo_reg_df,
              universe = geo_gene_universe,
              mode = "by_tf")

# ORA with KEGG for all targets combined
geo_ora_kegg_all_targets <- ora_kegg(regulon_df = geo_reg_df,
                 universe = geo_gene_universe,
                 mode = "all_targets")

#############################################
# 12. Meta analysis ORA: GO enrichment per TF 
#############################################

go_meta_byTF <- full_join(
  geo_ora_go_by_tf,
  tcga_ora_go_by_tf,
  by   = c("set", "ONTOLOGY", "ID", "Description"),
  suffix = c("_geo", "_tcga")
 ) %>%
 rowwise() %>%
 mutate(
  meta_pval = if (is.finite(pvalue_geo) && is.finite(pvalue_tcga)) {
   sumlog(c(pvalue_geo, pvalue_tcga))$p
  } else {
   NA_real_
  }
 ) %>%
 ungroup() %>%
 filter(!is.na(meta_pval)) %>%
 mutate(
  meta_padj   = p.adjust(meta_pval, method = "fdr"),
  log10_meta_pval = -log10(meta_pval)
 ) %>%
 left_join(
  sig_tmrs %>% dplyr::select(tf, mean_nes),
  by = c("set" = "tf")
 ) %>%
 arrange(desc(abs(mean_nes)))


filt_go_meta_byTF <- go_meta_byTF %>%
 filter(meta_padj < 0.05) %>%
 mutate(
  mean_enrichment = (FoldEnrichment_geo + FoldEnrichment_tcga) / 2
 ) %>%
 group_by(set, ONTOLOGY, ID, Description) %>%
 slice_min(order_by = meta_padj, n = 1, with_ties = FALSE) %>%
 ungroup() %>%
 arrange(abs(mean_nes), desc(set), mean_enrichment) %>%
 mutate(
  full_label = paste(set, ONTOLOGY, Description, sep = ": "),
  full_label = factor(full_label, levels = unique(full_label))
 )

write_tsv(filt_go_meta_byTF, paste0(outputsFolder, "meta_go_byTF.tsv"))
######################################################################################


meta_go_bar_plot_byTF <- ggplot(filt_go_meta_byTF, aes(x = full_label, y = mean_enrichment, fill = meta_padj)) +
 geom_col(width = 0.7, color = "grey20") +
 coord_flip() +
 scale_fill_gradient(
  low = "#2c7fb8", high = "#253494",
  name = "Meta-adjusted\np-value",
  trans = "reverse"
 ) +
 labs(
  title = "GO ontologies enrichment across TMR regulons",
  subtitle = "Top pathways with meta-adjusted p-value < 0.05",
  x = NULL,
  y = "Average fold enrichment",
  caption = "Data from GEO and TCGA"
 ) +
 theme_minimal(base_size = 14) +
 theme(
  axis.text.y = element_markdown(size = 12),
  axis.text.x = element_text(size = 12),
  legend.position = "right",
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_blank(),
  plot.title = element_text(face = "bold", size = 18),
  plot.subtitle = element_text(size = 13, margin = margin(b = 10)),
  plot.caption = element_text(size = 10, face = "italic")
 )

# LOLLIPOP para GO enrichment across TMR regulons
meta_go_lollipop_byTF <- ggplot(
  filt_go_meta_byTF,
  aes(x = full_label, y = mean_enrichment)
) +
  geom_segment(aes(xend = full_label, y = 0, yend = mean_enrichment),
               linewidth = 0.9, color = "grey60") +
  geom_point(aes(color = meta_padj),
             size = 3.4, stroke = 0.2) +
  coord_flip() +
  scale_color_gradient(
    low = "#2c7fb8", high = "#253494",
    name = "Meta-adjusted\np-value",
    trans = "reverse"
  ) +
  labs(
    title = "GO ontologies enrichment across TMR regulons",
    subtitle = "Top pathways with meta-adjusted p-value < 0.05",
    x = NULL,
    y = "Average fold enrichment",
    caption = "Data from GEO and TCGA"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_markdown(size = 8),
    axis.text.x = element_text(size = 10),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 13, margin = margin(b = 10)),
    plot.caption = element_text(size = 10, face = "italic")
  )


######################################################################################
ggsave(paste0(plotsFolder, "meta_go_bar_plot_byTF.pdf"), plot = meta_go_bar_plot_byTF, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_go_bar_plot_byTF.png"), plot = meta_go_bar_plot_byTF, 
    width = 12, height = 10, units = "in", dpi = 300)

ggsave(paste0(plotsFolder, "meta_go_lollipop_byTF.pdf"),
       plot = meta_go_lollipop_byTF, width = 12, height = 4, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_go_lollipop_byTF.png"),
       plot = meta_go_lollipop_byTF, width = 12, height = 4, units = "in", dpi = 300)


#############################################
# 13. Meta analysis ORA: KEGG enrichment per TF 
#############################################

kegg_meta_byTF <- full_join(
  geo_ora_kegg_by_tf,
  tcga_ora_kegg_by_tf,
  by   = c("set", "ID"),
  suffix = c("_geo", "_tcga")
 ) %>%
 mutate(
  category  = coalesce(category_geo,  category_tcga),
  subcategory = coalesce(subcategory_geo, subcategory_tcga),
  Description = coalesce(Description_geo, Description_tcga)
 ) %>%
 rowwise() %>%
 mutate(
  meta_pval = if (is.finite(pvalue_geo) && is.finite(pvalue_tcga)) {
   sumlog(c(pvalue_geo, pvalue_tcga))$p
  } else {
   NA_real_
  }
 ) %>%
 ungroup() %>%
 filter(!is.na(meta_pval)) %>%
 mutate(
  meta_padj   = p.adjust(meta_pval, method = "fdr"),
  log10_meta_pval = -log10(meta_pval)
 ) %>%
 left_join(
  sig_tmrs %>% dplyr::select(tf, mean_nes),
  by = c("set" = "tf")
 ) %>%
 arrange(desc(abs(mean_nes)))


filt_kegg_meta_byTF <- kegg_meta_byTF %>%
 filter(meta_padj < 0.05) %>%
 mutate(
  mean_enrichment = (FoldEnrichment_geo + FoldEnrichment_tcga) / 2
 ) %>%
 group_by(set, category, subcategory, ID, Description) %>%
 slice_min(order_by = meta_padj, n = 1, with_ties = FALSE) %>%
 ungroup() %>%
 arrange(abs(mean_nes), desc(set), mean_enrichment) %>%
 mutate(
  full_label = paste(set, Description, sep = ": "),
  full_label = factor(full_label, levels = unique(full_label))
 )

write_tsv(filt_kegg_meta_byTF, paste0(outputsFolder, "meta_kegg_byTF.tsv"))
######################################################################################

meta_kegg_bar_plot_byTF <- ggplot(filt_kegg_meta_byTF, aes(x = full_label, y = mean_enrichment, fill = meta_padj)) +
 geom_col(width = 0.7, color = "grey20") +
 coord_flip() +
 scale_fill_gradient(
  low = "#2c7fb8", high = "#253494",
  name = "Meta-adjusted\np-value",
  trans = "reverse"
 ) +
 labs(
  title = "KEGG phatways enrichment across TMR regulons",
  subtitle = "Top pathways with meta-adjusted p-value < 0.05",
  x = NULL,
  y = "Average fold enrichment",
  caption = "Data from GEO and TCGA"
 ) +
 theme_minimal(base_size = 14) +
 theme(
  axis.text.y = element_markdown(size = 12),
  axis.text.x = element_text(size = 12),
  legend.position = "right",
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_blank(),
  plot.title = element_text(face = "bold", size = 18),
  plot.subtitle = element_text(size = 13, margin = margin(b = 10)),
  plot.caption = element_text(size = 10, face = "italic")
 )

meta_kegg_lollipop_byTF <- ggplot(
  filt_kegg_meta_byTF,
  aes(x = full_label, y = mean_enrichment)
) +
  geom_segment(aes(xend = full_label, y = 0, yend = mean_enrichment),
               linewidth = 0.9, color = "grey60") +
  geom_point(aes(color = -log10(meta_padj), size = abs(mean_nes)),
             alpha = 0.95) +
  coord_flip() +
  scale_color_gradient(
    low = "#a6bddb", high = "#045a8d",
    name = expression(-log[10]("meta padj"))
  ) +
  scale_size_continuous(name = "|mean NES|") +
  labs(
    title = "KEGG pathways enrichment across TMR regulons",
    subtitle = "p.adj meta < 0.05 | color = -log10(p.adj), tamaño = |NES|",
    x = NULL,
    y = "Average fold enrichment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_markdown(size = 12),
    panel.grid.major.y = element_blank()
  )


######################################################################################

ggsave(paste0(plotsFolder, "meta_kegg_bar_plot_byTF.pdf"), plot = meta_kegg_bar_plot_byTF, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_kegg_bar_plot_byTF.png"), plot = meta_kegg_bar_plot_byTF, 
    width = 12, height = 10, units = "in", dpi = 300)

ggsave(paste0(plotsFolder, "meta_kegg_lollipop_byTF2.png"),
       plot = meta_kegg_lollipop_byTF2, width = 12, height = 10, units = "in", dpi = 300)


#############################################
# 14. Meta analysis ORA: GEO enrichment per all targets 
#############################################

# GEO did not have enriched ontologies

#############################################
# 15. Meta analysis ORA: KEGG enrichment per all targets 
#############################################

kegg_meta_all_targets <- full_join(
  geo_ora_kegg_all_targets,
  tcga_ora_kegg_all_targets,
  by   = c("set", "ID"),
  suffix = c("_geo", "_tcga")
 ) %>%
 mutate(
  category  = coalesce(category_geo,  category_tcga),
  subcategory = coalesce(subcategory_geo, subcategory_tcga),
  Description = coalesce(Description_geo, Description_tcga)
 ) %>%
 rowwise() %>%
 mutate(
  meta_pval = if (is.finite(pvalue_geo) && is.finite(pvalue_tcga)) {
   sumlog(c(pvalue_geo, pvalue_tcga))$p
  } else {
   NA_real_
  }
 ) %>%
 ungroup() %>%
 filter(!is.na(meta_pval)) %>%
 mutate(
  meta_padj   = p.adjust(meta_pval, method = "fdr"),
  log10_meta_pval = -log10(meta_pval)
 ) %>%
 left_join(
  sig_tmrs %>% dplyr::select(tf, mean_nes),
  by = c("set" = "tf")
 ) %>%
 arrange(desc(abs(mean_nes)))


filt_kegg_meta_all_targets <- kegg_meta_all_targets %>%
 filter(meta_padj < 0.05) %>%
 mutate(
  mean_enrichment = (FoldEnrichment_geo + FoldEnrichment_tcga) / 2
 ) %>%
 group_by(set, category, subcategory, ID, Description) %>%
 slice_min(order_by = meta_padj, n = 1, with_ties = FALSE) %>%
 ungroup() %>%
 arrange(abs(mean_nes), desc(set), mean_enrichment) %>%
 mutate(
  full_label = paste(set, Description, sep = ": "),
  full_label = factor(full_label, levels = unique(full_label))
 )

write_tsv(filt_kegg_meta_all_targets, paste0(outputsFolder, "meta_kegg_all.tsv"))
######################################################################################

meta_kegg_bar_plot_all_targets <- ggplot(filt_kegg_meta_all_targets, aes(x = full_label, y = mean_enrichment, fill = meta_padj)) +
 geom_col(width = 0.7, color = "grey20") +
 coord_flip() +
 scale_fill_gradient(
  low = "#2c7fb8", high = "#253494",
  name = "Meta-adjusted\np-value",
  trans = "reverse"
 ) +
 labs(
  title = "KEGG phatways enrichment across TMR regulons",
  subtitle = "Top pathways with meta-adjusted p-value < 0.05",
  x = NULL,
  y = "Average fold enrichment",
  caption = "Data from GEO and TCGA"
 ) +
 theme_minimal(base_size = 14) +
 theme(
  axis.text.y = element_markdown(size = 12),
  axis.text.x = element_text(size = 12),
  legend.position = "right",
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_blank(),
  plot.title = element_text(face = "bold", size = 18),
  plot.subtitle = element_text(size = 13, margin = margin(b = 10)),
  plot.caption = element_text(size = 10, face = "italic")
 )

meta_kegg_lollipop_all_targets <- ggplot(
  filt_kegg_meta_all_targets,
  aes(x = full_label, y = mean_enrichment)
) +
  geom_segment(aes(xend = full_label, y = 0, yend = mean_enrichment),
               linewidth = 0.9, color = "grey60") +
  geom_point(aes(color = -log10(meta_padj), size = abs(mean_nes)),
             alpha = 0.95) +
  coord_flip() +
  scale_color_gradient(
    low = "#a6bddb", high = "#045a8d",
    name = expression(-log[10]("meta padj"))
  ) +
  scale_size_continuous(name = "|mean NES|") +
  labs(
    title = "KEGG pathways enrichment across TMR regulons",
    subtitle = "p.adj meta < 0.05 | color = -log10(p.adj), tamaño = |NES|",
    x = NULL,
    y = "Average fold enrichment"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_markdown(size = 12),
    panel.grid.major.y = element_blank()
  )

######################################################################################

ggsave(paste0(plotsFolder, "meta_kegg_bar_plot_all_targets.pdf"), plot = meta_kegg_bar_plot_all_targets, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_kegg_bar_plot_all_targets.png"), plot = meta_kegg_bar_plot_all_targets, 
    width = 12, height = 10, units = "in", dpi = 300)

ggsave(paste0(plotsFolder, "meta_kegg_lollipop_all_targets2.png"),
       plot = meta_kegg_lollipop_all_targets2, width = 12, height = 10, units = "in", dpi = 300)


# ==== FIGURE 4 ====


p4a_kegg_all <- ggplot(
  filt_kegg_meta_all_targets,
  aes(x = full_label, y = mean_enrichment)
) +
  geom_segment(aes(xend = full_label, y = 0, yend = mean_enrichment),
               linewidth = 0.9, color = "grey60") +
  geom_point(aes(color = -log10(meta_padj)),
             size = 3.4, stroke = 0.2, alpha = 0.95) +
  coord_flip() +
  scale_color_gradient(
    low = "#a6bddb", high = "#045a8d",
    name = expression(-log[10]("meta padj"))
  ) +
  labs(
    title = "(a) KEGG pathways",
    subtitle = "p.adj meta < 0.05 (Fisher + FDR) | color = -log10(p.adj)",
    x = NULL, y = "Average fold enrichment",
    caption = "Data from GEO and TCGA"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_markdown(size = 8),
    panel.grid.major.y = element_blank()
  )

p4b_hallm_all <- ggplot(
  filt_hallm_meta_all,
  aes(x = full_label, y = mean_enrichment)
) +
  geom_segment(aes(xend = full_label, y = 0, yend = mean_enrichment),
               linewidth = 0.9, color = "grey60") +
  geom_point(aes(color = -log10(meta_padj)),
             size = 3.4, stroke = 0.2, alpha = 0.95) +
  coord_flip() +
  scale_color_gradient(
    low = "#a6bddb", high = "#045a8d",
    name = expression(-log[10]("meta padj"))
  ) +
  labs(
    title = "(b) MSigDB Hallmarks",
    subtitle = "p.adj meta < 0.05 (Fisher + FDR) | color = -log10(p.adj)",
    x = NULL, y = "Average fold enrichment",
    caption = "No significant GO enrichments were detected"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_markdown(size = 8),
    panel.grid.major.y = element_blank()
  )

figure4 <- (p4a_kegg_all / p4b_hallm_all) +
  plot_annotation(theme = theme(
      plot.title = element_text(face = "bold", size = 18, margin = margin(b = 8))
    )
  )

ggsave(paste0(plotsFolder, "Figure4_meta_enrichment_combined_targets.pdf"),
       plot = figure4, width = 15, height = 5, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "Figure4_meta_enrichment_combined_targets.png"),
       plot = figure4, width = 15, height = 5, units = "in", dpi = 300)
 
# ==== FIGURE 5 ====


p5a_hallm_byTF <- ggplot(
  filt_hallm_meta_byTF,
  aes(x = full_label, y = mean_enrichment)
) +
  geom_segment(aes(xend = full_label, y = 0, yend = mean_enrichment),
               linewidth = 0.9, color = "grey60") +
  geom_point(aes(color = -log10(meta_padj), size = abs(mean_nes)),
             alpha = 0.95) +
  coord_flip() +
  scale_color_gradient(low = "#a6bddb", high = "#045a8d",
                       name = expression(-log[10]("meta padj"))) +
  scale_size_continuous(name = "|mean NES|") +
  labs(title = "(A) Hallmark gene sets",
       subtitle = "Consistent links to estrogen response, inflammatory pathways, and EMT",
       x = NULL, y = "Average fold enrichment") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_markdown(size = 10),
    panel.grid.major.y = element_blank()
  )

p5b_kegg_byTF <- ggplot(
  filt_kegg_meta_byTF,
  aes(x = full_label, y = mean_enrichment)
) +
  geom_segment(aes(xend = full_label, y = 0, yend = mean_enrichment),
               linewidth = 0.9, color = "grey60") +
  geom_point(aes(color = -log10(meta_padj), size = abs(mean_nes)),
             alpha = 0.95) +
  coord_flip() +
  scale_color_gradient(low = "#a6bddb", high = "#045a8d",
                       name = expression(-log[10]("meta padj"))) +
  scale_size_continuous(name = "|mean NES|") +
  labs(title = "(B) KEGG pathways",
       subtitle = "TGF-β, TNF, adipocytokine, antigen presentation, cytokine interactions, thyroid hormone synthesis",
       x = NULL, y = "Average fold enrichment") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_markdown(size = 10),
    panel.grid.major.y = element_blank()
  )

p5c_go_byTF <- ggplot(
  filt_go_meta_byTF,
  aes(x = full_label, y = mean_enrichment)
) +
  geom_segment(aes(xend = full_label, y = 0, yend = mean_enrichment),
               linewidth = 0.9, color = "grey60") +
  geom_point(aes(color = -log10(meta_padj), size = abs(mean_nes)),
             alpha = 0.95) +
  coord_flip() +
  scale_color_gradient(low = "#a6bddb", high = "#045a8d",
                       name = expression(-log[10]("meta padj"))) +
  scale_size_continuous(name = "|mean NES|") +
  labs(title = "(C) Gene Ontology (GO)",
       subtitle = "Mainly associated with RUNX2 regulons (actin filaments, ECM components)",
       x = NULL, y = "Average fold enrichment") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_markdown(size = 10),
    panel.grid.major.y = element_blank()
  )

figure5 <- (p5a_hallm_byTF / p5b_kegg_byTF / p5c_go_byTF) +
  plot_layout(heights = c(3, 2, 0.6)) +  # A > B > C
  plot_annotation(theme = theme())

ggsave(paste0(plotsFolder, "Figure5_meta_enrichment_TMR_regulons.pdf"),
       plot = figure5, width = 15, height = 16, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "Figure5_meta_enrichment_TMR_regulons.png"),
       plot = figure5, width = 15, height = 16, units = "in", dpi = 300)

# ==== FIGURA 6 (versión mejorada)  ============================================
# ----------------------------
### heatmap plot
# ----------------------------

# — 1. Prepare data —
tf_hallmark_long <- filt_hallm_meta_byTF %>% 
 distinct(set, label) %>% 
 mutate(present = 1) %>% 
 pivot_wider(names_from = label, values_from = present, values_fill = 0) %>% 
 pivot_longer(-set, names_to = "hallmark", values_to = "present") %>% 
 filter(present == 1)

hallmark_counts <- tf_hallmark_long %>% 
 count(hallmark) %>% 
 arrange(desc(n))

tf_counts <- tf_hallmark_long %>% 
 count(set) %>% 
 arrange(n)

hallmark_levels <- hallmark_counts$hallmark
tf_levels    <- tf_counts$set

tf_hallmark_long <- tf_hallmark_long %>% 
 mutate(
  hallmark = factor(hallmark, levels = hallmark_levels),
  set    = factor(set,    levels = tf_levels)
 )
hallmark_counts <- hallmark_counts %>% 
 mutate(hallmark = factor(hallmark, levels = hallmark_levels))
tf_counts    <- tf_counts    %>% 
 mutate(set    = factor(set, levels = tf_levels))

# Consolida por (set, hallmark): color = -log10(padj) más extremo; size = promedio |NES|
df_bubbles <- filt_hallm_meta_byTF %>%
  dplyr::transmute(
    set,
    hallmark = label,
    logp  = -log10(meta_padj),
    nes   = mean_nes
  ) %>%
  dplyr::group_by(set, hallmark) %>%
  dplyr::summarise(
    logp = max(logp, na.rm = TRUE),
    size = mean(abs(nes), na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    hallmark = factor(hallmark, levels = hallmark_levels),
    set      = factor(set,      levels = tf_levels)
  )

# — 2. Three panels —

# 2a) Dot matrix (heatmap‐like)

# =========================
# 1) Panel central (burbujas)
# =========================
# df_bubbles ya debe existir y tener: set, hallmark, logp (-log10 padj), size (=|NES| promedio)
# Asegura el mismo orden que tus barras:
df_bubbles <- df_bubbles %>%
  mutate(
    hallmark = factor(hallmark, levels = hallmark_levels),
    set      = factor(set,      levels = tf_levels)
  )

main_plot <- ggplot(df_bubbles, aes(x = hallmark, y = set)) +
  geom_point(
    aes(fill = logp, size = size),
    shape = 21, color = "grey20", stroke = 0.25, alpha = 0.95
  ) +
  scale_fill_gradient(
    low  = "#a6bddb",
    high = "#045a8d",
    name = expression(-log[10]("meta padj"))
  ) +
  scale_size_continuous(
    name  = "Mean NES",
    range = c(2, 6)
  ) +
  scale_x_discrete(expand = expansion(add = c(0.5, 0.5))) +
  scale_y_discrete(expand = expansion(add = c(0.5, 0.5))) +
  theme_minimal(base_size = 12) +
  theme(
    axis.title   = element_blank(),
    axis.text.x  = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid   = element_blank()
  )

# =========================
# 2) Barras superior e izquierda (tu mismo código)
# =========================
top_bar <- ggplot(hallmark_counts, aes(x = hallmark, y = n)) +
  geom_col(fill = "grey40", width = 0.4) +
  geom_text(
    data  = subset(hallmark_counts, n >= 4),
    aes(label = n),
    vjust  = 2, color = "white", size = 3
  ) +
  geom_text(
    data  = subset(hallmark_counts, n <= 3),
    aes(label = n),
    vjust  = -0.3, color = "black", size = 3
  ) +
  scale_x_discrete(limits = hallmark_levels, expand = c(0, 0)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text  = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

right_bar <- ggplot(tf_counts, aes(x = n, y = set)) +
  geom_col(fill = "grey40", width = 0.6) +
  geom_text(
    data = dplyr::filter(tf_counts, n >= 2),
    aes(label = n),
    hjust = 2, color = "white", size = 3
  ) +
  geom_text(
    data = dplyr::filter(tf_counts, n == 1),
    aes(label = n),
    hjust = -0.2, color = "black", size = 3
  ) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  scale_y_discrete(limits = tf_levels, expand = c(0, 0)) +
  theme_minimal() +
  theme(
    axis.title  = element_blank(),
    axis.text   = element_blank(),
    axis.ticks  = element_blank(),
    panel.grid  = element_blank(),
    plot.margin = margin(r = 6)
  )

# =========================
# 3) Extraer leyenda y ensamblar con cowplot
# =========================

# Leyenda horizontal (colecta fill + size)
legend_grob <- cowplot::get_legend(
  main_plot +
    theme(
      legend.position   = "bottom",
      legend.direction  = "horizontal",
      legend.box        = "horizontal",
      legend.title      = element_text(size = 10),
      legend.text       = element_text(size = 9)
    ) +
    guides(
      fill = guide_colorbar(order = 1, barwidth = unit(120, "pt"), barheight = unit(6, "pt")),
      size = guide_legend(order = 2, override.aes = list(alpha = 1), direction = "horizontal")
    )
)

# Panel central sin leyenda
main_plot_noleg <- main_plot + theme(legend.position = "none")

# Fila superior (barras arriba)
upper <- cowplot::plot_grid(
  NULL, top_bar, NULL,
  ncol = 3,
  rel_widths = c(0.3, 3.9, 1.15),
  align = "h"
)

# Fila intermedia (matriz + barras derecha)
middle <- cowplot::plot_grid(
  main_plot_noleg, right_bar,
  ncol = 2,
  rel_widths = c(4, 1),
  align = "h"
)

# Fila inferior (leyenda)
bottom <- cowplot::plot_grid(
  legend_grob, ncol = 1
)

# Ensamble final: arriba / medio / leyenda
co_ocurrence <- cowplot::plot_grid(
  upper, middle, bottom,
  ncol = 1,
  rel_heights = c(0.5, 4, 0.45),  # ajusta la altura de la leyenda aquí
  align = "v"
)

# (opcional) guardar
# ggsave(paste0(plotsFolder, "Figure6_hallmarks_byTF_cooccurrence_with_legend_bottom.png"),
#        co_ocurrence, width = 16, height = 9, units = "in", dpi = 300)

co_ocurrence

ggsave(paste0(plotsFolder,"meta_tmrs-hallamrks_co-ocurrence-plot.pdf"), plot = co_ocurrence, 
    width = 15, height = 9, units = "in", dpi = 300)
ggsave(paste0(plotsFolder,"meta_tmrs-hallamrks_co-ocurrence-plot.png"), plot = co_ocurrence, 
    width = 15, height = 9, units = "in", dpi = 300)

# save.image("meta_analysis_session.RData")