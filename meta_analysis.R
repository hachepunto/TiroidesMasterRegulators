# meta_analysis.R
# Meta-analysis of master regulator results from TCGA and GEO
# Hugo Tovar, National Institute of Genomic Medicine, Mexico

#############################################
# 1. Load master regulator results and perform meta-analysis
#############################################

suppressPackageStartupMessages({
 library(tidyverse)
 library(metap)
 library(EnhancedVolcano)
 library(viper)
 library(igraph)
 library(clusterProfiler)
 library(ggtext)
 library(ComplexUpset)
 library(ggplot2)
 library(cowplot)
})

# Load helper functions
source("scripts/helpers.R")

# Create output directories
outputsFolder <- paste0(getwd(), "/meta_results/")
plotsFolder <- paste0(getwd(), "/meta_results/plots/")
dir.create(outputsFolder, showWarnings = FALSE)
dir.create(plotsFolder, showWarnings = FALSE)

#############################################
# 2. Join TCGA and GEO MRA results compute combined p-values using Fisher's method
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
  y = 'meta_pval',
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

# Cargar objetos VIPER con regulones completos
tcga_mrs <- readRDS("rds/tcga_mrs.rds")
geo_mrs <- readRDS("rds/geo_mrs.rds")

# Regulones por dataset
regulon_tcga <- tcga_mrs$regulon
regulon_geo <- geo_mrs$regulon

# TMRs validados por meta-análisis
tmrs_vector <- meta_results %>%
 filter(meta_padj < 0.05) %>%
 pull(tf)


# Aplicar a todos los TMRs validados por meta-análisis
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

# Construimos el regulón combinado solo para los TMRs validados
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

# Load Hallmark data (GSEA Broad GMT)
hallmark_gmt <- read.gmt("extra_data/h.all.v2025.1.Hs.symbols.gmt")
hallmark_names <- read_tsv("extra_data/hallmarks_names.txt")

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

meta_hallm_byTF_barplot

######################################################################################

meta_hallm_byTF_dotplot <- ggplot(filt_hallm_meta_byTF, aes(x = mean_enrichment, y = full_label)) +
 geom_point(aes(size = abs(mean_nes), color = meta_padj)) +
 scale_color_gradient(
  low = "#2c7fb8", high = "#253494",
  name = "Meta-adjusted\np-value",
  trans = "reverse"
 ) +
 scale_size_continuous(name = "Abs. Mean NES") +
 labs(
  title = "Hallmark pathway enrichment across TMR regulons",
  x = "Average fold enrichment",
  y = NULL
 ) +
 theme_minimal(base_size = 14) +
 theme(
  axis.text.y = element_text(size = 11),
  axis.text.x = element_text(size = 12),
  plot.title = element_text(face = "bold", size = 18),
  legend.position = "right"
 )

meta_hallm_byTF_dotplot

######################################################################################

ggsave(paste0(plotsFolder, "meta_hallm_byTF_barplot.pdf"), plot = meta_hallm_byTF_barplot, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_hallm_byTF_barplot.png"), plot = meta_hallm_byTF_barplot, 
    width = 12, height = 10, units = "in", dpi = 300)

ggsave(paste0(plotsFolder, "meta_hallm_byTF_dotplot.pdf"), plot = meta_hallm_byTF_dotplot, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_hallm_byTF_dotplot.png"), plot = meta_hallm_byTF_dotplot, 
    width = 12, height = 10, units = "in", dpi = 300)



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

meta_hallm_all_barplot

######################################################################################

meta_hallm_all_dotplot <- ggplot(filt_hallm_meta_all, aes(x = mean_enrichment, y = full_label)) +
 geom_point(aes(size = abs(mean_enrichment), color = meta_padj)) +
 scale_color_gradient(
  low = "#2c7fb8", high = "#253494",
  name = "Meta-adjusted\np-value",
  trans = "reverse"
 ) +
 scale_size_continuous(name = "Abs. Mean NES") +
 labs(
  title = "Hallmark pathway enrichment across TMR regulons",
  x = "Average fold enrichment",
  y = NULL
 ) +
 theme_minimal(base_size = 14) +
 theme(
  axis.text.y = element_text(size = 11),
  axis.text.x = element_text(size = 12),
  plot.title = element_text(face = "bold", size = 18),
  legend.position = "right"
 )

meta_hallm_all_dotplot

######################################################################################

ggsave(paste0(plotsFolder, "meta_hallm_all_barplot.pdf"), plot = meta_hallm_all_barplot, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_hallm_all_barplot.png"), plot = meta_hallm_all_barplot, 
    width = 12, height = 10, units = "in", dpi = 300)

ggsave(paste0(plotsFolder, "meta_hallm_all_dotplot.pdf"), plot = meta_hallm_all_dotplot, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_hallm_all_dotplot.png"), plot = meta_hallm_all_dotplot, 
    width = 12, height = 10, units = "in", dpi = 300)

# ----------------------------
### Upset plot
# ----------------------------


# Convertir a formato wide con presencia/ausencia
upset_hallmark_matrix <- filt_hallm_meta_byTF %>%
 select(set, label) %>% 
 mutate(value = 1) %>%
 pivot_wider(names_from = label, values_from = value, values_fill = 0)

# Obtener los nombres de las columnas de hallmark
upset_hallmarks <- setdiff(names(upset_hallmark_matrix), "set")

# Graficar
upset(upset_hallmark_matrix, 
   intersect = upset_hallmarks,
   name = "Shared Hallmarks",
   width_ratio = 0.2,
   min_size = 1,
   sort_sets = "descending",
   sort_intersections_by = "cardinality")

# ----------------------------
### heatmap plot
# ----------------------------

# — 1. Prepara tus datos igual que antes —
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

# — 2. Los tres paneles —

# 2a) Matriz de puntos (heatmap‐like)
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

# 2b) Barras arriba con etiquetas mixtas
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

# 2c) Barras a la derecha con etiquetas condicionadas
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

# 3. Ensámblalo como antes
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
meta_go_bar_plot_byTF
######################################################################################
meta_go_dot_plot_byTF <- ggplot(filt_go_meta_byTF, aes(x = mean_enrichment, y = full_label)) +
 geom_point(aes(size = mean_enrichment, color = meta_padj)) +
 scale_color_gradient(
  low = "#2c7fb8", high = "#253494",
  name = "Meta-adjusted\np-value",
  trans = "reverse"
 ) +
 scale_size_continuous(name = "Mean enrichment") +
 labs(
  title = "GO ontologies enrichment across TMR regulons",
  x = "Average fold enrichment",
  y = NULL
 ) +
 theme_minimal(base_size = 14) +
 theme(
  axis.text.y = element_text(size = 11),
  axis.text.x = element_text(size = 12),
  plot.title = element_text(face = "bold", size = 18),
  legend.position = "right"
 )
meta_go_dot_plot_byTF
######################################################################################
ggsave(paste0(plotsFolder, "meta_go_dot_plot_byTF.pdf"), plot = meta_go_dot_plot_byTF, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_go_dot_plot_byTF.png"), plot = meta_go_dot_plot_byTF, 
    width = 12, height = 10, units = "in", dpi = 300)

ggsave(paste0(plotsFolder, "meta_go_bar_plot_byTF.pdf"), plot = meta_go_bar_plot_byTF, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_go_bar_plot_byTF.png"), plot = meta_go_bar_plot_byTF, 
    width = 12, height = 10, units = "in", dpi = 300)


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
meta_kegg_bar_plot_byTF
######################################################################################
meta_kegg_dot_plot_byTF <- ggplot(filt_kegg_meta_byTF, aes(x = mean_enrichment, y = full_label)) +
 geom_point(aes(size = mean_enrichment, color = meta_padj)) +
 scale_color_gradient(
  low = "#2c7fb8", high = "#253494",
  name = "Meta-adjusted\np-value",
  trans = "reverse"
 ) +
 scale_size_continuous(name = "Mean enrichment") +
 labs(
  title = "KEGG pathway enrichment across TMR regulons",
  x = "Average fold enrichment",
  y = NULL
 ) +
 theme_minimal(base_size = 14) +
 theme(
  axis.text.y = element_text(size = 11),
  axis.text.x = element_text(size = 12),
  plot.title = element_text(face = "bold", size = 18),
  legend.position = "right"
 )
meta_kegg_dot_plot_byTF
######################################################################################
ggsave(paste0(plotsFolder, "meta_kegg_dot_plot_byTF.pdf"), plot = meta_kegg_dot_plot_byTF, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_kegg_dot_plot_byTF.png"), plot = meta_kegg_dot_plot_byTF, 
    width = 12, height = 10, units = "in", dpi = 300)

ggsave(paste0(plotsFolder, "meta_kegg_bar_plot_byTF.pdf"), plot = meta_kegg_bar_plot_byTF, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_kegg_bar_plot_byTF.png"), plot = meta_kegg_bar_plot_byTF, 
    width = 12, height = 10, units = "in", dpi = 300)


#############################################
# 14. Meta analysis ORA: GEO enrichment per all targets 
#############################################

# GEO No tuvo ontologías enriquecidas

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
meta_kegg_bar_plot_all_targets
######################################################################################
meta_kegg_dot_plot_all_targets <- ggplot(filt_kegg_meta_all_targets, aes(x = mean_enrichment, y = full_label)) +
 geom_point(aes(size = mean_enrichment, color = meta_padj)) +
 scale_color_gradient(
  low = "#2c7fb8", high = "#253494",
  name = "Meta-adjusted\np-value",
  trans = "reverse"
 ) +
 scale_size_continuous(name = "Mean enrichment") +
 labs(
  title = "KEGG pathway enrichment across TMR regulons",
  x = "Average fold enrichment",
  y = NULL
 ) +
 theme_minimal(base_size = 14) +
 theme(
  axis.text.y = element_text(size = 11),
  axis.text.x = element_text(size = 12),
  plot.title = element_text(face = "bold", size = 18),
  legend.position = "right"
 )
meta_kegg_dot_plot_all_targets
######################################################################################
ggsave(paste0(plotsFolder, "meta_kegg_dot_plot_all_targets.pdf"), plot = meta_kegg_dot_plot_all_targets, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_kegg_dot_plot_all_targets.png"), plot = meta_kegg_dot_plot_all_targets, 
    width = 12, height = 10, units = "in", dpi = 300)

ggsave(paste0(plotsFolder, "meta_kegg_bar_plot_all_targets.pdf"), plot = meta_kegg_bar_plot_all_targets, 
    width = 12, height = 10, units = "in", dpi = 300)
ggsave(paste0(plotsFolder, "meta_kegg_bar_plot_all_targets.png"), plot = meta_kegg_bar_plot_all_targets, 
    width = 12, height = 10, units = "in", dpi = 300)



# save.image("meta_analysis_session.RData")