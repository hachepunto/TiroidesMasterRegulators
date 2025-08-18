
# tmr_motifs.R
# Transcriptomic analysis of TCGA-THCA for the subtype "Papillary adenocarcinoma, NOS"
# Hugo Tovar, National Inmstitute of Genomic Medicine, Mexico hatovar@inmegen.gob.mx


#############################################
# 1) Load required libraries and set up folders
#############################################

# Motif databases and tools
suppressPackageStartupMessages({
  library(JASPAR2022)
  library(TFBSTools)
  
  # Genomics utilities
  library(biomaRt)
  library(GenomicRanges)
  library(BSgenome.Hsapiens.UCSC.hg38)
  library(Biostrings)
  
  # Data wrangling & I/O (load only what's needed)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(tibble)
  
  # Visualization
  library(circlize)
  library(RColorBrewer)
  library(metap)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(grid)
})

# Create output directories (idempotent)
outputsFolder <- file.path(getwd(), "TF_motives_results")
plotsFolder   <- file.path(outputsFolder, "plots/")
dir.create(outputsFolder, showWarnings = FALSE, recursive = TRUE)
dir.create(plotsFolder,   showWarnings = FALSE, recursive = TRUE)

#############################################
# 2) Obtaining motifs
#############################################

# Load the list of significant TMRs (adjusted p < 0.05)
# Expect a TSV with at least columns: 'tf' and 'meta_padj'
meta_tmrs <- vroom::vroom("meta_results/meta_mrs_results.tsv", .name_repair = janitor::make_clean_names) 

sign_tmrs <- meta_tmrs %>%
 filter(meta_padj < 0.05) %>% 
 pull(tf)
 
# Helper: fetch the latest PFMatrix version for a given TF name from JASPAR
get_pfm_for_tf <- function(tf_name) {
  ms <- TFBSTools::getMatrixSet(
    JASPAR2022,
    opts = list(all_versions = TRUE, name = tf_name)
  )
  if (length(ms) == 0) {
    warning("No motif found in JASPAR for: ", tf_name)
    return(NULL)
  }
  # Choose the latest version by ID suffix (e.g., MA0593.1 -> 1)
  ids <- vapply(ms, TFBSTools::ID, character(1))
  vers <- suppressWarnings(as.numeric(sub(".*\\.(\\d+)$", "\\1", ids)))
  ms[[which.max(vers)]]
}

# Build a PFM list for TMRs present in JASPAR (by TF name)
pfm_list <- purrr::map(sign_tmrs, get_pfm_for_tf) %>%
  rlang::set_names(sign_tmrs) %>%
  purrr::compact()

found_tfs   <- names(pfm_list)
missing_tfs <- setdiff(sign_tmrs, found_tfs)

# (Documented manual addition)
# Some motifs were retrieved manually from JASPAR (web UI) because they did not
# resolve by TF name. We add them here by explicit JASPAR IDs.
# recovered <- vroom::vroom("tfs_found_jaspar.txt")

manual_ids <- c("MA0145.2", "UN0262.1", "MA0040.1", "MA0117.2", "MA0840.1", "MA0643.1", "UN0333.1")

manual_pfms_core <- getMatrixSet(
  JASPAR2022,
  opts = list(ID = manual_ids)
)

# Coerce PFMatrixList to a plain list and name by the intended TF symbols
manual_pfms <- as.list(manual_pfms_core)
manual_pfms <- set_names(
  manual_pfms,
  recovered$tf[ match(names(manual_pfms), recovered$id_jaspar) ]
)

# Merge automatic and manual PFMs
# If a TF appears in both, prefer the manually specified entry.
pfm_list_extended <- c(
  pfm_list[setdiff(names(pfm_list), names(manual_pfms))],
  manual_pfms
)

# Convert PFMs to PWMs
pwm_list_extended <- map(pfm_list_extended, toPWM, pseudocounts = 0.1)

# Build a compact metadata table for bookkeeping and reproducibility
metadata_table <- purrr::imap_dfr(pwm_list_extended, function(pwm, tf) {
  tibble::tibble(
    tf         = tf,
    id_jaspar  = TFBSTools::ID(pwm),
    type       = pwm@tags$type %||% NA_character_,
    collection = pwm@tags$collection %||% NA_character_,
    species    = if (!is.null(pwm@tags$species)) paste(pwm@tags$species, collapse = ", ") else NA_character_
  )
}) %>%
  dplyr::arrange(tf)

readr::write_tsv(metadata_table, file.path(outputsFolder, "motif_metadata.tsv"))
print(metadata_table, n = nrow(metadata_table))
# # A tibble: 28 × 5
#    tf      id_jaspar type               collection  species
#    <chr>   <chr>     <chr>              <chr>       <chr>
#  1 FOXP2   MA0593.1  ChIP-seq           CORE        Homo sapiens
#  2 TCFL5   MA0632.2  HT-SELEX           CORE        Homo sapiens
#  3 TFCP2   MA1968.1  HT-SELEX           CORE        Homo sapiens
#  4 SREBF1  MA0595.1  ChIP-seq           CORE        Homo sapiens
#  5 PROX1   MA0794.1  HT-SELEX           CORE        Homo sapiens
#  6 ETV4    MA0764.1  HT-SELEX           CORE        Homo sapiens
#  7 TEAD4   MA0809.1  HT-SELEX           CORE        Homo sapiens
#  8 RARA    MA0729.1  HT-SELEX           CORE        Homo sapiens
#  9 RUNX2   MA0511.1  ChIP-seq           CORE        Homo sapiens
# 10 ETV1    MA0761.1  HT-SELEX           CORE        Homo sapiens
# 11 RXRG    MA0856.1  HT-SELEX           CORE        Homo sapiens
# 12 RARB    MA1552.1  HT-SELEX           CORE        Homo sapiens
# 13 BHLHE40 MA0464.2  HT-SELEX           CORE        Homo sapiens
# 14 HEY2    MA0649.1  HT-SELEX           CORE        Homo sapiens
# 15 ETV5    MA0765.1  HT-SELEX           CORE        Homo sapiens
# 16 GLIS3   MA0737.1  HT-SELEX           CORE        Homo sapiens
# 17 PLAG1   MA0163.1  bacterial 1-hybrid CORE        Homo sapiens
# 18 ZKSCAN3 MA1973.1  NA                 CORE        Homo sapiens
# 19 ZKSCAN5 MA1652.1  ChIP-seq           CORE        Homo sapiens
# 20 PKNOX2  MA0783.1  HT-SELEX           CORE        Homo sapiens
# 21 FOXE1   MA1487.1  HT-SELEX           CORE        Homo sapiens
# 22 TFCP2L1 MA0145.2  ChIP-seq           CORE        Mus musculus
# 23 SALL4   UN0262.1  ChIP-seq           UNVALIDATED Mus musculus
# 24 FOXQ1   MA0040.1  SELEX              CORE        Rattus norvegicus
# 25 MAFB    MA0117.2  HT-SELEX           CORE        Mus musculus
# 26 CREB5   MA0840.1  HT-SELEX           CORE        Mus musculus
# 27 ESRRG   MA0643.1  HT-SELEX           CORE        Mus musculus
# 28 ZNF548  UN0333.1  ChIP-seq           UNVALIDATED Homo sapiens

#############################################
# 3) Obtaining promoter regions
#############################################

# Query Ensembl for TSS data
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

tss_data <- getBM(
  attributes = c(
    "hgnc_symbol",
    "ensembl_transcript_id",
    "chromosome_name",
    "strand",
    "transcription_start_site"
  ),
  filters = "hgnc_symbol",
  values = sign_tmrs,
  mart = ensembl
)

# Report genes not found
not_found <- setdiff(sign_tmrs, unique(tss_data$hgnc_symbol))
if(length(not_found)) {
  warning("These genes were not found in Ensembl:\n", paste(not_found, collapse=", "))
}

# Select one representative TSS per gene
# Rule: for + strand take the minimal TSS (most upstream); for - strand take the maximal TSS.
# If a gene appears on both strands (rare), the first strand encountered is used and a warning is issued.
tss_upstream <- tss_data %>%
  group_by(hgnc_symbol) %>%
  dplyr::slice(if (unique(strand) == 1) which.min(transcription_start_site)
        else                 which.max(transcription_start_site)) %>%
  ungroup() %>%
  mutate(
    tss = transcription_start_site,
    tss_promoter_start = if_else(strand == 1, tss - 2000, tss - 200),
    tss_promoter_end   = if_else(strand == 1, tss + 200, tss + 2000)
  )

# Build GRanges for promoters (UCSC hg38 uses "chr" prefix)
promoter_gr <- GRanges(
  seqnames = paste0("chr", tss_upstream$chromosome_name),
  ranges   = IRanges(start = pmax(tss_upstream$tss_promoter_start, 1),
                     end   = tss_upstream$tss_promoter_end),
  strand   = ifelse(tss_upstream$strand == 1, "+", "-"),
  gene     = tss_upstream$hgnc_symbol
)

# Extract promoter sequences from BSgenome (hg38)
promoter_seqs <- getSeq(BSgenome.Hsapiens.UCSC.hg38, promoter_gr)

# Name sequences as TF_chr:start-end(strand)
names(promoter_seqs) <- paste0(
  tss_upstream$hgnc_symbol, "_",
  seqnames(promoter_gr), ":",
  start(promoter_gr), "-", end(promoter_gr),
  strand(promoter_gr)
)

# (Optional) Save promoters as FASTA for external use
# Biostrings::writeXStringSet(promoter_seqs, filepath = file.path(outputsFolder, "promoters_50tmrs_hg38.fa"))

#############################################
# 4) Motif hit search
#############################################

# Build all TF–target combinations
combos <- expand_grid(
  TF_from = names(pwm_list_extended),
  TF_to   = names(promoter_seqs)
)

# Search helper: run PWM scan for a (TF_from, TF_target) pair
search_pair <- function(TF_from, TF_to) {
  pwm <- pwm_list_extended[[TF_from]]
  seq <- promoter_seqs[[TF_to]]
  
  hits <- searchSeq(pwm, seq, min.score = "85%", strand = "*")
  if (length(hits) == 0) return(NULL)
  
  df <- as.data.frame(hits)
  df$TF_from <- TF_from
  df$TF_to   <- TF_to
  df
}

# 3. Execute scans over all pairs and bind results
binding_hits <- combos %>%
  pmap_dfr(~ search_pair(..1, ..2))

# binding_hits es un tibble con columnas:
#   start, end, strand, score, ... plus TF_from y TF_to

tss_up2 <- tss_upstream
colnames(tss_up2)[ colnames(tss_up2) == "strand" ] <- "strand_tss"

binding_hits2 <- binding_hits %>%
  mutate(
    TF_target = str_remove(TF_to, "_chr.*")
  )

# annotate genomic coordinates
binding_hits_annot <- binding_hits2 %>%
  left_join(
    tss_up2 %>% dplyr::select(hgnc_symbol, chromosome_name, strand_tss,
                       tss_promoter_start, tss_promoter_end),
    by = c("TF_target" = "hgnc_symbol")
  ) %>%
  mutate(
    genomic_start = if_else(strand_tss == 1L,
                            tss_promoter_start + start,
                            tss_promoter_end   - end),
    genomic_end   = if_else(strand_tss == 1L,
                            tss_promoter_start + end,
                            tss_promoter_end   - start),
    chromosome    = paste0("chr", chromosome_name)
  ) %>%
  dplyr::select(TF_from, TF_target, chromosome, genomic_start, genomic_end, everything())

# Save tabular outputs
readr::write_tsv(binding_hits_annot, paste0(outputsFolder,"binding_all50x28.tsv"))

#############################################
# 5) Summary table
#############################################

# 50 TFs list
tf_list <- tibble(tf = unique(c(binding_hits_annot$TF_from, binding_hits_annot$TF_target)))

# For each TF as emitter (TF_from): total hits and number of distinct targets
tf_from_summary <- binding_hits_annot %>%
  group_by(tf = TF_from) %>%
  summarise(
    n_hits    = n(),
    n_targets = n_distinct(TF_target),
    .groups   = "drop"
  )

# For each TF as target (TF_target): total incoming hits
tf_to_summary <- binding_hits_annot %>%
  group_by(tf = TF_target) %>%
  summarise(
    n_strikes = n(),
    .groups   = "drop"
  )

# "Beaters": number of distinct emitters that target each TF
beaters_summary <- binding_hits_annot %>%
  group_by(tf = TF_target) %>%
  summarise(
    n_beaters = n_distinct(TF_from),
    .groups   = "drop"
  )

tf_summary <- tf_list %>%
  left_join(tf_from_summary,   by = "tf") %>%
  left_join(tf_to_summary,     by = "tf") %>%
  left_join(beaters_summary,   by = "tf") %>%
  replace_na(list(
    n_hits    = 0,
    n_targets = 0,
    n_strikes = 0,
    n_beaters = 0
  ))

# Assemble final summary table
tf_summary <- tf_summary %>%
  left_join(
    tss_upstream %>%
      dplyr::select(
        hgnc_symbol,
        chromosome      = chromosome_name,
        genomic_start  = tss_promoter_start,
        genomic_end    = tss_promoter_end,
        strand          = strand
      ),
    by = c("tf" = "hgnc_symbol")
  ) %>%
  dplyr::select(
    tf,
    chromosome,
    genomic_start,
    genomic_end,
    strand,
    n_hits,
    n_targets,
    n_strikes,
    n_beaters
  ) %>% 
  arrange(desc(n_targets), desc(n_hits))

# Save and preview
print(tf_summary, n = 50)
readr::write_tsv(tf_summary, paste0(outputsFolder, "interactions_summary.tsv"))

#############################################
# 6) Circos plot
#############################################

# Build interaction matrix in the exact TF order used in the summary
tf_order <- tf_summary$tf

interaction_matrix <- binding_hits_annot %>%
  count(TF_from, TF_target, name = "n") %>%
  complete(TF_from = tf_order, TF_target = tf_order, fill = list(n = 0)) %>%
  pivot_wider(names_from = TF_target, values_from = n) %>%
  tibble::column_to_rownames("TF_from") %>%
  as.matrix()

# Ensure rows/cols are in the same order
interaction_matrix <- interaction_matrix[tf_order, tf_order, drop = FALSE]
stopifnot(identical(rownames(interaction_matrix), tf_order))
stopifnot(identical(colnames(interaction_matrix), tf_order))

## 2) Colors and sector labels

base_pal   <- colorRampPalette(brewer.pal(11, "Spectral"))(length(tf_order))
set.seed(42)  # reproducible shuffle
tf_colors  <- setNames(sample(base_pal), tf_order)
names(tf_colors) <- tf_order

sector_labels <- paste0(tf_summary$tf, " (", tf_summary$n_targets, "/", tf_summary$n_beaters, ")")
names(sector_labels) <- tf_summary$tf   

## Per-link color: inherit emitter color; transparent if no link
color_mat <- matrix(tf_colors[rownames(interaction_matrix)],
                    nrow = nrow(interaction_matrix), ncol = ncol(interaction_matrix),
                    dimnames = dimnames(interaction_matrix))
color_mat[interaction_matrix == 0] <- "transparent"

## Drawing helper (used for on-screen and file outputs)
circos.clear()
circos.par(start.degree = 90, gap.degree = 1, track.margin = c(0.02, 0.02))

chordDiagram(
  x               = interaction_matrix,
  order           = tf_order,
  grid.col        = tf_colors,
  col             = color_mat,
  link.arr.col    = color_mat,
  link.border     = NA,
  directional     = TRUE,
  direction.type  = "arrows",
  link.arr.length = 0.05,
  link.arr.width  = 0.05,
  transparency    = 0.25,
  annotationTrack = "grid"
)

# Safety check
stopifnot(identical(get.all.sector.index(), tf_order))

## Add outer labels (outside the ring)
circos.trackPlotRegion(
  track.index = 1,
  bg.border   = NA,
  panel.fun   = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    xlim <- get.cell.meta.data("xlim")
    ylim <- get.cell.meta.data("ylim")
    circos.text(
      x       = mean(xlim),
      y       = ylim[2] + mm_y(2.5),
      labels  = sector_labels[sec],
      facing  = "clockwise",
      niceFacing = FALSE,
      adj     = c(0, 0.5),
      cex     = 0.75
    )
  }
)
title("Circos plot of directional interactions among 50 thyroid Transcriptional Master Regulators (Outgoing vs. Incoming TFs)")


#############################################
# 7) Expresión analysis
#############################################

deg_tcga <- vroom::vroom("results/DEGs_limma_TCGA.tsv", .name_repair = janitor::make_clean_names)
deg_geo <- vroom::vroom("results/DEGs_limma_GEO.tsv", .name_repair = janitor::make_clean_names)

tmr_df <- deg_tcga %>%
  dplyr::select(gene_name, t_tcga = t, log_fc_tcga = log_fc, p_tcga = p_value) %>%
  inner_join(
    deg_geo %>% select(gene_name, t_geo = t, log_fc_geo = log_fc, p_geo = p_value),
    by = "gene_name"
  ) %>%
  filter(gene_name %in% tf_order) %>%
  mutate(
    direction_tcga = sign(log_fc_tcga),
    direction_geo  = sign(log_fc_geo),
    same_direction = direction_tcga == direction_geo
  )

table(tmr_df$same_direction)
# TRUE
#   50

tmr_meta_fisher <- tmr_df %>%
  rowwise() %>%
  mutate(p_meta_fisher = sumlog(c(p_tcga, p_geo))$p) %>%
  ungroup() %>%
  mutate(p_meta_fisher_adj = p.adjust(p_meta_fisher, "BH"))

tmr_meta <- tmr_df %>%
  mutate(
    se_tcga = abs(log_fc_tcga) / pmax(abs(t_tcga), .Machine$double.eps),
    se_geo  = abs(log_fc_geo)  / pmax(abs(t_geo),  .Machine$double.eps),
    w_tcga  = 1 / (se_tcga^2),
    w_geo   = 1 / (se_geo^2),

    lfc_meta = (w_tcga*log_fc_tcga + w_geo*log_fc_geo) / (w_tcga + w_geo),
    se_meta  = sqrt(1/(w_tcga + w_geo)),
    z_meta_iv = lfc_meta / se_meta,
    p_meta_iv = 2*pnorm(abs(z_meta_iv), lower.tail = FALSE),
    p_meta_iv_adj = p.adjust(p_meta_iv, "BH"),
    dir = if_else(lfc_meta > 0, "up", "down")
  )


meta_sig <- tmr_meta %>%
  transmute(gene_name, dir, sig = p_meta_iv_adj < 0.05)

tf_mode <- binding_hits_annot %>%
  semi_join(meta_sig, by = c("TF_from"   = "gene_name")) %>%
  semi_join(meta_sig, by = c("TF_target" = "gene_name")) %>%
  left_join(rename(meta_sig, dir_from = dir, sig_from = sig),   by = c("TF_from"   = "gene_name")) %>%
  left_join(rename(meta_sig, dir_to   = dir, sig_to   = sig),   by = c("TF_target" = "gene_name")) %>%
  mutate(
    mode = case_when(
      sig_from & sig_to & dir_from == dir_to ~ "possible_activation",
      sig_from & sig_to & dir_from != dir_to ~ "possible_repression",
      TRUE                                   ~ "uncertain"
    )
  )



# ========= 1) META direction by TF =========
meta_dir <- tmr_meta %>%
  transmute(tf = gene_name, dir, sig = p_meta_iv_adj < 0.05)

stopifnot(all(tf_order %in% meta_dir$tf))

# ========= 2) Classify links: activation/repression =========
edge_df <- binding_hits_annot %>%
  count(TF_from, TF_target, name = "n_hits") %>%
  left_join(rename(meta_dir, dir_from = dir, sig_from = sig), by = c("TF_from" = "tf")) %>%
  left_join(rename(meta_dir, dir_to   = dir, sig_to   = sig), by = c("TF_target" = "tf")) %>%
  mutate(
    mode = case_when(
      sig_from & sig_to & dir_from == dir_to ~ "possible_activation",
      sig_from & sig_to & dir_from != dir_to ~ "possible_repression",
      TRUE                                   ~ "uncertain"
    )
  )

# ========= 3) ARACNe network support =========

aracne_tcga <- read_tsv("tcga_tumor_network.txt") %>% 
      rename(tf = Regulator, target = Target, mi = MI)
aracne_geo <- read_tsv("geo_tumor_network.txt") %>% 
      rename(tf = Regulator, target = Target, mi = MI)


aracne_union <- bind_rows(
  aracne_tcga %>% mutate(ds = "tcga"),
  aracne_geo  %>% mutate(ds = "geo")
) %>%
  count(tf, target, name = "aracne_support") # 1 o 2

edge_df2 <- edge_df %>%
  left_join(aracne_union, by = c("TF_from" = "tf", "TF_target" = "target")) %>%
  mutate(
    aracne_support = replace_na(aracne_support, 0L)
  )

# Filter by minimum support
min_support <- 1
edge_plot <- edge_df2 %>%
  filter(aracne_support >= min_support)

# ========= 4) Matrix for chord and colors =========

edge_clean <- edge_plot %>%
  transmute(TF_from, TF_target, weight = n_hits * aracne_support) %>%
  filter(!is.na(TF_from), !is.na(TF_target))

# “observed” matrix (only for present pairs)
mat_obs <- xtabs(weight ~ TF_from + TF_target,
                data = edge_clean) %>%
            as.matrix()

full_mat <- matrix(0, nrow = length(tf_order), ncol = length(tf_order),
                   dimnames = list(tf_order, tf_order))
common_r <- intersect(rownames(mat_obs), tf_order)
common_c <- intersect(colnames(mat_obs), tf_order)
full_mat[common_r, common_c] <- mat_obs[common_r, common_c]

interaction_mat <- full_mat

base_pal <- colorRampPalette(brewer.pal(11, "Spectral"))(length(tf_order))
set.seed(42)
tf_colors <- setNames(sample(base_pal), tf_order)

# Link color matrix inheriting from the emisor
color_mat <- matrix(tf_colors[rownames(interaction_mat)],
                    nrow = nrow(interaction_mat), ncol = ncol(interaction_mat),
                    dimnames = dimnames(interaction_mat))
color_mat[interaction_mat == 0] <- "transparent"

grid_cols <- tf_colors

# ========= 5) Circus drawing =========
circos.clear()
circos.par(start.degree = 90, gap.degree = 1, track.margin = c(0.02, 0.02))

chordDiagram(
  x               = interaction_mat,
  order           = tf_order,
  grid.col        = tf_colors,
  col             = color_mat,
  link.arr.col    = color_mat,
  link.border     = NA,
  directional     = TRUE,
  direction.type  = "arrows",
  link.arr.length = 0.05,
  link.arr.width  = 0.05,
  transparency    = 0.25,
  annotationTrack = "grid",
  preAllocateTracks = list(
    list(track.height = uh(4, "mm")),
    list(track.height = uh(3, "mm"))
  )
)

# ========= 6) Outer ring: meta-direction (up/down) =========

dir_cols <- c(up = "#E41A1C", down = "#377EB8")  # red=up, blue=down
names(dir_cols) <- names(dir_cols)

old_pad <- circos.par("cell.padding")
circos.par(cell.padding = c(0, 0, 0, 0)) 

circos.trackPlotRegion(
  track.index = 2,
  bg.border   = NA,
  panel.fun   = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    dir  <- meta_dir$dir[match(sec, meta_dir$tf)]
    col  <- ifelse(is.na(dir), "grey85", dir_cols[dir])
    xl   <- get.cell.meta.data("xlim")
    circos.rect(xl[1], 0, xl[2], 1, col = col, border = NA)
  }
)

# restores padding for the rest
circos.par(cell.padding = old_pad)

# ========= 7) Labels outside the circle with (n_targets/n_beaters) =========
label_df <- tibble::tibble(tf = tf_order) %>%
  mutate(
    n_targets_mi  = purrr::map_int(tf, ~ sum(interaction_mat[.x, ] > 0, na.rm = TRUE)),
    n_beaters_mi  = purrr::map_int(tf, ~ sum(interaction_mat[, .x] > 0, na.rm = TRUE)),
    label         = paste0(tf, " (", n_targets_mi, "/", n_beaters_mi, ")")
  )
sector_labels <- setNames(label_df$label, label_df$tf)

circos.trackPlotRegion(
  track.index = 1,
  bg.border   = NA,
  panel.fun   = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    xl   <- get.cell.meta.data("xlim")
    circos.text(
      x       = mean(xl),
      y       = 0.5, 
      labels  = sector_labels[sec],
      facing  = "clockwise",
      niceFacing = FALSE,
      adj     = c(0, 0.5),
      cex     = 0.75
    )
  }
)


# Legends
lg_dir <- Legend(
  labels = c("Up-regulated", "Down-regulated"),
  legend_gp = gpar(fill = c("#D62728", "#1F77B4")),
  title = "Expression direction",
  grid_width = unit(3, "mm"),
  grid_height = unit(3, "mm")
)

lg_aracne <- Legend(
  labels = c("High support", "Low support"),
  legend_gp = gpar(lwd = c(2, 0.5), col = "grey20"),
  title = "ARACNe support",
  type = "lines"
)

lg_pack <- packLegend(lg_dir, lg_aracne, direction = "horizontal")

draw(
  lg_pack,
  x    = unit(0.08, "npc"), 
  y    = unit(0.08, "npc"), 
  just = c("left", "bottom")
)

################################################################################
################################################################################

# Save to PDF and PNG
pdf_file <- paste0(plotsFolder, "TF_binding_circos.pdf")
png_file <- paste0(plotsFolder, "TF_binding_circos.png")

# Dimensions in inches or pixels
pdf_wh  <- 16
png_wh  <- 12400
png_res    <- 600 
dir_cols <- c(up = "#E41A1C", down = "#377EB8")  # rojo=up, azul=down
names(dir_cols) <- names(dir_cols)



# 1) PDF save
pdf(pdf_file, width = pdf_wh+0.5, height = pdf_wh)
circos.clear()
circos.par(start.degree = 90, gap.degree = 1, track.margin = c(0.02, 0.02))
chordDiagram(
  x               = interaction_mat,
  order           = tf_order,
  grid.col        = tf_colors,
  col             = color_mat,
  link.arr.col    = color_mat,
  link.border     = NA,
  directional     = TRUE,
  direction.type  = "arrows",
  link.arr.length = 0.05,
  link.arr.width  = 0.05,
  transparency    = 0.25,
  annotationTrack = "grid",
  preAllocateTracks = list(
    list(track.height = uh(4, "mm")),
    list(track.height = uh(3, "mm"))
  ))
old_pad <- circos.par("cell.padding")
circos.par(cell.padding = c(0, 0, 0, 0)) 
circos.trackPlotRegion(
  track.index = 2,
  bg.border   = NA,
  panel.fun   = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    dir  <- meta_dir$dir[match(sec, meta_dir$tf)]
    col  <- ifelse(is.na(dir), "grey85", dir_cols[dir])
    xl   <- get.cell.meta.data("xlim")
    circos.rect(xl[1], 0, xl[2], 1, col = col, border = NA)
  })
circos.par(cell.padding = old_pad)
circos.trackPlotRegion(
  track.index = 1,
  bg.border   = NA,
  panel.fun   = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    xl   <- get.cell.meta.data("xlim")
    circos.text(
      x       = mean(xl),
      y       = 0.5,
      labels  = sector_labels[sec],
      facing  = "clockwise",
      niceFacing = FALSE,
      adj     = c(0, 0.5),
      cex     = 0.75
    )
  })
draw(
  lg_pack,
  x    = unit(0.06, "npc"), 
  y    = unit(0.06, "npc"), 
  just = c("left", "bottom"))
title("Circos plot of directional interactions among 50 thyroid Transcriptional Master Regulators (Outgoing vs. Incoming TFs)")
dev.off()

# 2) PNG save
png(png_file, width = png_wh, height = png_wh, res = png_res)
circos.clear()
circos.par(start.degree = 90, gap.degree = 1, track.margin = c(0.02, 0.02))
chordDiagram(
  x               = interaction_mat,
  order           = tf_order,
  grid.col        = tf_colors,
  col             = color_mat,
  link.arr.col    = color_mat,
  link.border     = NA,
  directional     = TRUE,
  direction.type  = "arrows",
  link.arr.length = 0.05,
  link.arr.width  = 0.05,
  transparency    = 0.25,
  annotationTrack = "grid",
  preAllocateTracks = list(
    list(track.height = uh(4, "mm")),
    list(track.height = uh(3, "mm"))
  ))
old_pad <- circos.par("cell.padding")
circos.par(cell.padding = c(0, 0, 0, 0))
circos.trackPlotRegion(
  track.index = 2,
  bg.border   = NA,
  panel.fun   = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    dir  <- meta_dir$dir[match(sec, meta_dir$tf)]
    col  <- ifelse(is.na(dir), "grey85", dir_cols[dir])
    xl   <- get.cell.meta.data("xlim")
    circos.rect(xl[1], 0, xl[2], 1, col = col, border = NA)
  })
circos.par(cell.padding = old_pad)
circos.trackPlotRegion(
  track.index = 1,
  bg.border   = NA,
  panel.fun   = function(x, y) {
    sec  <- get.cell.meta.data("sector.index")
    xl   <- get.cell.meta.data("xlim")
    circos.text(
      x       = mean(xl),
      y       = 0.5, 
      labels  = sector_labels[sec], 
      facing  = "clockwise",
      niceFacing = FALSE,
      adj     = c(0, 0.5),
      cex     = 0.75
    )
  })
draw(
  lg_pack,
  x    = unit(0.08, "npc"),  
  y    = unit(0.08, "npc"),   
  just = c("left", "bottom"))
title("Circos plot of directional interactions among 50 thyroid Transcriptional Master Regulators (Outgoing vs. Incoming TFs)")
dev.off()


