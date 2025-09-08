# scripts/helpers.R
# Utility functions for thyroid project analyses
# Hugo Tovar · July 2025

# Helper: Get the mode (most frequent value)
getmode <- function(v) {
 uniqv <- unique(v)
 uniqv[which.max(tabulate(match(v, uniqv)))]
}

# Helper: Download data if not already present
download_if_needed <- function(query, directory = "GDCdata") {
 if (!dir.exists(directory)) {
  message("➤ Downloading data to ", directory)
  GDCdownload(query)
 } else {
  message("✔ Data directory '", directory, "' already exists. Skipping download.")
 }
}

# Helper: Extract top transcriptional regulators from msviper object
viper_mrsTopTable <- function(mrs, n = sum(mrs$es$p.value < p_threshold), p_threshold = 0.05){
 mrs_tops <- tibble::tibble(
  TF = names(mrs$es$nes[mrs$es$p.value < p_threshold]),
  nes = mrs$es$nes[mrs$es$p.value < p_threshold],
  p.value = mrs$es$p.value[mrs$es$p.value < p_threshold]
 )
 mrs_tops <- mrs_tops[order(abs(mrs_tops$nes), decreasing = TRUE),]
 return(mrs_tops[1:n,])
}

# Helper: Get regulon as tidy table
getregulon <- function(regulon, 
            tf_list = names(regulon), 
            likelihood_threshold = 0, 
            export = FALSE, 
            filename = "regulon.tsv") {
 
 tsv_list <- lapply(tf_list, function(name) {
  if (is.null(regulon[[name]])) return(NULL)

  tf <- rep(name, length(regulon[[name]]$tfmode))
  likelihoods <- regulon[[name]]$likelihood
  targets <- names(regulon[[name]]$tfmode)
  correlations <- regulon[[name]]$tfmode
  
  keep <- likelihoods >= likelihood_threshold
  if (sum(keep) == 0) return(NULL)
  
  tibble::tibble(
   tf = tf[keep],
   target = targets[keep],
   correlation = correlations[keep],
   likelihood = likelihoods[keep]
  ) |> 
   dplyr::arrange(tf, dplyr::desc(likelihood))
 })

 tsv <- dplyr::bind_rows(tsv_list)
 
 if (export) {
  readr::write_tsv(tsv, filename)
  message("TSV exported to: ", filename)
 }

 return(tsv)
}


run_gsea_go <- function(df, organism_db = org.Hs.eg.db) {
 # Cargar librerías necesarias
 suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(dplyr)
 })

 # Validación de columnas
 if (!all(c("gene_name", "logFC") %in% colnames(df))) {
  stop("The data frame most have 'gene_name' and 'logFC' columns")
 }

 # Filtrar posibles NA en logFC
 df <- df %>% filter(!is.na(logFC))

 message("→ Maping gene symbols to Entrez IDs...")

 # Mapear símbolos a Entrez IDs
 gene_mapping <- bitr(df$gene_name,
            fromType = "SYMBOL",
            toType = "ENTREZID",
            OrgDb = organism_db)

 if (nrow(gene_mapping) == 0) {
  stop("None symbol could be mapped to EntrezID")
 }

 # Unir logFC con los IDs mapeados
 df_mapped <- df %>%
  inner_join(gene_mapping, by = c("gene_name" = "SYMBOL")) %>%
  group_by(ENTREZID) %>%
  summarize(logFC = mean(logFC, na.rm = TRUE), .groups = "drop")

 # Crear el vector de genes ordenado
 gene_list <- df_mapped$logFC
 names(gene_list) <- df_mapped$ENTREZID
 gene_list <- sort(gene_list, decreasing = TRUE)

 message("→ Running GSEA with GO ontology (BP and MF)...")

 # Ejecutar GSEA con ontología ALL
 gsea_all <- gseGO(geneList = gene_list,
          OrgDb = organism_db,
          ont = "ALL",
          keyType = "ENTREZID",
          minGSSize = 10,
          maxGSSize = 500,
          pvalueCutoff = 1,
          verbose = FALSE)

 # Filtrar solo BP y MF
 gsea_filtered <- gsea_all
 gsea_filtered@result <- gsea_all@result %>%
  filter(ONTOLOGY %in% c("MF", "BP"))

 message("✓ GSEA complited")

 return(list(
  gsea_result = gsea_filtered,
  gene_list = gene_list,
  gene_mapping = df_mapped
 ))
}


plot_ridge_panels <- function(gsea_result, top_n = 10, title = "GO ridgeplot") {
 library(dplyr)
 library(clusterProfiler)
 library(patchwork)

 res <- gsea_result@result

 # Top términos BP y MF por separado
 top_bp_ids <- res %>%
  filter(ONTOLOGY == "BP") %>%
  arrange(p.adjust) %>%
  slice_head(n = top_n) %>%
  pull(ID)

 top_mf_ids <- res %>%
  filter(ONTOLOGY == "MF") %>%
  arrange(p.adjust) %>%
  slice_head(n = top_n) %>%
  pull(ID)

 # Subsets de resultados
 gsea_bp <- gsea_result
 gsea_bp@result <- res %>% filter(ID %in% top_bp_ids)

 gsea_mf <- gsea_result
 gsea_mf@result <- res %>% filter(ID %in% top_mf_ids)

 # Crear gráficas individuales
 p1 <- ridgeplot(gsea_bp, showCategory = top_n) +
  ggtitle("GO: Biological Process") +
  theme_minimal()

 p2 <- ridgeplot(gsea_mf, showCategory = top_n) +
  ggtitle("GO: Molecular Function") +
  theme_minimal()

 # Combinar con título global
 (p1 / p2) + 
  patchwork::plot_annotation(title = title)
}


ora_hallmarks <- function(regulon_df, hallmark_gmt, universe, mode = c("by_tf", "all_targets")) {
 library(clusterProfiler)
 library(dplyr)

 mode <- match.arg(mode)

 if (!all(c("tf", "target") %in% colnames(regulon_df))) {
  stop("The regulon_df object must contain 'tf' and 'target' columns")
 }

 # Verifica universo en mayúsculas para símbolos (por si acaso)
 universe <- toupper(universe)
 regulon_df <- regulon_df %>%
  mutate(target = toupper(target), tf = toupper(tf))

 # función auxiliar para cada conjunto de genes
 run_fora <- function(gene_set, set_name) {
  gene_set <- intersect(gene_set, universe)
  if (length(gene_set) < 5) {
   message("Skipping ", set_name, ": only ", length(gene_set), " genes.")
   return(NULL)
  }
  res <- enricher(gene = gene_set,
          TERM2GENE = hallmark_gmt,
          universe = universe,
          pAdjustMethod = "BH",
          minGSSize = 5,
          maxGSSize = 500)

  if (is.null(res)) return(NULL)

  res_df <- as.data.frame(res@result)
  res_df$set <- set_name
  return(res_df)
 }

 if (mode == "by_tf") {
  tf_list <- unique(regulon_df$tf)
  results <- lapply(tf_list, function(tf) {
   targets_tf <- regulon_df %>%
    clusterProfiler::filter(tf == !!tf) %>%
    pull(target)
   run_fora(targets_tf, tf)
  })
 } else {
  all_targets <- unique(regulon_df$target)
  results <- list(run_fora(all_targets, "ALL_TARGETS"))
 }

 # 🧹 Limpia los NULL antes de unir
 results <- results[!sapply(results, is.null)]

 final <- bind_rows(results)
 return(final)
}

# --- ORA de GO (ALL) ---
ora_go <- function(regulon_df, universe, mode = c("by_tf", "all_targets"), organism_db = org.Hs.eg.db) {
 library(clusterProfiler)
 library(org.Hs.eg.db)
 library(dplyr)

 mode <- match.arg(mode)
 regulon_df <- regulon_df %>% mutate(target = toupper(target), tf = toupper(tf))
 universe <- toupper(universe)

 run_go <- function(gene_set, set_name) {
  gene_set <- intersect(gene_set, universe)
  if (length(gene_set) < 5) return(NULL)

  ego <- tryCatch({
   enrichGO(gene = gene_set,
        universe = universe,
        OrgDb = organism_db,
        keyType = "SYMBOL",
        ont = "ALL",
        pAdjustMethod = "BH",
        minGSSize = 5,
        maxGSSize = 500,
        readable = TRUE)
  }, error = function(e) NULL)

  if (!is.null(ego)) {
   as.data.frame(ego@result) %>% mutate(set = set_name)
  } else NULL
 }

 if (mode == "by_tf") {
  tf_list <- unique(regulon_df$tf)
  results <- lapply(tf_list, function(tf) {
   targets_tf <- regulon_df %>% filter(tf == !!tf) %>% pull(target)
   run_go(targets_tf, tf)
  })
 } else {
  all_targets <- unique(regulon_df$target)
  results <- list(run_go(all_targets, "ALL_TARGETS"))
 }

 bind_rows(results)
}

# --- ORA de KEGG ---
ora_kegg <- function(regulon_df, universe, mode = c("by_tf", "all_targets")) {
 library(clusterProfiler)
 library(dplyr)

 mode <- match.arg(mode)
 regulon_df <- regulon_df %>% mutate(target = toupper(target), tf = toupper(tf))
 universe <- toupper(universe)
 mapping <- bitr(universe, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
 symbol2entrez <- setNames(mapping$ENTREZID, mapping$SYMBOL)
 entrez_universe <- mapping$ENTREZID

 run_kegg <- function(gene_set, set_name) {
  gene_set <- intersect(gene_set, names(symbol2entrez))
  gene_ids <- symbol2entrez[gene_set]
  gene_ids <- gene_ids[!is.na(gene_ids)]

  if (length(gene_ids) < 5) return(NULL)

  ekegg <- tryCatch({
   enrichKEGG(gene = gene_ids,
         universe = entrez_universe,
         organism = "hsa",
         pAdjustMethod = "BH",
         minGSSize = 5,
         maxGSSize = 500)
  }, error = function(e) NULL)

  if (!is.null(ekegg)) {
   as.data.frame(ekegg@result) %>% mutate(set = set_name)
  } else NULL
 }

 if (mode == "by_tf") {
  tf_list <- unique(regulon_df$tf)
  results <- lapply(tf_list, function(tf) {
   targets_tf <- regulon_df %>% filter(tf == !!tf) %>% pull(target)
   run_kegg(targets_tf, tf)
  })
 } else {
  all_targets <- unique(regulon_df$target)
  results <- list(run_kegg(all_targets, "ALL_TARGETS"))
 }

 bind_rows(results)
}

# Function to find consensus TMR targets shared across both regulons
 find_consensus_targets <- function(tf, tmrs_vector, regulon1, regulon2) {
 if (!(tf %in% names(regulon1)) || !(tf %in% names(regulon2))) {
  return(NULL)
 }

 targets1 <- names(regulon1[[tf]]$tfmode)
 targets2 <- names(regulon2[[tf]]$tfmode)

 consensus_targets <- intersect(targets1, targets2)

 # Filter only targets that are also TMRs
 intersect(consensus_targets, tmrs_vector)
}


# Function to build a regulon for a TF based on shared targets and geometric mean of likelihood
build_meta_regulon <- function(tf, reg1, reg2) {
 if (!(tf %in% names(reg1)) || !(tf %in% names(reg2))) {
  return(NULL)
 }
 
 tf_reg1 <- reg1[[tf]]
 tf_reg2 <- reg2[[tf]]

 # Commun targets
 targets_common <- intersect(names(tf_reg1$tfmode), names(tf_reg2$tfmode))
 if (length(targets_common) == 0) return(NULL)

 idx1 <- match(targets_common, names(tf_reg1$tfmode))
 idx2 <- match(targets_common, names(tf_reg2$tfmode))
 
 # Geometric mean of likelihood
 lik1 <- tf_reg1$likelihood[idx1]
 lik2 <- tf_reg2$likelihood[idx2]
 geom_mean_lik <- sqrt(lik1 * lik2)
 
 # Use the simple average of tfmode as an approximation
 tfmode1 <- tf_reg1$tfmode[idx1]
 tfmode2 <- tf_reg2$tfmode[idx2]
 mean_tfmode <- (tfmode1 + tfmode2) / 2
 
 list(
  tfmode = setNames(mean_tfmode, targets_common),
  likelihood = geom_mean_lik
 )
}

regulon2sif_tbl <- function(regulon){
 sif <- tibble::tibble(source = character(), target = character(), interaction = numeric()) # Crear un tibble vacío
 for (i in names(regulon)) {
  source = rep(i, length(regulon[[i]]$tfmod))
  psif = tibble::tibble(source, target = names(regulon[[i]]$tfmod), interaction = regulon[[i]]$tfmod)
  sif <- dplyr::bind_rows(sif, psif)
 }
 return(sif)
}


#' Read a PFM file from CIS-BP and convert it to a PFMatrix object
#'
#' This function attempts to parse a Position Frequency Matrix (PFM) file from
#' the CIS-BP bulk download. It supports both the "standard" format with
#' A/C/G/T headers and fallback parsing when headers are missing or malformed.
#' The function automatically distinguishes between frequency matrices (values <= 1)
#' and count matrices (values > 1), and scales frequencies to integer counts.
#'
#' @param motif_id  Character. Motif identifier (e.g. "M05089_3.00").
#' @param pwms_dir  Character. Path to the directory containing CIS-BP PFM files.
#' @param scale     Integer. Scaling factor for frequency matrices (default = 100).
#' @param verbose   Logical. If TRUE, prints messages when files cannot be read.
#'
#' @return A TFBSTools::PFMatrix object, or NULL if parsing fails.
#'
read_cisbp_pfm <- function(motif_id, pwms_dir, scale = 100, verbose = TRUE) {
  # Build file path (try both "motif_id.txt" and "motif_id")
  f <- file.path(pwms_dir, paste0(motif_id, ".txt"))
  if (!file.exists(f)) {
    f <- file.path(pwms_dir, motif_id)
    if (!file.exists(f)) {
      if (verbose) message("File not found: ", motif_id)
      return(NULL)
    }
  }
  
  mat <- NULL
  
  # --- Attempt to read with proper headers (A/C/G/T) ---
  suppressWarnings({
    df <- try(readr::read_tsv(f, col_types = readr::cols(.default="d"), progress = FALSE),
              silent = TRUE)
  })
  
  if (!inherits(df, "try-error") && all(c("A","C","G","T") %in% names(df))) {
    mat <- t(as.matrix(df[, c("A","C","G","T")]))
    rownames(mat) <- c("A","C","G","T")
  }
  
  # --- Fallback: parse as raw numeric text without headers ---
  if (is.null(mat)) {
    raw <- readLines(f, warn = FALSE)
    raw <- raw[!grepl("^\\s*(#|;|$)", raw)]   # remove comments and empty lines
    if (!length(raw)) return(NULL)
    if (grepl("[A-Za-z]", raw[1])) raw <- raw[-1]  # drop header line if it has letters
    if (!length(raw)) return(NULL)
    
    # Tokenize rows into numeric values
    tok <- strsplit(raw, "\\s+")
    maxlen <- max(lengths(tok))
    tok <- lapply(tok, function(x){ length(x) <- maxlen; x })
    dfx <- as.data.frame(do.call(rbind, tok), stringsAsFactors = FALSE)
    
    # Detect numeric columns
    is_num <- function(x) grepl("^[-+]?\\d*\\.?\\d+(e[-+]?\\d+)?$", x, ignore.case = TRUE)
    num_frac <- vapply(dfx, function(col) mean(is_num(col), na.rm = TRUE), numeric(1))
    ord <- order(num_frac, decreasing = TRUE)
    idx4 <- if (length(ord) >= 5) tail(ord, 4L) else head(ord, 4L)
    
    if (length(idx4) < 4 || any(num_frac[idx4] < 0.8)) return(NULL)
    
    mat_raw <- apply(dfx[, idx4, drop = FALSE], 2, as.numeric)
    if (ncol(mat_raw) == 4) mat_raw <- t(mat_raw)
    if (nrow(mat_raw) != 4 || ncol(mat_raw) < 1) return(NULL)
    
    rownames(mat_raw) <- c("A","C","G","T")
    mat <- mat_raw
  }
  
  # --- Final checks ---
  if (is.null(mat) || nrow(mat) != 4 || ncol(mat) < 1 || !all(is.finite(mat))) return(NULL)
  maxv <- suppressWarnings(max(mat))
  if (!is.finite(maxv)) return(NULL)
  
  # --- Frequency vs counts ---
  if (maxv <= 1.000001) {
    # Normalize columns to sum 1 if slightly off
    if (any(abs(colSums(mat) - 1) > 1e-6)) {
      mat <- apply(mat, 2, function(col) col / sum(col))
    }
    mc <- round(mat * scale)
    storage.mode(mc) <- "integer"
    dtype <- "freq->counts"
  } else {
    mc <- round(mat)
    storage.mode(mc) <- "integer"
    dtype <- "counts"
  }
  
  # --- Build PFMatrix object ---
  pfm <- TFBSTools::PFMatrix(ID = motif_id, name = motif_id,
                             strand = "+", profileMatrix = mc)
  pfm@tags$source      <- "CIS-BP"
  pfm@tags$file        <- basename(f)
  pfm@tags$data_type   <- dtype
  pfm@tags$count_scale <- if (dtype == "freq->counts") scale else NA_integer_
  pfm
}

# Minimal version: return raw numeric matrix instead of PFMatrix
read_cisbp_pfm_debug <- function(motif_id, pwms_dir, scale = 100, verbose = TRUE) {
  f <- file.path(pwms_dir, paste0(motif_id, ".txt"))
  if (!file.exists(f)) f <- file.path(pwms_dir, motif_id)
  if (!file.exists(f)) {
    if (verbose) message("File not found: ", motif_id)
    return(NULL)
  }
  
  df <- try(readr::read_tsv(f, col_types = readr::cols(.default="d"), progress = FALSE), silent = TRUE)
  
  if (!inherits(df, "try-error") && all(c("A","C","G","T") %in% names(df))) {
    mat <- t(as.matrix(df[, c("A","C","G","T")]))
    rownames(mat) <- c("A","C","G","T")
  } else {
    raw <- readLines(f, warn = FALSE)
    raw <- raw[!grepl("^\\s*(#|;|$)", raw)]
    if (length(raw) == 0) return(NULL)
    if (grepl("[A-Za-z]", raw[1])) raw <- raw[-1]
    if (length(raw) == 0) return(NULL)
    tok <- strsplit(raw, "\\s+")
    maxlen <- max(lengths(tok))
    tok <- lapply(tok, function(x){ length(x) <- maxlen; x })
    dfx <- as.data.frame(do.call(rbind, tok), stringsAsFactors = FALSE)
    mat <- apply(dfx, 2, as.numeric)
    if (ncol(mat) == 4) mat <- t(mat)
    rownames(mat) <- c("A","C","G","T")
  }
  
  if (is.null(mat) || !all(is.finite(mat))) return(NULL)
  mat
}
