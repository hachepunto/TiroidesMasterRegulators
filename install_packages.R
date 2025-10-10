# install_packages.R
# Utility for install packages used in this repository
# José Daniel López Mendiola · october 2025 danolopez0614@gmail.com

# Usage
# RScript install_packages.R

if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# List of Bioconductor packages (including dependencies)
bioc_packages <- c( "GO.db", "GOSemSim", "DOSE", "enrichplot",  # clusterProfiler dependencies
                    "biomaRt", "Biostrings", "BSgenome.Hsapiens.UCSC.hg38", 
                    "clusterProfiler", "ComplexHeatmap", "DESeq2", "edgeR", 
                    "GenomicRanges", "GEOquery", "JASPAR2022", "org.Hs.eg.db", 
                    "TCGAbiolinks", "TFBSTools", "viper")

# List of CRAN packages
cran_packages <- c( "affy", "circlize", "ComplexUpset", "cowplot", 
                    "dplyr", "EnhancedVolcano", "forcats", "ggplot2", "ggtext", 
                    "gplots", "gprofiler2", "grid", "igraph", "janitor", "limma", 
                    "magrittr", "metap", "patchwork", "pheatmap", "purrr", 
                    "RColorBrewer", "readr", "scales", "stringr", "tibble", "tidyr",
                    "tidyverse", "vroom")

# Install CRAN packages
cat("Installing CRAN packages...\n")
for (pkg in cran_packages) {
  if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
    cat("Installing:", pkg, "\n")
    install.packages(pkg, dependencies = TRUE)
  }
}

# Install Bioconductor packages (with all dependencies)
cat("Installing Bioconductor packages...\n")
for (pkg in bioc_packages) {
  if (!require(pkg, quietly = TRUE, character.only = TRUE)) {
    cat("Installing:", pkg, "\n")
    BiocManager::install(pkg, dependencies = TRUE, ask = FALSE, update = TRUE)
  }
}

# Specific installation of clusterProfiler with all its dependencies
if (!require("clusterProfiler", quietly = TRUE)) {
  cat("Installing clusterProfiler and its dependencies...\n")
  BiocManager::install("clusterProfiler", dependencies = TRUE, ask = FALSE)
}

# Verify installation of all packages
cat("\nVerifying package installation...\n")

for (pkg in c(bioc_packages, cran_packages)) {
  if (require(pkg, quietly = TRUE, character.only = TRUE)) {
    cat("✓", pkg, "loaded correctly\n")
  } else {
    cat("✗ Error loading:", pkg, "\n")
  }
}

# Try to load clusterProfiler specifically
cat("\nTesting clusterProfiler...\n")
if (!require("clusterProfiler", quietly = TRUE)) {
  # Install dependencies manually if necessary
  BiocManager::install(c("GO.db", "GOSemSim", "DOSE", "enrichplot"), ask = FALSE)
  BiocManager::install("clusterProfiler", ask = FALSE)
}
