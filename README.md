# Two Cohorts, One Network: Consensus Master Regulators Orchestrating Papillary Thyroid Carcinoma

**Author:** Hugo Tovar  
**Institution:** National Institute of Genomic Medicine (INMEGEN), Mexico  
**Contact:** hatovar@inmegen.gob.mx

## Overview

This repository contains the scripts, helper functions, and reproducible pipelines used in the study:

*Two Cohorts, One Network: Consensus Master Regulators Orchestrating Papillary Thyroid Carcinoma*

**Published article**  
Tapia-Carrillo, D.; Zambada-Moreno, O.; Hernández-Lemus, E.; Tovar, H. (2025). *Two Cohorts, One Network: Consensus Master Regulators Orchestrating Papillary Thyroid Carcinoma*. **International Journal of Molecular Sciences**, 26(22): 11231.  
DOI: http://doi.org/10.3390/ijms262211231

This analysis integrates transcriptomic data from two independent cohorts — TCGA-THCA (RNA-seq) and GEO GSE33630 (microarray) — to identify Transcriptional Master Regulators (TMRs) through a unified network-based framework.

## Repository contents

- Reproducible analysis scripts
- Helper functions
- Pipeline components used in the manuscript
- Submodules required to reproduce the workflow

## How to clone this repository

```bash
git clone https://github.com/hachepunto/TiroidesMasterRegulators.git
cd TiroidesMasterRegulators
git submodule update --init --recursive
Rscript install_packages.R  # optional
```

## Citation

If you use this repository or its code in your work, please cite the published article above.

> “Networks reveal what lists conceal.”  
> — Systems Biology proverb
