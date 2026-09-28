# PANORAMIC
**P**ooled **AN**alysis **O**f Va**R**iance-**A**ware **M**odeling and **I**nference of **C**olocalization (**PANORAMIC**)

[![DOI](https://zenodo.org/badge/1086842429.svg)](https://doi.org/10.5281/zenodo.19927197)

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/plevritis-lab/panoramic/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/plevritis-lab/panoramic/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

PANORAMIC provides variance-aware multi-sample spatial colocalization
analysis for spatial omics studies. It estimates sample-level
colocalization effects, quantifies within-sample uncertainty with spatial
bootstrap procedures, and propagates that uncertainty through multilevel
random-effects meta-analysis to compare biological groups while accounting
for heterogeneity across samples and patients.

Its default statistic, `local_comp_enrichment`, measures edge-corrected local
competition enrichment for ordered cell-type pairs. The optional
`local_comp_global_enrichment` statistic instead uses a random-label null
conditional on observed cell locations; this can help study label mixing in
density-heterogeneous tissue, but it does not establish direct cell-cell
interaction.

The method is described in Chang et al., *Bioinformatics* 42(8): btag546
(published July 31, 2026):
https://doi.org/10.1093/bioinformatics/btag546

![alt text](https://github.com/plevritis-lab/panoramic/blob/main/panoramic_graphic2.png)

## Installation
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install("panoramic")
library(panoramic)
```

To install the development version from GitHub (optional):

```r
# install.packages("remotes")
remotes::install_github("plevritis-lab/panoramic")
```

## Documentation

The package vignette walks through sample preparation, spatial-statistic
estimation, multilevel pooling, differential testing, and visualization:

```r
browseVignettes("panoramic")
```

## Citation

If you use PANORAMIC in published work, please cite the Bioinformatics
article and, when appropriate, the Zenodo software archive:

- https://doi.org/10.1093/bioinformatics/btag546
- https://doi.org/10.5281/zenodo.19927197
