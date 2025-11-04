
# localEnrichment <img src="man/figures/logo.png" align="right" width="120" />

**Local Enrichment Analysis Tools for Omics Data**

[![R-CMD-check](https://github.com/aspresearch/localEnrichmentR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alexsanchezbio/localEnrichmentR/actions/workflows/R-CMD-check.yaml)

------------------------------------------------------------------------

## 📦 Overview

`localEnrichment` provides flexible tools for performing
**Over-Representation Analysis (ORA)** and **Quantitative Enrichment
Analysis (QEA)** on omics data, using *local, user-defined mappings*
(e.g. pathways, chemical classes, ontologies).

Unlike other enrichment packages, it is designed to work **without
relying on online databases**, allowing reproducible, standalone
analyses for metabolomics, transcriptomics, or proteomics data.

------------------------------------------------------------------------

## 🧭 Installation

You can install the development version directly from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("aspresearch/localEnrichmentR")
```

# Quickstart

``` r

library(localEnrichment)

# Example metadata and sets
meta <- list(
  mapping_name = "chemical_classes",
  feature_id_type = "HMDB_ID",
  set_source = "MEGA",
  version = "v2025.1",
  description = "Chemical classes mapping"
)

sets <- data.frame(
  set_id = c("CL001", "CL002", "CL003"),
  set_name = c("Amino Acids", "Biogenic Amines", "Lipids"),
  feature_ids = c(
    "HMDB0000161;HMDB0000191;HMDB0000243;HMDB0000284",
    "HMDB0001432;HMDB0001448;HMDB0001820",
    "HMDB0000221;HMDB0000562;HMDB0000625;HMDB0000933"
  )
)

Eset <- EnrichmentSet(sets, meta)

# Perform Over-Representation Analysis (ORA)
selected <- c("HMDB0000161", "HMDB0000284", "HMDB0001820")
res <- ora_test(Eset, selected)
plot_enrichment(res, type = "bar")
```

# Main features

EnrichmentSet() — class to handle enrichment mappings and metadata

ora_test() — over-representation analysis using Fisher’s test

qea_test() — quantitative enrichment analysis

plot_enrichment() — bar or dot plots for enriched sets

plot_qea_profile() — enrichment profiles (ES curves) for QEA results

filter_sets_by_features() — reduce large mappings to relevant sets

prepare_features_list() — filter and extract feature lists for
enrichment

# Authors

Alex Sanchez-Pla (<asanchez@ub.edu>)

# License

MIT + file LICENSE
