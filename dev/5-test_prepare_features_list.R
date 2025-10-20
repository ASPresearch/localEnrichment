features_example <- data.frame(
  Feature_ID = paste0("HMDB0000", 161:170),
  correlation = c(0.35, 0.42, 0.12, 0.29, -0.25, -0.32, -0.18, 0.10, 0.05, -0.02),
  p.adj = c(0.01, 0.04, 0.30, 0.10, 0.25, 0.60, 0.80, 0.05, 0.90, 0.40),
  chemical_families = rep(c("Amino Acids", "Biogenic Amines"), each = 5)
)
require(dplyr)
prepare_features_list(features_example, fdr_cutoff = 0.25, cor_cutoff = 0.1)
