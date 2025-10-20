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

# Vector de scores (p.ex. correlacions o loadings)
scores <- c(
  HMDB0000161 = 0.35, HMDB0000191 = 0.42,
  HMDB0000243 = 0.12, HMDB0000284 = 0.29,
  HMDB0001432 = -0.25, HMDB0001448 = -0.32,
  HMDB0001820 = -0.18, HMDB0000221 = 0.10,
  HMDB0000562 = 0.05, HMDB0000625 = -0.02,
  HMDB0000933 = -0.11
)

res_qsea <- qea_test(Eset, scores, nperm = 1000)
print(res_qsea)
# plot_enrichment(res_qsea, type = "dot")
