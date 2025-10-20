meta <- list(
  mapping_name = "chemical_classes",
  feature_id_type = "HMDB_ID",
  feature_species = "Homo sapiens",
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
selected <- c("HMDB0000161", "HMDB0000284", "HMDB0001820")

res <- ora_test(Eset, selected)
plot_enrichment(res, type = "bar")
plot_enrichment(res, type = "dot")
