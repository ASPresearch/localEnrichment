meta <- list(
  mapping_name = "smpdb_pathways",
  feature_id_type = "HMDB_ID",
  feature_species = "Homo sapiens",
  set_source = "SMPDB",
  version = "SMPDB v2.0 (2024-05)",
  description = "Metabolic pathways curated from SMPDB"
)

sets <- data.frame(
  set_id = c("PW001", "PW002"),
  set_name = c("Glycolysis", "TCA Cycle"),
  feature_ids = c("HMDB0000060;HMDB0000122", "HMDB0000209;HMDB0000210"),
  description = c("Energy metabolism", "Mitochondrial cycle")
)

Eset <- EnrichmentSet(sets, meta)
show(Eset)
summary(Eset)

# Convert to list for ORA/QEA functions
sets_list <- as(Eset, "list")
str(sets_list)

# Coerce to data.frame
df <- as.data.frame(Eset)

head(df)

# Convert to list again
lst <- as(EnrichmentSet, "list")
length(lst)
#> number of sets

# Coerce into MetaboliteSet
# There are three possible outcomes for pathway names

as.MetaboliteSetDataFrame(Eset)
as.MetaboliteSetDataFrame(Eset, id_type ="id")
as.MetaboliteSetDataFrame(Eset, id_type ="both")
