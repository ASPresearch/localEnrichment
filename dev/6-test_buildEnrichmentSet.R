meta_df <- data.frame(
  HMDB = c("HMDB0000060", "HMDB0000122", "HMDB0000209", "HMDB0000210"),
  Pathway = c("Glycolysis", "Glycolysis", "TCA Cycle", "TCA Cycle")
)

Eset <- buildEnrichmentSet(
  data = meta_df,
  id_col = "HMDB",
  category_col = "Pathway",
  set_name = "SMPDB_pathways",
  source = "SMPDB",
  species = "Homo sapiens",
  version = "2024-05",
  description = "Metabolic pathways curated from SMPDB"
)

show(Eset)
summary(Eset)
as.data.frame(Eset)


# --- 0. Get mapping data

load(file="./inst/extdata/myProject_map.Rda")
head(myProject_map[,1:6])

# 1️⃣ KEGG-based enrichment
KEGGset <- buildEnrichmentSet(
  data = myProject_map,
  id_col = "HMDB",
  category_col = "KEGG",
  set_name = "KEGG_pathways",
  source = "MetaboAnalyst mapping"
)

metabsChemicalClasses<- openxlsx::read.xlsx("./inst/extdata/MEGA-r_name-HMDB-ChemicalClasses.xlsx")
# 2️⃣ Chemical classes enrichment
ChemClassSet <- buildEnrichmentSet(
  data = metabsChemicalClasses,
  id_col = "HMDB",
  category_col = "ChemicalClass",
  set_name = "ChemicalClasses",
  source = "MEGA annotation"
)

# 3️⃣ Visualitzar resum
summary(KEGGset)
summary(ChemClassSet)

## SMPDB pathways. Generated as a nested list. Casted into an EnrichmentSet

load("./inst/extdata/smpdb_pathway.rda")
SMPDBset <- as.EnrichmentSet.list(smpdb_pathway)
summary(SMPDBset)

## Filter it to keep only the sets you are interested in

hmdb_study <- unique(myProject_map$HMDB)

# Filter the SMPDB enrichment set
SMPDBset_filtered <- filterEnrichmentSet(SMPDBset, hmdb_study)
summary(SMPDBset_filtered)


SMPDBset_filtered<- as.data.frame(SMPDBset_filtered)
head(SMPDBset_filtered)

