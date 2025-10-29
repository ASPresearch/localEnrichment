#' #' Build an EnrichmentSet object from a mapping table
#' #'
#' #' @description
#' #' Create an \code{EnrichmentSet} from a mapping between metabolite IDs
#' #' and category annotations (e.g., KEGG pathways, chemical classes).
#' #'
#' #' @param data A data frame containing mapping information.
#' #' @param id_col Character. Column name with metabolite IDs (e.g. "HMDB").
#' #' @param category_col Character. Column name with category names (e.g. "KEGG").
#' #' @param set_name Character. Optional descriptive name for the enrichment set.
#' #' @param source Character. Optional database or mapping source.
#' #' @param species Character. Species name (default = "Homo sapiens").
#' #'
#' #' @return An object of class \code{EnrichmentSet}.
#' #' @export
#' buildEnrichmentSet <- function(data,
#'                                id_col,
#'                                category_col,
#'                                set_name = NULL,
#'                                source = NULL,
#'                                species = "Homo sapiens") {
#'   if (!requireNamespace("dplyr", quietly = TRUE))
#'     stop("Package 'dplyr' is required.")
#'   library(dplyr)
#'
#'   stopifnot(id_col %in% names(data), category_col %in% names(data))
#'
#'   df <- data %>%
#'     dplyr::select(all_of(c(id_col, category_col))) %>%
#'     dplyr::filter(!is.na(.data[[category_col]]), .data[[category_col]] != "") %>%
#'     dplyr::distinct()
#'
#'   # --- Collapse IDs per category correctly ---
#'   set_df <- df %>%
#'     dplyr::group_by(.data[[category_col]]) %>%
#'     dplyr::summarise(feature_ids = paste(unique(.data[[id_col]]), collapse = ";"),
#'                      .groups = "drop") %>%
#'     dplyr::mutate(
#'       set_id = .data[[category_col]],
#'       set_name = .data[[category_col]]
#'     ) %>%
#'     dplyr::select(set_id, set_name, feature_ids)
#'
#'   # --- Metadata list ---
#'   meta <- list(
#'     mapping_name = ifelse(is.null(set_name), category_col, set_name),
#'     feature_id_type = id_col,
#'     feature_species = species,
#'     set_source = ifelse(is.null(source), "user_defined", source),
#'     version = as.character(Sys.Date()),
#'     description = paste("Mapping between", id_col, "and", category_col)
#'   )
#'
#'   # --- Create the EnrichmentSet ---
#'   EnrichmentSet(data = set_df, metadata = meta)
#' }


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

