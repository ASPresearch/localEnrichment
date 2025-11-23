# --- 1. Definició de l’EnrichmentSet ---------------------------------------

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


# --- 2. Scores de les característiques -------------------------------------

scores <- c(
  HMDB0000161 = 0.35, HMDB0000191 = 0.42,
  HMDB0000243 = 0.12, HMDB0000284 = 0.29,
  HMDB0001432 = -0.25, HMDB0001448 = -0.32,
  HMDB0001820 = -0.18, HMDB0000221 = 0.10,
  HMDB0000562 = 0.05, HMDB0000625 = -0.02,
  HMDB0000933 = -0.11
)

# --- 3. Funció auxiliar per generar el perfil ------------------------------

qea_profile <- function(eset, scores, set_name) {
  if (!inherits(eset, "EnrichmentSet"))
    stop("Input must be an EnrichmentSet.")

  sets_list <- as(eset, "list")   # clau: set_id

  # 👉 permetre noms humans
  if (!(set_name %in% names(sets_list))) {
    lookup <- eset@data
    if (set_name %in% lookup$set_name) {
      set_name <- lookup$set_id[match(set_name, lookup$set_name)]
    } else {
      stop("Specified set_name not found in EnrichmentSet.")
    }
  }

  hits <- sets_list[[set_name]]

  df <- data.frame(
    feature = names(sort(scores, decreasing = TRUE)),
    score = sort(scores, decreasing = TRUE)
  )
  df$hit <- df$feature %in% hits

  Nh <- sum(df$hit)
  N <- nrow(df)
  if (Nh == 0) stop("No hits for this set in scores.")

  running_score <- cumsum(df$hit / Nh - (!df$hit) / (N - Nh))

  data.frame(
    position = seq_len(N),
    running_score = running_score,
    hit = df$hit
  )
}


# --- 5 SImplifiquem l'Eset

Eset_small <- filter_sets_by_features(Eset, names(scores))

# 2. QEA complet amb perfils
res_qea <- qea_test(Eset_small, scores,
                    min_set_size = 2,
                    return_profiles = TRUE,
                    nperm = 1000)

# 3. Plot per un set concret (per ID)
plot_qea_profile(res_qea, sets = "CL001")

# O per nom humà:
plot_qea_profile(res_qea, sets = "Amino Acids")

# 4. Tots els significatius
plot_qea_profile(res_qea, show = "significant", p_cutoff = 0.25)
