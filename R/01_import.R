# 01 - Import des données
#
# Étape 1 : télécharge les fichiers sources dans data/1_brut/ (une seule fois,
#           les fichiers déjà présents sont réutilisés).
# Étape 2 : les lit, les met au même format et écrit 4 fichiers lisibles dans
#           data/2_intermediaire/. Voir data/LISEZMOI.md pour le détail.
#
# À lancer depuis le dossier projet_velo/.

library(data.table)

options(timeout = 1200)
brut  <- "data/1_brut"
inter <- "data/2_intermediaire"
for (d in file.path(brut, c("velo", "meteo", "calendrier"))) dir.create(d, recursive = TRUE, showWarnings = FALSE)
dir.create(inter, recursive = TRUE, showWarnings = FALSE)

telecharger <- function(url, fichier) {
  if (!file.exists(fichier)) download.file(url, fichier, mode = "wb", quiet = TRUE)
  fichier
}

# ===========================================================================
# 1. Comptages vélo 2018-2025 : archives annuelles de la Ville de Paris
# ===========================================================================
# Le format change selon les années (séparateur, guillemets, format de date),
# mais l'ordre des colonnes est stable :
#   1 identifiant du compteur, 2 nom du compteur, 3 identifiant du site,
#   5 comptage horaire, 6 date et heure locale.
# L'archive 2025 fait près de 2 Go décompressée : on ne lit que ces 5 colonnes.

url_archives <- paste0(
  "https://opendata.paris.fr/api/explore/v2.1/catalog/datasets/",
  "comptage-velo-historique-donnees-compteurs/attachments/"
)
archives <- c(
  "2018" = "2018_comptage_velo_donnees_compteurs_csv_zip",
  "2019" = "2019_comptage_velo_donnees_compteurs_csv_zip",
  "2020" = "2020_comptage_velo_donnees_compteurs_csv_zip",
  "2021" = "2021_comptage_velo_donnees_compteurs_csv_zip",
  "2022" = "2022_comptage_velo_donnees_compteurs_zip",
  "2023" = "2023_comptage_velo_donnees_compteurs_zip",
  "2024" = "2024_comptage_velo_donnees_compteurs_zip",
  "2025" = "2025_comptage_velo_donnees_compteurs_zip"
)

lire_archive <- function(annee) {
  zip <- telecharger(paste0(url_archives, archives[[annee]]),
                     file.path(brut, "velo", paste0("velo_archive_", annee, ".zip")))
  csv <- grep("^[^_].*\\.csv$", unzip(zip, list = TRUE)$Name, value = TRUE)[1]
  unzip(zip, files = csv, exdir = tempdir(), overwrite = TRUE)
  chemin <- file.path(tempdir(), csv)
  on.exit(unlink(chemin))

  h <- fread(chemin, select = c(1, 2, 3, 5, 6), encoding = "UTF-8",
             colClasses = list(character = c(1, 2, 3, 5, 6)),
             col.names = c("id_compteur", "nom_compteur", "id_site", "comptage", "date_heure"),
             showProgress = FALSE)
  # Formats de comptage rencontrés : "12", "12.0" (2020, 2021, 2024) et, en 2025,
  # un séparateur de milliers en espace insécable ("1 022"). On garde seulement
  # les chiffres et le point décimal, octet par octet, puis on convertit.
  # Certains noms (2025) sont encodés en Latin-1 au lieu d'UTF-8 : on les convertit.
  mal_encodes <- !validUTF8(h$nom_compteur)
  if (any(mal_encodes)) {
    h[mal_encodes, nom_compteur := iconv(nom_compteur, "latin1", "UTF-8")]
  }
  corriges <- sum(grepl("[^0-9.]", h$comptage, useBytes = TRUE))
  h[, comptage := as.numeric(gsub("[^0-9.]", "", comptage, useBytes = TRUE))]

  # Agrégation au jour. Les deux formats de date commencent par AAAA-MM-JJ local.
  j <- h[, .(nb_velos = sum(comptage, na.rm = TRUE),
             nb_heures_mesurees = .N,
             nom_compteur = first(nom_compteur),
             id_site = first(id_site)),
         by = .(id_compteur, date = as.IDate(substr(date_heure, 1, 10)))]

  cat(sprintf("Vélo %s : %8d mesures horaires -> %6d lignes compteur x jour%s\n",
              annee, nrow(h), nrow(j),
              if (corriges > 0) sprintf(" (%d comptages '1 022' corrigés)", corriges) else ""))
  j
}

velo_archives <- rbindlist(lapply(names(archives), lire_archive))

# ===========================================================================
# 2. Comptages 2026 : jeu de données des 13 derniers mois
# ===========================================================================
# L'agrégation au jour est faite par le serveur, en heure de Paris.
url_courant <- paste0(
  "https://opendata.paris.fr/api/explore/v2.1/catalog/datasets/",
  "comptage-velo-donnees-compteurs/exports/csv?",
  "select=", URLencode("sum(sum_counts) as n, count(*) as nb_heures", reserved = TRUE),
  "&group_by=", URLencode('id_compteur, nom_compteur, id as id_site, date_format(date, "yyyy-MM-dd") as jour',
                          reserved = TRUE),
  "&timezone=Europe%2FParis&limit=-1&delimiter=%3B"
)
velo_courant <- fread(telecharger(url_courant, file.path(brut, "velo", "velo_13_derniers_mois_par_jour.csv")),
                      colClasses = list(character = c("id_compteur", "nom_compteur", "id_site")))
setnames(velo_courant, c("jour", "n", "nb_heures"), c("date", "nb_velos", "nb_heures_mesurees"))
velo_courant[, date := as.IDate(date)]
cat(sprintf("Vélo 13 derniers mois : %d lignes compteur x jour (%s -> %s)\n",
            nrow(velo_courant), min(velo_courant$date), max(velo_courant$date)))

# L'archive 2025 couvre toute l'année : on ne prend que 2026 dans ce jeu.
velo <- rbindlist(list(velo_archives[date < as.IDate("2026-01-01")],
                       velo_courant[date >= as.IDate("2026-01-01")]),
                  use.names = TRUE)
setcolorder(velo, c("date", "id_compteur", "nom_compteur", "id_site", "nb_velos", "nb_heures_mesurees"))
setorder(velo, id_compteur, date)
fwrite(velo, file.path(inter, "velo_par_compteur_et_par_jour.csv"))

# ===========================================================================
# 3. Météo quotidienne : Météo-France, station Paris-Montsouris (75114001)
# ===========================================================================
url_meteo <- "https://meteofrance.s3.sbg.io.cloud.ovh.net/data/synchro_ftp/BASE/QUOT/"
sources_meteo <- c(
  "meteo_paris_1950-2024.csv.gz" = "Q_75_previous-1950-2024_RR-T-Vent.csv.gz",
  "meteo_paris_2025-2026.csv.gz" = "Q_75_latest-2025-2026_RR-T-Vent.csv.gz"
)
meteo <- rbindlist(lapply(names(sources_meteo), function(nom) {
  fichier <- telecharger(paste0(url_meteo, sources_meteo[[nom]]), file.path(brut, "meteo", nom))
  fread(fichier, select = c("NUM_POSTE", "AAAAMMJJ", "RR", "TN", "TX", "TM", "FFM"))
}))
# Codes Météo-France : RR pluie, TN / TX / TM températures min / max / moyenne, FFM vent moyen.
meteo <- meteo[NUM_POSTE == 75114001 & AAAAMMJJ >= 20180101,
               .(date = as.IDate(as.character(AAAAMMJJ), "%Y%m%d"),
                 pluie_mm = RR, temperature_min = TN, temperature_max = TX,
                 temperature_moyenne = TM, vent_moyen_m_s = FFM)]
meteo <- unique(meteo, by = "date")[order(date)]
fwrite(meteo, file.path(inter, "meteo_paris_par_jour.csv"))
cat(sprintf("Météo : %d jours (%s -> %s)\n", nrow(meteo), min(meteo$date), max(meteo$date)))

# ===========================================================================
# 4. Calendrier : jours fériés et vacances scolaires de Paris (zone C)
# ===========================================================================
feries <- jsonlite::fromJSON(telecharger(
  "https://calendrier.api.gouv.fr/jours-feries/metropole.json",
  file.path(brut, "calendrier", "jours_feries.json")
))
feries <- data.table(date = as.IDate(names(feries)), nom_jour_ferie = unlist(feries))
fwrite(feries[order(date)], file.path(inter, "jours_feries.csv"))

url_vacances <- paste0(
  "https://data.education.gouv.fr/api/explore/v2.1/catalog/datasets/",
  "fr-en-calendrier-scolaire/exports/csv?",
  "select=", URLencode("description, start_date, end_date, annee_scolaire", reserved = TRUE),
  "&where=", URLencode('location="Paris" and population!="Enseignants"', reserved = TRUE),
  "&timezone=Europe%2FParis&limit=-1&delimiter=%3B"
)
# Dates lues comme du texte : sinon R les convertit en UTC et recule d'un jour.
vacances <- fread(telecharger(url_vacances, file.path(brut, "calendrier", "vacances_scolaires_paris.csv")),
                  colClasses = "character")
# En heure de Paris : start_date = premier jour de vacances, end_date = jour de la rentrée.
vacances <- unique(vacances[, .(nom_vacances = description,
                                annee_scolaire,
                                premier_jour = as.IDate(substr(start_date, 1, 10)),
                                jour_de_rentree = as.IDate(substr(end_date, 1, 10)))])
fwrite(vacances[order(premier_jour)], file.path(inter, "vacances_scolaires.csv"))

# ===========================================================================
# Bilan
# ===========================================================================
cat("\n=== Fichiers écrits dans", inter, "===\n")
cat(sprintf("velo_par_compteur_et_par_jour.csv : %d lignes, %d compteurs, %s -> %s\n",
            nrow(velo), uniqueN(velo$id_compteur), min(velo$date), max(velo$date)))
cat(sprintf("meteo_paris_par_jour.csv          : %d jours\n", nrow(meteo)))
cat(sprintf("jours_feries.csv                  : %d jours fériés\n", nrow(feries)))
cat(sprintf("vacances_scolaires.csv            : %d périodes\n", nrow(vacances)))
