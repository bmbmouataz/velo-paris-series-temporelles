# 02 - Nettoyage : des comptages bruts à une série journalière fiable
#
# Entrées  : data/2_intermediaire/*.csv          (produits par R/01_import.R)
# Sorties  : data/3_final/velo_paris_par_jour.csv     la série à modéliser
#            data/3_final/panel_compteurs.csv          qualité de chaque compteur, retenu ou non
#            data/3_final/jours_corriges.csv           chaque valeur remplacée, avec le motif
#            figures/02_nettoyage/*.png                une figure par décision
#
# Démarche :
#   1. Écarter la série sans identifiant (2024-2025).
#   2. Suivre chaque compteur dans le temps grâce à son nom normalisé.
#   3. Additionner les canaux d'un même compteur.
#   4. Détecter les jours anormaux (absents, incomplets, zéros, valeurs aberrantes).
#   5. Choisir un panel de compteurs stables depuis 2019.
#   6. Remplacer les jours anormaux du panel par une estimation.
#   7. Vérifier la robustesse avec un panel élargi depuis 2020.
#   8. Ajouter météo, calendrier et événements, puis écrire la série finale.
#
# À lancer depuis le dossier projet_velo/.

library(data.table)
library(ggplot2)
library(patchwork)
library(scales)

inter <- "data/2_intermediaire"
final <- "data/3_final"
fig   <- "figures/02_nettoyage"
dir.create(final, recursive = TRUE, showWarnings = FALSE)
dir.create(fig, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_minimal(base_size = 11) +
            theme(plot.title = element_text(face = "bold"),
                  plot.subtitle = element_text(colour = "grey35"),
                  panel.grid.minor = element_blank(),
                  legend.position = "bottom"))
couleurs_motif <- c("valide" = "grey80", "jour absent" = "#4C72B0",
                    "journée incomplète" = "#55A868", "zéro" = "#C44E52",
                    "anormalement bas" = "#DD8452", "anormalement haut" = "#8172B3")
sauver <- function(g, nom, l = 10, h = 5.5) {
  ggsave(file.path(fig, nom), g, width = l, height = h, dpi = 130, bg = "white")
}

debut_panel  <- as.IDate("2019-01-01")   # panel principal
debut_elargi <- as.IDate("2020-01-01")   # panel de contrôle
seuil_valide <- 0.90                     # part minimale de jours valides
seuil_derive <- 3                        # écart de niveau maximal toléré face aux autres compteurs

# Événements qui perturbent le trafic (servent aux figures et à la série finale)
evenements <- data.table(
  evenement = c("Grève des transports", "Confinement 1", "Confinement 2",
                "Confinement 3", "Jeux olympiques"),
  debut = as.IDate(c("2019-12-05", "2020-03-17", "2020-10-30", "2021-04-03", "2024-07-26")),
  fin   = as.IDate(c("2020-01-19", "2020-05-10", "2020-12-14", "2021-05-02", "2024-08-11"))
)

# ===========================================================================
# 0. Lecture
# ===========================================================================
velo <- fread(file.path(inter, "velo_par_compteur_et_par_jour.csv"),
              colClasses = list(character = c("id_compteur", "id_site")))
meteo    <- fread(file.path(inter, "meteo_paris_par_jour.csv"))
feries   <- fread(file.path(inter, "jours_feries.csv"))
vacances <- fread(file.path(inter, "vacances_scolaires.csv"))

# Le dernier jour du fichier est incomplet (données arrêtées en cours de journée).
dernier_jour <- max(velo$date) - 1L
velo <- velo[date <= dernier_jour]
cat("Période des comptages :", format(min(velo$date)), "->", format(dernier_jour), "\n")

# ===========================================================================
# 1. La série sans identifiant (2024-2025)
# ===========================================================================
sans_id <- velo[id_compteur == ""]
velo    <- velo[id_compteur != ""]

comparaison <- merge(
  sans_id[, .(date, sans_identifiant = nb_velos)],
  velo[year(date) %in% 2024:2025, .(compteurs_identifies = sum(nb_velos)), by = date],
  by = "date"
)
g <- ggplot(melt(comparaison, id.vars = "date", variable.name = "serie"),
            aes(date, value, colour = serie)) +
  geom_line(linewidth = 0.3, alpha = 0.8) +
  scale_colour_manual(values = c(sans_identifiant = "#C44E52", compteurs_identifies = "grey30"),
                      labels = c("Série sans identifiant", "Somme des compteurs identifiés"),
                      name = NULL) +
  scale_y_continuous(labels = label_number(big.mark = " ")) +
  labs(title = "1. Une série sans identifiant en 2024-2025",
       subtitle = sprintf("Aucun nom, aucune adresse, aucune coordonnée ; elle pèse %.0f %% du trafic identifié. Elle est écartée.",
                          100 * mean(comparaison$sans_identifiant / comparaison$compteurs_identifies)),
       x = NULL, y = "Vélos par jour")
sauver(g, "01_serie_sans_identifiant.png", h = 4.5)

# ===========================================================================
# 2. Suivre un compteur dans le temps : le nom normalisé
# ===========================================================================
# Les identifiants changent d'une année à l'autre et la casse des noms aussi
# ("97 Avenue Denfert..." puis "97 avenue Denfert..."). On normalise le nom :
# minuscules, sans accents, sans ponctuation, espaces uniques.
normaliser_nom <- function(x) {
  x <- tolower(stringi::stri_trans_general(x, "Latin-ASCII"))
  trimws(gsub(" +", " ", gsub("[^a-z0-9-]+", " ", x)))
}
velo[, compteur := normaliser_nom(nom_compteur)]
libelles <- velo[order(-date), .(libelle = first(nom_compteur)), by = compteur]

cat("Identifiants :", uniqueN(velo$id_compteur), "| noms bruts :", uniqueN(velo$nom_compteur),
    "| compteurs après normalisation :", uniqueN(velo$compteur), "\n")

actifs <- velo[, .(`Identifiants distincts` = uniqueN(id_compteur),
                   `Compteurs (nom normalisé)` = uniqueN(compteur)),
               by = .(mois = as.IDate(format(date, "%Y-%m-01")))]
g <- ggplot(melt(actifs, id.vars = "mois"), aes(mois, value, colour = variable)) +
  geom_step(linewidth = 0.9) +
  scale_colour_manual(values = c("#DD8452", "#4C72B0"), name = NULL) +
  labs(title = "2. Le parc de compteurs grandit, et les identifiants surestiment sa taille",
       subtitle = "En 2021-2022 et 2024-2025, un même compteur apparaît sous plusieurs identifiants (un par canal)",
       x = NULL, y = "Nombre actif dans le mois")
sauver(g, "02_compteurs_actifs_par_mois.png", h = 4.5)

# ===========================================================================
# 3. Les canaux : additionner ou garder le principal ?
# ===========================================================================
# Exemple : un compteur suivi sur toute la période, avec plusieurs canaux certaines années.
exemple <- "18 quai de l hotel de ville no-se"
canaux <- velo[compteur == exemple,
               .(nb_velos = mean(nb_velos)),
               by = .(mois = as.IDate(format(date, "%Y-%m-01")), id_compteur)]
canaux[, canal := paste("Canal", as.integer(factor(id_compteur)))]
# Barres empilées (et non geom_area, qui relierait les mois où un canal est absent).
g1 <- ggplot(canaux, aes(mois, nb_velos, fill = canal)) +
  geom_col(width = 28, alpha = 0.9) +
  scale_fill_brewer(palette = "Set2", name = NULL) +
  labs(title = sprintf("3a. « %s » : un canal principal et un canal secondaire",
                       libelles[compteur == exemple, libelle]),
       x = NULL, y = "Vélos par jour (moyenne mensuelle)")

par_jour <- velo[, .(somme = sum(nb_velos), principal = max(nb_velos), nb_canaux = .N),
                 by = .(compteur, date)]
stables <- par_jour[, .(nb_annees = uniqueN(year(date))), by = compteur][nb_annees >= 7, compteur]
comparatif <- par_jour[compteur %in% stables,
                       .(`Somme des canaux` = mean(somme), `Canal principal seul` = mean(principal),
                         plusieurs_canaux = mean(nb_canaux) > 1.05),
                       by = .(annee = year(date))]
g2 <- ggplot(melt(comparatif, id.vars = c("annee", "plusieurs_canaux")),
             aes(annee, value, colour = variable)) +
  geom_rect(data = comparatif[plusieurs_canaux == TRUE],
            aes(xmin = annee - 0.5, xmax = annee + 0.5, ymin = -Inf, ymax = Inf),
            inherit.aes = FALSE, fill = "#FFF3D6") +
  geom_line(linewidth = 0.9) + geom_point(size = 2) +
  scale_colour_manual(values = c("#4C72B0", "#C44E52"), name = NULL) +
  scale_x_continuous(breaks = 2018:2026) +
  labs(title = sprintf("3b. Moyenne des %d compteurs présents 7 ans ou plus", length(stables)),
       subtitle = "En jaune : années à plusieurs canaux. Le canal principal seul crée de faux creux ;\nla somme est continue avec les années à un seul canal (2023, 2026). On additionne.",
       x = NULL, y = "Vélos par jour et par compteur")
sauver(g1 / g2, "03_canaux_somme_ou_principal.png", h = 8.5)

# Une ligne par compteur et par jour, sur une grille complète (les jours absents apparaissent).
jour <- velo[, .(nb_velos = sum(nb_velos), nb_heures = min(nb_heures_mesurees)),
             by = .(compteur, date)]
# La grille va jusqu'au dernier jour des données : un compteur déposé en cours de
# route apparaît alors en « jour absent » au lieu de disparaître sans bruit.
grille <- jour[, .(date = seq(min(date), dernier_jour, by = 1L)), by = compteur]
jour <- merge(grille, jour, by = c("compteur", "date"), all.x = TRUE)
setorder(jour, compteur, date)

# ===========================================================================
# 4. Détecter les jours anormaux
# ===========================================================================
# Chaque compteur est vu comme : nb_velos = niveau propre x facteur commun du jour x bruit.
#   - niveau propre : médiane glissante sur 29 jours des jours sains du compteur ;
#   - facteur commun : médiane, sur tous les compteurs, du rapport nb_velos / niveau.
# Un jour est aberrant si le compteur s'écarte du facteur commun d'un facteur > 3
# (bas ou haut). Un confinement ou un 25 décembre touchent tous les compteurs à la
# fois : ils modifient le facteur commun et ne sont donc pas signalés à tort.
niveau_glissant <- function(x, ok) {
  n <- frollapply(fifelse(ok, x, NA_real_), 29, median, na.rm = TRUE, align = "center")
  zoo::na.approx(n, na.rm = FALSE, rule = 2)
}
jour[, base_ok := !is.na(nb_velos) & nb_heures >= 23 & nb_velos > 0]
jour[, sain := base_ok]
for (passage in 1:2) {   # 2e passage : le niveau est recalculé sans les aberrations trouvées
  jour[, niveau := niveau_glissant(nb_velos, sain), by = compteur]
  jour[, rapport := nb_velos / niveau]
  jour[, facteur_commun := median(rapport[sain & niveau >= 30], na.rm = TRUE), by = date]
  jour[, ecart := rapport / facteur_commun]
  jour[, aberrant := base_ok & niveau >= 30 & (ecart < 1 / 3 | ecart > 3)]
  jour[, sain := base_ok & !aberrant]
}
jour[, motif := fcase(is.na(nb_velos), "jour absent",
                      nb_heures < 23, "journée incomplète",
                      nb_velos == 0, "zéro",
                      aberrant & ecart < 1, "anormalement bas",
                      aberrant, "anormalement haut",
                      default = "valide")]
jour[, motif := factor(motif, levels = names(couleurs_motif))]
cat(sprintf("Jours anormaux : %.1f %% de l'ensemble\n", 100 * mean(jour$motif != "valide")))

# 4a. Répartition des anomalies par année
repartition <- jour[motif != "valide", .N, by = .(annee = year(date), motif)]
g <- ggplot(repartition, aes(factor(annee), N, fill = motif)) +
  geom_col() +
  scale_fill_manual(values = couleurs_motif, name = NULL, drop = TRUE) +
  scale_y_continuous(labels = label_number(big.mark = " ")) +
  labs(title = "4a. Les anomalies par année et par type (tous compteurs)",
       subtitle = sprintf("%.1f %% des jours compteur sont anormaux. Les jours absents dominent : pannes longues, compteurs déposés.",
                          100 * mean(jour$motif != "valide")),
       x = NULL, y = "Jours compteur")
sauver(g, "04a_anomalies_par_annee.png", h = 4.8)

# 4b. Le facteur commun : il capte les chocs qui touchent tout Paris
fc <- unique(jour[, .(date, facteur_commun)])[!is.na(facteur_commun)]
g <- ggplot(fc, aes(date, facteur_commun)) +
  geom_rect(data = evenements, aes(xmin = debut, xmax = fin, ymin = -Inf, ymax = Inf, fill = evenement),
            inherit.aes = FALSE, alpha = 0.25) +
  geom_line(linewidth = 0.25, colour = "grey30") +
  geom_hline(yintercept = 1, linetype = 2, colour = "grey50") +
  scale_fill_brewer(palette = "Set1", name = NULL) +
  labs(title = "4b. Le facteur commun du jour : trafic de Paris rapporté au niveau habituel de chaque compteur",
       subtitle = "Week-ends, fériés, pluie, confinements font varier tous les compteurs ensemble : ce ne sont pas des anomalies",
       x = NULL, y = "Facteur commun (1 = jour habituel)")
sauver(g, "04b_facteur_commun.png", h = 4.8)

# 4c. Exemple : le compteur avec le plus de jours aberrants détectés
pire <- jour[date >= debut_panel, .(n = sum(motif %in% c("anormalement bas", "anormalement haut", "zéro"))),
             by = compteur][order(-n)][1, compteur]
ex <- jour[compteur == pire & date >= debut_panel]
g <- ggplot(ex, aes(date, nb_velos)) +
  geom_line(linewidth = 0.25, colour = "grey40") +
  geom_line(aes(y = niveau * facteur_commun), colour = "#4C72B0", linewidth = 0.3, alpha = 0.7) +
  geom_point(data = ex[motif != "valide" & !is.na(nb_velos)], aes(colour = motif), size = 1.1) +
  scale_colour_manual(values = couleurs_motif, name = NULL, drop = TRUE) +
  labs(title = sprintf("4c. Exemple de détection : « %s »", libelles[compteur == pire, libelle]),
       subtitle = "Gris : comptage observé. Bleu : valeur attendue (niveau du compteur x facteur commun). Points : jours signalés.",
       x = NULL, y = "Vélos par jour")
sauver(g, "04c_exemple_detection.png", h = 4.8)

# ===========================================================================
# 5. Choisir le panel de compteurs
# ===========================================================================
qualite_panel <- function(debut) {
  periode <- seq(debut, dernier_jour, by = 1L)
  q <- jour[date >= debut, .(part_valide = sum(motif == "valide") / length(periode),
                             premier_jour = min(date), dernier_valide = max(date[motif == "valide"]),
                             niveau_moyen = mean(nb_velos[motif == "valide"])),
            by = compteur]
  premiers <- jour[, .(premier_jour = min(date)), by = compteur]
  q[premiers, on = "compteur", premier_jour := i.premier_jour]
  # Dérive : le niveau annuel du compteur, rapporté à son propre niveau moyen, ne doit
  # pas s'écarter de plus d'un facteur 3 de celui des autres compteurs (déplacement,
  # capteur à moitié en panne...).
  annuel <- jour[date >= debut & motif == "valide", .(niveau = median(nb_velos)),
                 by = .(compteur, annee = year(date))]
  annuel[, forme := niveau / exp(mean(log(niveau))), by = compteur]
  annuel[compteur %in% q[part_valide >= seuil_valide, compteur],
         forme_mediane := median(forme), by = annee]
  annuel[, forme_mediane := forme_mediane[!is.na(forme_mediane)][1], by = annee]
  annuel[, ecart_relatif := forme / forme_mediane]
  derive <- annuel[, .(derive_max = max(pmax(ecart_relatif, 1 / ecart_relatif))), by = compteur]
  q <- merge(q, derive, by = "compteur", all.x = TRUE)
  # Le compteur doit aussi fonctionner jusqu'à la fin : les derniers mois serviront
  # de période de test aux modèles, on n'y veut pas de valeurs estimées en bloc.
  q[, actif_a_la_fin := dernier_valide >= dernier_jour - 30L]
  q[, retenu := premier_jour <= debut & part_valide >= seuil_valide & derive_max <= seuil_derive & actif_a_la_fin]
  q[, raison := fcase(premier_jour > debut, "installé trop tard",
                      !actif_a_la_fin, "arrêté avant la fin",
                      part_valide < seuil_valide, sprintf("seulement %.0f %% de jours valides", 100 * part_valide),
                      derive_max > seuil_derive, sprintf("niveau qui dérive (x%.1f)", derive_max),
                      default = "retenu")]
  list(qualite = q, annuel = annuel)
}

principal <- qualite_panel(debut_panel)
elargi    <- qualite_panel(debut_elargi)
panel     <- principal$qualite[retenu == TRUE, compteur]
panel_el  <- elargi$qualite[retenu == TRUE, compteur]
cat(sprintf("Panel principal (depuis %s) : %d compteurs | panel élargi (depuis %s) : %d compteurs\n",
            debut_panel, length(panel), debut_elargi, length(panel_el)))

# 5a. Calendrier de présence et de qualité de tous les compteurs
presence <- jour[, .(part_valide = mean(motif == "valide")),
                 by = .(compteur, mois = as.IDate(format(date, "%Y-%m-01")))]
ordre <- jour[, .(debut = min(date)), by = compteur][order(debut, compteur), compteur]
presence[, compteur := factor(compteur, levels = rev(ordre))]
presence[, statut := factor(fifelse(compteur %in% panel, sprintf("Panel (%d)", length(panel)), "Écartés"),
                            levels = c(sprintf("Panel (%d)", length(panel)), "Écartés"))]
g <- ggplot(presence, aes(mois + 14L, compteur, fill = part_valide)) +
  geom_tile(width = 31, height = 0.9) +
  geom_vline(xintercept = debut_panel, colour = "#C44E52", linewidth = 0.6) +
  facet_grid(statut ~ ., scales = "free_y", space = "free_y") +
  scale_fill_viridis_c(option = "mako", direction = -1, labels = percent, name = "Jours valides dans le mois") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y", expand = c(0, 0)) +
  labs(title = sprintf("5a. Présence et qualité des %d compteurs, mois par mois", length(ordre)),
       subtitle = sprintf("Une ligne par compteur, triés par date d'installation. Clair = panne. Trait rouge : début du panel (%s).",
                          format(debut_panel, "%d/%m/%Y")),
       x = NULL, y = NULL) +
  theme(axis.text.y = element_blank(), panel.grid = element_blank(),
        strip.text.y = element_text(angle = 0, face = "bold"), legend.key.width = unit(1.5, "cm"))
sauver(g, "05a_presence_compteurs.png", l = 11, h = 11)

# 5b. Pourquoi les compteurs sont retenus ou non
raisons <- principal$qualite[, .N, by = .(raison = fifelse(grepl("jours valides", raison), "trop de jours anormaux",
                                                   fifelse(grepl("dérive", raison), "niveau qui dérive", raison)))]
g1 <- ggplot(raisons, aes(reorder(raison, N), N, fill = raison == "retenu")) +
  geom_col(show.legend = FALSE) + coord_flip() +
  geom_text(aes(label = N), hjust = -0.2) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("grey60", "#4C72B0")) +
  labs(title = "5b. Sélection du panel principal", x = NULL, y = "Nombre de compteurs")
g2 <- ggplot(principal$qualite[premier_jour <= debut_panel & actif_a_la_fin == TRUE],
             aes(part_valide, derive_max, colour = retenu)) +
  annotate("rect", xmin = seuil_valide, xmax = 1, ymin = 1, ymax = seuil_derive, fill = "#E8F0FA") +
  geom_point(size = 2.2, alpha = 0.85) +
  geom_vline(xintercept = seuil_valide, linetype = 2) +
  geom_hline(yintercept = seuil_derive, linetype = 2) +
  scale_x_continuous(labels = percent) + scale_y_log10() +
  scale_colour_manual(values = c("grey55", "#4C72B0"), labels = c("écarté", "retenu"), name = NULL) +
  labs(subtitle = "Compteurs installés avant 2019 et encore actifs : zone bleue = critères respectés",
       x = "Part de jours valides depuis 2019", y = "Dérive maximale du niveau (échelle log)")
sauver(g1 | g2, "05b_selection_panel.png", l = 11, h = 4.8)

# 5c. Niveaux annuels : les compteurs retenus évoluent ensemble
an <- principal$annuel[compteur %in% principal$qualite[premier_jour <= debut_panel, compteur]]
an[, statut := fifelse(compteur %in% panel, "retenu", "écarté")]
g <- ggplot(an, aes(annee, forme, group = compteur, colour = statut)) +
  geom_line(data = an[statut == "écarté"], alpha = 0.6) +
  geom_line(data = an[statut == "retenu"], alpha = 0.5) +
  stat_summary(data = an[statut == "retenu"], aes(group = 1), fun = median,
               geom = "line", colour = "black", linewidth = 1.2) +
  scale_y_log10() + scale_x_continuous(breaks = 2019:2026) +
  scale_colour_manual(values = c("écarté" = "#DD8452", "retenu" = "#4C72B0"), name = NULL) +
  labs(title = "5c. Niveau annuel de chaque compteur, rapporté à sa propre moyenne",
       subtitle = "Noir : médiane du panel. Les compteurs écartés (orange) ont des sauts sans rapport avec la tendance commune.",
       x = NULL, y = "Niveau relatif (échelle log)")
sauver(g, "05c_niveaux_annuels.png", h = 5)

# ===========================================================================
# 6. Remplacer les jours anormaux du panel
# ===========================================================================
# Valeur de remplacement = niveau du compteur x facteur commun du jour, ce dernier
# étant recalculé sur les seuls compteurs du panel.
construire_serie <- function(compteurs, debut) {
  p <- jour[compteur %in% compteurs & date >= debut]
  p[, niveau := niveau_glissant(nb_velos, motif == "valide"), by = compteur]
  p[, facteur_panel := median((nb_velos / niveau)[motif == "valide"], na.rm = TRUE), by = date]
  # Jours où aucun compteur du panel ne fonctionne (coupure générale du système) :
  # facteur = moyenne du même jour de la semaine, 7 jours avant et 7 jours après.
  f <- unique(p[, .(date, facteur_panel)])[order(date)]
  f[, coupure_generale := is.na(facteur_panel)]
  f[, facteur_panel := fcoalesce(facteur_panel,
                                 rowMeans(cbind(shift(facteur_panel, 7), shift(facteur_panel, -7)), na.rm = TRUE))]
  p[, c("facteur_panel") := NULL]
  p <- merge(p, f, by = "date")
  p[, estime := round(niveau * facteur_panel)]
  p[, retenu_val := fifelse(motif == "valide", as.numeric(nb_velos), estime)]
  p
}
p  <- construire_serie(panel, debut_panel)
pe <- construire_serie(panel_el, debut_elargi)

serie <- p[, .(velos = sum(retenu_val), part_estimee = sum(retenu_val[motif != "valide"]) / sum(retenu_val),
               compteurs_valides = sum(motif == "valide"), coupure_generale = first(coupure_generale)), by = date]
cat("Jours de coupure générale (aucun compteur du panel valide) :",
    paste(format(serie[coupure_generale == TRUE, date]), collapse = ", "), "\n")

# 6a. Calendrier de qualité du panel, jour par jour
p[, compteur_f := factor(compteur, levels = p[, .(m = mean(retenu_val)), by = compteur][order(m), compteur])]
g <- ggplot(p, aes(date, compteur_f, fill = motif)) +
  geom_tile() +
  scale_fill_manual(values = couleurs_motif, name = NULL, drop = TRUE) +
  scale_y_discrete(labels = function(x) libelles[match(x, compteur), libelle]) +
  labs(title = sprintf("6a. Les %d compteurs du panel, jour par jour", length(panel)),
       subtitle = sprintf("Chaque case colorée est remplacée par une estimation : %.1f %% des jours compteur.",
                          100 * mean(p$motif != "valide")),
       x = NULL, y = NULL) +
  theme(axis.text.y = element_text(size = 7), panel.grid = element_blank())
sauver(g, "06a_panel_jour_par_jour.png", l = 11, h = 7)

# 6b. Exemple de remplacement sur la plus longue panne du panel
p[, bloc := rleid(motif != "valide"), by = compteur]
pannes <- p[motif != "valide", .(debut = min(date), fin = max(date), duree = .N), by = .(compteur, bloc)][order(-duree)]
pp <- pannes[1]
fenetre <- p[compteur == pp$compteur & date >= pp$debut - 45 & date <= pp$fin + 45]
g <- ggplot(fenetre, aes(date)) +
  annotate("rect", xmin = pp$debut, xmax = pp$fin, ymin = -Inf, ymax = Inf, fill = "#FBE3E3") +
  geom_line(aes(y = nb_velos, colour = "Comptage observé"), linewidth = 0.5, na.rm = TRUE) +
  geom_line(aes(y = estime, colour = "Estimation"), linewidth = 0.5, linetype = 2) +
  geom_point(data = fenetre[motif != "valide"], aes(y = retenu_val, colour = "Valeur retenue"), size = 1.2) +
  scale_colour_manual(values = c("Comptage observé" = "grey25", "Estimation" = "#4C72B0", "Valeur retenue" = "#C44E52"), name = NULL) +
  labs(title = sprintf("6b. Remplacement d'une panne de %d jours : « %s »", pp$duree, libelles[compteur == pp$compteur, libelle]),
       subtitle = "L'estimation suit le rythme des autres compteurs (week-ends, météo) au niveau propre de ce compteur",
       x = NULL, y = "Vélos par jour")
sauver(g, "06b_exemple_remplacement.png", h = 4.8)

# 6c. Part du trafic estimée, mois par mois
mensuel <- serie[, .(part_estimee = sum(velos * part_estimee) / sum(velos)), by = .(mois = as.IDate(format(date, "%Y-%m-01")))]
g <- ggplot(mensuel, aes(mois, part_estimee)) +
  geom_col(fill = "#4C72B0") +
  geom_hline(yintercept = mean(serie$part_estimee), linetype = 2, colour = "#C44E52") +
  scale_y_continuous(labels = percent) +
  labs(title = "6c. Part du trafic du panel qui est estimée, et non observée",
       subtitle = sprintf("Moyenne : %.1f %% (pointillé). Le reste est mesuré directement.", 100 * mean(serie$part_estimee)),
       x = NULL, y = "Part estimée")
sauver(g, "06c_part_estimee.png", h = 4.3)

# ===========================================================================
# 7. Robustesse : le panel élargi raconte-t-il la même histoire ?
# ===========================================================================
serie_el <- pe[, .(velos = sum(retenu_val)), by = date]
robustesse <- merge(serie[, .(date, principal = velos)], serie_el[, .(date, elargi = velos)], by = "date")
base <- robustesse[year(date) == 2021, lapply(.SD, mean), .SDcols = c("principal", "elargi")]
robustesse[, `:=`(principal = 100 * principal / base$principal, elargi = 100 * elargi / base$elargi)]
robustesse[, `:=`(principal_28j = frollmean(principal, 28, align = "center"),
                  elargi_28j = frollmean(elargi, 28, align = "center"))]
correlation <- robustesse[, cor(principal, elargi)]
g1 <- ggplot(robustesse, aes(date)) +
  geom_line(aes(y = principal_28j, colour = sprintf("Panel principal (%d compteurs, depuis 2019)", length(panel))), linewidth = 0.7, na.rm = TRUE) +
  geom_line(aes(y = elargi_28j, colour = sprintf("Panel élargi (%d compteurs, depuis 2020)", length(panel_el))), linewidth = 0.7, na.rm = TRUE) +
  scale_colour_manual(values = c("#4C72B0", "#DD8452"), name = NULL) +
  labs(title = "7. Robustesse : les deux panels ont la même dynamique",
       subtitle = sprintf("Indice base 100 = moyenne 2021, moyenne glissante 28 jours. Corrélation journalière : %.3f.\nLa hausse est un peu plus forte pour le panel élargi : ses compteurs, plus récents, sont sur de nouveaux axes cyclables.", correlation),
       x = NULL, y = "Indice")
g2 <- ggplot(robustesse, aes(elargi, principal)) +
  geom_point(alpha = 0.15, size = 0.6, colour = "#4C72B0") +
  geom_abline(linetype = 2) +
  labs(x = "Panel élargi (indice)", y = "Panel principal (indice)")
sauver(g1 + g2 + plot_layout(widths = c(2.3, 1)), "07_robustesse_panels.png", l = 11, h = 4.8)

# ===========================================================================
# 8. Série finale : trafic + météo + calendrier + événements
# ===========================================================================
fin_serie <- min(dernier_jour, max(meteo$date))
finale <- serie[date <= fin_serie]
finale <- merge(finale, meteo, by = "date", all.x = TRUE)

manquants_meteo <- finale[, lapply(.SD, function(x) sum(is.na(x))),
                          .SDcols = c("pluie_mm", "temperature_moyenne", "temperature_min", "temperature_max", "vent_moyen_m_s")]
cat("Valeurs météo manquantes avant correction :\n"); print(manquants_meteo)
# Température et vent : interpolation linéaire (variables lisses). Pluie : 0 si manquante
# n'est pas défendable, on interpole aussi mais on garde la trace.
finale[, meteo_interpolee := is.na(pluie_mm) | is.na(temperature_moyenne)]
for (v in c("pluie_mm", "temperature_moyenne", "temperature_min", "temperature_max", "vent_moyen_m_s")) {
  set(finale, j = v, value = zoo::na.approx(finale[[v]], na.rm = FALSE, rule = 2))
}

jours_fr <- c("lundi", "mardi", "mercredi", "jeudi", "vendredi", "samedi", "dimanche")
# wday() : 1 = dimanche ... 7 = samedi ; on décale pour que 1 = lundi.
finale[, jour_semaine := factor(jours_fr[((wday(date) + 5) %% 7) + 1], levels = jours_fr)]
finale[, est_week_end := jour_semaine %in% c("samedi", "dimanche")]
finale <- merge(finale, feries, by = "date", all.x = TRUE)
finale[, est_ferie := !is.na(nom_jour_ferie)]
finale[vacances, on = .(date >= premier_jour, date < jour_de_rentree), nom_vacances := i.nom_vacances]
finale[, est_vacances := !is.na(nom_vacances)]
finale[evenements, on = .(date >= debut, date <= fin), evenement := i.evenement]
setorder(finale, date)

setcolorder(finale, c("date", "velos", "jour_semaine", "est_week_end", "est_ferie", "nom_jour_ferie",
                      "est_vacances", "nom_vacances", "evenement",
                      "temperature_moyenne", "temperature_min", "temperature_max", "pluie_mm", "vent_moyen_m_s",
                      "meteo_interpolee", "part_estimee", "compteurs_valides", "coupure_generale"))
fwrite(finale, file.path(final, "velo_paris_par_jour.csv"))

# 8a. La série finale et ses événements
finale[, velos_28j := frollmean(velos, 28, align = "center")]
g <- ggplot(finale, aes(date)) +
  geom_rect(data = evenements, aes(xmin = debut, xmax = fin, ymin = -Inf, ymax = Inf, fill = evenement),
            inherit.aes = FALSE, alpha = 0.3) +
  geom_line(aes(y = velos), linewidth = 0.2, colour = "grey55") +
  geom_line(aes(y = velos_28j), linewidth = 0.9, colour = "black", na.rm = TRUE) +
  scale_fill_brewer(palette = "Set1", name = NULL) +
  scale_y_continuous(labels = label_number(big.mark = " ")) +
  labs(title = sprintf("8a. Série finale : vélos comptés chaque jour par les %d compteurs du panel", length(panel)),
       subtitle = sprintf("%s -> %s, %d jours. Gris : jour par jour. Noir : moyenne glissante 28 jours.",
                          format(min(finale$date), "%d/%m/%Y"), format(max(finale$date), "%d/%m/%Y"), nrow(finale)),
       x = NULL, y = "Vélos par jour")
sauver(g, "08a_serie_finale.png", l = 11, h = 5)

# 8b. Calendrier du trafic : une case par jour
cal <- copy(finale)[, `:=`(annee = year(date), semaine = as.integer(format(date, "%W")),
                           jour = factor(jour_semaine, levels = rev(jours_fr)))]
g <- ggplot(cal, aes(semaine, jour, fill = velos)) +
  geom_tile(colour = "white", linewidth = 0.15) +
  facet_wrap(~annee, ncol = 1, strip.position = "left") +
  scale_fill_viridis_c(option = "inferno", labels = label_number(big.mark = " "), name = "Vélos par jour") +
  scale_x_continuous(breaks = c(1, 9, 18, 27, 36, 45), labels = c("janv.", "mars", "mai", "juil.", "sept.", "nov."), expand = c(0, 0)) +
  labs(title = "8b. Le calendrier du vélo à Paris",
       subtitle = "Une case = un jour. On lit les week-ends, l'été, Noël, les confinements de 2020 et la hausse d'année en année.",
       x = NULL, y = NULL) +
  theme(axis.text.y = element_text(size = 6), panel.grid = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold"), legend.key.width = unit(1.5, "cm"))
sauver(g, "08b_calendrier_trafic.png", l = 11, h = 10)

# 8c. Contrôle de la météo et du calendrier joints à la série
g1 <- ggplot(finale, aes(date, temperature_moyenne)) +
  geom_line(linewidth = 0.25, colour = "#C44E52") +
  geom_point(data = finale[meteo_interpolee == TRUE], colour = "black", size = 0.8) +
  labs(title = "8c. Variables jointes à la série", subtitle = "Points noirs : jours de météo interpolés", x = NULL, y = "Température moyenne (°C)")
g2 <- ggplot(finale, aes(date, pluie_mm)) +
  geom_col(width = 1, fill = "#4C72B0") +
  labs(x = NULL, y = "Pluie (mm)")
g3 <- ggplot(finale, aes(date, 1)) +
  geom_tile(data = finale[est_vacances == TRUE], aes(fill = "Vacances scolaires"), height = 0.8) +
  geom_tile(data = finale[est_ferie == TRUE], aes(fill = "Jour férié"), height = 0.8, width = 3) +
  scale_fill_manual(values = c("Vacances scolaires" = "#55A868", "Jour férié" = "#C44E52"), name = NULL) +
  labs(x = NULL, y = NULL) + theme(axis.text.y = element_blank(), panel.grid = element_blank())
sauver(g1 / g2 / g3 + plot_layout(heights = c(2, 2, 0.6)), "08c_meteo_calendrier.png", l = 11, h = 7)

# ===========================================================================
# Sorties de contrôle
# ===========================================================================
panel_info <- merge(principal$qualite, libelles, by = "compteur")
panel_info[, dans_panel_elargi := compteur %in% panel_el]
setcolorder(panel_info, c("libelle", "compteur", "retenu", "raison", "premier_jour", "part_valide", "derive_max", "niveau_moyen"))
fwrite(panel_info[order(-retenu, libelle)], file.path(final, "panel_compteurs.csv"))

corriges <- merge(p[motif != "valide", .(date, compteur, motif, valeur_observee = nb_velos, valeur_retenue = retenu_val)],
                  libelles, by = "compteur")
setcolorder(corriges, c("date", "libelle", "compteur"))
fwrite(corriges[order(compteur, date)], file.path(final, "jours_corriges.csv"))

cat("\n=== Bilan du nettoyage ===\n")
cat(sprintf("Série finale         : %d jours, %s -> %s\n", nrow(finale), min(finale$date), max(finale$date)))
cat(sprintf("Panel principal      : %d compteurs (sur %d installés avant 2019)\n",
            length(panel), principal$qualite[premier_jour <= debut_panel, .N]))
cat(sprintf("Jours compteur remplacés : %d (%.1f %%), trafic estimé : %.1f %%\n",
            nrow(corriges), 100 * mean(p$motif != "valide"), 100 * mean(serie$part_estimee)))
cat(sprintf("Corrélation panel principal / élargi : %.3f\n", correlation))
cat(sprintf("Vélos par jour : moyenne %.0f, min %.0f (%s), max %.0f (%s)\n",
            mean(finale$velos), min(finale$velos), finale[which.min(velos), date],
            max(finale$velos), finale[which.max(velos), date]))
cat("Figures :", length(list.files(fig, pattern = "png$")), "dans", fig, "\n")
