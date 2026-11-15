# 03 - Exploration de la série
#
# Entrée  : data/3_final/velo_paris_par_jour.csv
# Sorties : figures/03_exploration/*.png, resultats/03_*.csv
#
# Questions posées :
#   1. À quoi ressemble la série (tendance, saisons, chocs) ?
#   2. Faut-il passer au log (variance qui croît avec le niveau, comptages) ?
#   3. Quels cycles : semaine, année ? (profils, ACF, périodogramme)
#   4. Quel est l'effet de la météo, du calendrier et des événements ?
#
# Pour isoler un effet (pluie, férié...) sans être gêné par la tendance et la
# semaine, on calcule un « trafic relatif » : le trafic du jour divisé par la
# médiane des mêmes jours de semaine des 3 semaines avant et des 3 semaines après.
# Un trafic relatif de 0,8 = 20 % de vélos en moins qu'un jour comparable.

source("R/00_outils.R")
preparer_figures("03_exploration")
d <- charger_serie()

# Trafic relatif (voir en-tête)
voisins <- c(-21L, -14L, -7L, 7L, 14L, 21L)
ref <- do.call(cbind, lapply(voisins, function(k) shift(d$velos, -k)))
ref[ref <= 0] <- NA
d[, reference := apply(ref, 1, median, na.rm = TRUE)]
d[, relatif := velos / reference]
d_normal <- d[evenement == "" & !coupure_generale]   # jours sans événement exceptionnel

# ===========================================================================
# 1. Vue d'ensemble
# ===========================================================================
d[, moyenne_7  := frollmean(velos, 7,  align = "center")]
d[, moyenne_365 := frollmean(velos, 365, align = "center")]
g <- ggplot(d, aes(as.Date(date))) +
  bandes_evenements() +
  geom_line(aes(y = velos), colour = "grey70", linewidth = 0.25) +
  geom_line(aes(y = moyenne_7, colour = "Moyenne mobile 7 jours"), linewidth = 0.5) +
  geom_line(aes(y = moyenne_365, colour = "Moyenne mobile 365 jours"), linewidth = 1.1) +
  scale_colour_manual(values = c("Moyenne mobile 7 jours" = bleu, "Moyenne mobile 365 jours" = rouge), name = NULL) +
  scale_fill_manual(values = couleurs_evenement, name = NULL) +
  scale_y_continuous(labels = milliers) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Vélos comptés chaque jour à Paris (27 compteurs), 2019-2026",
       subtitle = "Gris : jour par jour. La moyenne sur 7 jours efface la semaine, celle sur 365 jours ne garde que la tendance",
       x = NULL, y = "Vélos par jour") +
  guides(fill = guide_legend(nrow = 1), colour = guide_legend(nrow = 1)) +
  theme(legend.box = "vertical")
sauver(g, "01_serie_complete.png", h = 6)

# Une année zoomée : le motif hebdomadaire saute aux yeux
annee_zoom <- d[date >= as.IDate("2024-09-01") & date <= as.IDate("2024-12-31")]
g <- ggplot(annee_zoom, aes(as.Date(date), velos)) +
  geom_line(colour = "grey50") +
  geom_point(aes(colour = est_week_end), size = 1.8) +
  geom_point(data = annee_zoom[est_ferie == TRUE], shape = 21, size = 4, colour = rouge, stroke = 1) +
  scale_colour_manual(values = c("FALSE" = bleu, "TRUE" = orange), labels = c("En semaine", "Week-end"), name = NULL) +
  scale_y_continuous(labels = milliers) +
  scale_x_date(date_breaks = "2 weeks", date_labels = "%d %b") +
  labs(title = "Zoom sur l'automne 2024 : un cycle de 7 jours très régulier",
       subtitle = "Cercles rouges : jours fériés (1er et 11 novembre, Noël). Baisse nette pendant les vacances de Noël",
       x = NULL, y = "Vélos par jour")
sauver(g, "02_zoom_automne_2024.png")

# ===========================================================================
# 2. Pourquoi le log ? modèle multiplicatif, comptages type Poisson
# ===========================================================================
# Si la dispersion augmente avec le niveau, le modèle est multiplicatif :
# X_t = m_t × s_t × z_t, et log X_t = log m_t + log s_t + log z_t est additif.
lambda <- BoxCox.lambda(ts(d_normal$velos, frequency = 7), method = "guerrero")
d[, semaine := as.IDate(cut(as.Date(date), "week"))]
par_semaine <- d_normal[, .(moyenne = mean(velos), ecart_type = sd(velos),
                            moyenne_log = mean(log_velos), ecart_type_log = sd(log_velos), n = .N),
                        by = .(semaine = as.IDate(cut(as.Date(date), "week")))][n == 7]
pente_brut <- coef(lm(log(ecart_type) ~ log(moyenne), par_semaine))[2]
g1 <- ggplot(par_semaine, aes(moyenne, ecart_type)) +
  geom_point(alpha = 0.5, colour = bleu) + geom_smooth(method = "lm", formula = y ~ x, colour = rouge) +
  scale_x_continuous(labels = milliers) + scale_y_continuous(labels = milliers) +
  labs(title = "Série brute", subtitle = sprintf("L'écart-type croît avec la moyenne (élasticité %.2f)", pente_brut),
       x = "Moyenne de la semaine", y = "Écart-type dans la semaine")
g2 <- ggplot(par_semaine, aes(moyenne_log, ecart_type_log)) +
  geom_point(alpha = 0.5, colour = vert) + geom_smooth(method = "lm", formula = y ~ x, colour = rouge) +
  labs(title = "Série en log", subtitle = "Légère sur-correction : la dispersion baisse un peu avec le niveau",
       x = "Moyenne du log", y = "Écart-type du log")
g3 <- ggplot(d, aes(velos)) + geom_histogram(bins = 70, fill = bleu, alpha = 0.8) +
  scale_x_continuous(labels = milliers) + labs(title = "Distribution brute", x = "Vélos par jour", y = NULL)
g4 <- ggplot(d, aes(log_velos)) + geom_histogram(bins = 70, fill = vert, alpha = 0.8) +
  labs(title = "Distribution du log", subtitle = "Queue gauche : confinements", x = "log(vélos)", y = NULL)
g <- (g1 | g2) / (g3 | g4) +
  plot_annotation(title = "Variance et niveau : pourquoi modéliser log(vélos)",
                  subtitle = sprintf("Une semaine = un point. Élasticité 0,5 = Poisson (racine), 1 = log. Box-Cox optimal : lambda = %.2f, entre les deux ; le log est retenu car ses effets se lisent en %%", lambda))
sauver(g, "03_justification_log.png", h = 8)

# Choix de Box-Cox par la méthode de Guerrero : lambda proche de 0 = log
cat(sprintf("Box-Cox (Guerrero) : lambda = %.2f (0 = log)\n", lambda))

# ===========================================================================
# 3. La saison hebdomadaire
# ===========================================================================
g1 <- ggplot(d_normal[est_ferie == FALSE], aes(jour_semaine, velos, fill = est_week_end)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.4) +
  scale_fill_manual(values = c("FALSE" = bleu, "TRUE" = orange), guide = "none") +
  scale_y_continuous(labels = milliers) +
  labs(title = "Par jour de la semaine", x = NULL, y = "Vélos par jour")
profil <- d_normal[!est_ferie & !(type_vacances %in% c("Été", "Noël")),
                   .(velos = mean(velos)), by = .(annee, jour_semaine)]
profil[, indice := velos / mean(velos), by = annee]
g2 <- ggplot(profil, aes(jour_semaine, indice, colour = factor(annee), group = annee)) +
  geom_line(linewidth = 0.9) + geom_point() +
  scale_colour_viridis_d(name = NULL, end = 0.9) +
  labs(title = "Profil hebdomadaire par année (1 = moyenne de l'année)",
       subtitle = "Hors vacances d'été et de Noël. Profil presque identique chaque année : mardi-jeudi au sommet, dimanche -30 %",
       x = NULL, y = "Indice") + guides(colour = guide_legend(nrow = 1))
ratio <- d_normal[est_ferie == FALSE, .(semaine = mean(velos[!est_week_end]), week_end = mean(velos[est_week_end])),
                  by = .(mois = as.IDate(format(date, "%Y-%m-01")))]
g3 <- ggplot(ratio, aes(as.Date(mois), week_end / semaine)) +
  geom_line(colour = violet) + geom_smooth(method = "loess", formula = y ~ x, span = 0.4, colour = rouge, se = FALSE) +
  scale_y_continuous(labels = percent) + scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Trafic du week-end rapporté à la semaine, chaque mois",
       subtitle = "Stable autour de 65 % : le week-end compte environ deux tiers d'un jour ouvré (pic en 2020 : confinements)",
       x = NULL, y = "Week-end / semaine")
g <- (g1 | g2) / g3 + plot_annotation(title = "Saison hebdomadaire : période 7")
sauver(g, "04_saison_hebdomadaire.png", h = 9)

# ===========================================================================
# 4. La saison annuelle
# ===========================================================================
d[, lisse_7 := frollmean(velos, 7, align = "center")]
g1 <- ggplot(d[annee >= 2019], aes(jour_annee, lisse_7, colour = factor(annee), group = annee)) +
  geom_line(linewidth = 0.7) +
  scale_colour_viridis_d(name = NULL, end = 0.92) +
  scale_x_continuous(breaks = cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30)), labels = mois_fr) +
  scale_y_continuous(labels = milliers) +
  labs(title = "Une courbe par année (moyenne mobile 7 jours)",
       subtitle = "Même forme chaque année : creux en août et à Noël, pic en septembre et en juin",
       x = NULL, y = "Vélos par jour") + guides(colour = guide_legend(nrow = 1))
g2 <- ggplot(d_normal, aes(mois, velos)) +
  geom_boxplot(fill = "#A1C9F4", outlier.size = 0.4, outlier.alpha = 0.4) +
  scale_y_continuous(labels = milliers) +
  labs(title = "Par mois (hors événements)", x = NULL, y = "Vélos par jour")
g <- g1 / g2 + plot_layout(heights = c(1.3, 1)) + plot_annotation(title = "Saison annuelle : période 365 jours")
sauver(g, "05_saison_annuelle.png", h = 9)

# Carte de chaleur jour de semaine x semaine de l'année, une ligne par année
d[, semaine_annee := pmin(isoweek(date), 52L)]
chaleur <- d[, .(velos = mean(velos)), by = .(annee, semaine_annee, jour_semaine)]
g <- ggplot(chaleur, aes(semaine_annee, jour_semaine, fill = velos)) +
  geom_tile() + facet_wrap(~annee, ncol = 2) +
  scale_fill_viridis_c(option = "magma", labels = milliers, name = "Vélos par jour") +
  scale_y_discrete(limits = rev(jours_fr)) +
  scale_x_continuous(breaks = c(1, 10, 20, 30, 40, 50)) +
  labs(title = "Les deux saisons en même temps : semaine (verticale) et année (horizontale)",
       subtitle = "Bandes sombres : week-ends, août, Noël, confinements. Bandes claires : semaines de rentrée",
       x = "Semaine de l'année", y = NULL) +
  theme(legend.key.width = unit(2, "cm"))
sauver(g, "06_carte_semaine_annee.png", h = 9)

# ===========================================================================
# 5. Tendance : croissance d'une année sur l'autre
# ===========================================================================
annuel <- d[, .(velos = mean(velos), n = .N), by = annee][n >= 300]
annuel[, croissance := velos / shift(velos) - 1]
# Comparaison sur la même période : janvier -> 11 novembre (2026 est incomplète)
meme_periode <- d[jour_annee <= yday(as.IDate("2026-11-11")), .(velos = mean(velos)), by = annee]
meme_periode[, croissance := velos / shift(velos) - 1]
g1 <- ggplot(meme_periode, aes(factor(annee), velos)) +
  geom_col(fill = bleu, width = 0.7) +
  geom_text(aes(label = ifelse(is.na(croissance), "", sprintf("%+.0f %%", 100 * croissance))),
            vjust = -0.4, size = 3.5) +
  scale_y_continuous(labels = milliers, expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Moyenne du 1er janvier au 11 novembre, et évolution sur un an",
       x = NULL, y = "Vélos par jour")
glissement <- d[, .(date, velos_365 = frollmean(velos, 365))]
glissement[, evolution := velos_365 / shift(velos_365, 365) - 1]
g2 <- ggplot(glissement[!is.na(evolution)], aes(as.Date(date), evolution)) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_area(fill = vert, alpha = 0.4) + geom_line(colour = vert) +
  scale_y_continuous(labels = percent) + scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Évolution sur un an de la moyenne des 365 derniers jours",
       subtitle = "Rebond après 2020, puis croissance plus lente mais toujours positive", x = NULL, y = NULL)
g <- g1 / g2 + plot_annotation(title = "Une tendance à la hausse, non linéaire")
sauver(g, "07_tendance_croissance.png", h = 8)

g <- ggplot(d_normal[annee >= 2019], aes(velos, factor(annee), fill = after_stat(x))) +
  ggridges::geom_density_ridges_gradient(scale = 1.8, rel_min_height = 0.01, colour = "white") +
  scale_fill_viridis_c(option = "viridis", guide = "none") +
  scale_x_continuous(labels = milliers) +
  labs(title = "Distribution du trafic journalier, année par année",
       subtitle = "Deux bosses : les week-ends et jours de vacances (à gauche) et les jours ouvrés (à droite)",
       x = "Vélos par jour", y = NULL)
sauver(g, "08_distributions_annuelles.png", h = 6)

# ===========================================================================
# 6. Dépendance temporelle : nuages décalés, ACF, périodogramme
# ===========================================================================
retards <- c(1, 6, 7, 14, 182, 364)
decale <- rbindlist(lapply(retards, function(k)
  data.table(retard = sprintf("Retard %d jour%s", k, ifelse(k > 1, "s", "")),
             x = shift(d$log_velos, k), y = d$log_velos)))
decale[, retard := factor(retard, levels = unique(retard))]
cors <- decale[, .(r = cor(x, y, use = "complete.obs")), by = retard]
g <- ggplot(decale, aes(x, y)) + geom_point(alpha = 0.15, size = 0.5, colour = bleu) +
  geom_abline(colour = rouge) + facet_wrap(~retard, nrow = 2) +
  geom_text(data = cors, aes(x = -Inf, y = Inf, label = sprintf("r = %.2f", r)),
            hjust = -0.2, vjust = 1.5, fontface = "bold") +
  labs(title = "Nuages de points décalés : log(vélos) du jour contre celui d'il y a k jours",
       subtitle = "La corrélation est plus forte à 7 jours qu'à 1 jour : c'est la semaine qui structure la série ; 364 = 52 semaines",
       x = "log(vélos) en t - k", y = "log(vélos) en t")
sauver(g, "09_nuages_decales.png", h = 6.5)

apprentissage <- d[periode == "Apprentissage"]
g1 <- gg_acf(apprentissage$log_velos, 400, "acf", "ACF de log(vélos) jusqu'à 400 jours") +
  scale_x_continuous(breaks = c(0, 7, 91, 182, 273, 365)) +
  labs(subtitle = "Pics tous les 7 jours, décroissance très lente et remontée vers 365 : non stationnaire")
g2 <- gg_acf(apprentissage$log_velos, 35, "acf", "Zoom sur les 35 premiers retards")
g <- g1 / g2 + plot_annotation(caption = "Période d'apprentissage (juin 2021 - novembre 2025). En rouge : multiples de 7")
sauver(g, "10_acf_serie.png", h = 8)

# Périodogramme : quelles fréquences portent la variance ? (estimation de la densité spectrale)
x <- apprentissage$log_velos - fitted(lm(log_velos ~ t, apprentissage))
sp <- spec.pgram(x, taper = 0, detrend = FALSE, plot = FALSE)
spectre <- data.table(frequence = sp$freq, densite = sp$spec)
pics <- data.table(frequence = c(1 / 365.25, 1 / 7, 2 / 7, 3 / 7),
                   nom = c("1 an", "7 jours", "3,5 jours (2/7)", "2,3 jours (3/7)"))
g <- ggplot(spectre, aes(frequence, densite)) +
  geom_line(colour = "grey30", linewidth = 0.3) +
  geom_vline(data = pics, aes(xintercept = frequence), colour = rouge, linetype = 2) +
  geom_label(data = pics, aes(x = frequence, y = max(spectre$densite) * 3, label = nom),
             size = 3, colour = rouge) +
  scale_y_log10() +
  scale_x_sqrt(breaks = c(1 / 365.25, 1 / 30, 1 / 7, 2 / 7, 3 / 7, 0.5),
               labels = c("1/365", "1/30", "1/7", "2/7", "3/7", "1/2")) +
  labs(title = "Périodogramme de log(vélos) sans tendance",
       subtitle = "Les pics à 1/7, 2/7 et 3/7 sont les harmoniques du cycle hebdomadaire ; le pic à gauche est l'année",
       x = "Fréquence (cycles par jour, échelle racine)", y = "Densité spectrale (échelle log)")
sauver(g, "11_periodogramme.png")

# ===========================================================================
# 7. Météo
# ===========================================================================
g1 <- ggplot(d_normal[!est_ferie & !est_week_end], aes(temperature_moyenne, velos)) +
  geom_point(aes(colour = annee), alpha = 0.4, size = 0.8) +
  geom_smooth(method = "loess", formula = y ~ x, colour = rouge) +
  scale_colour_viridis_c(name = NULL, end = 0.9) + scale_y_continuous(labels = milliers) +
  labs(title = "Trafic des jours ouvrés selon la température",
       subtitle = "Mélange l'effet de la température et celui de la saison", x = "Température moyenne (°C)", y = "Vélos par jour") +
  theme(legend.key.width = unit(1.2, "cm"))
g2 <- ggplot(d_normal[!is.na(relatif)], aes(temperature_moyenne, relatif)) +
  geom_point(alpha = 0.15, size = 0.6, colour = bleu) +
  geom_smooth(method = "loess", formula = y ~ x, colour = rouge) +
  geom_hline(yintercept = 1, linetype = 2) +
  coord_cartesian(ylim = c(0.4, 1.4)) +
  labs(title = "Trafic relatif selon la température",
       subtitle = "Saison retirée : effet propre du froid et des fortes chaleurs", x = "Température moyenne (°C)", y = "Trafic relatif")
classes_pluie <- c("Sec (0 mm)", "Bruine (< 1 mm)", "1 à 5 mm", "5 à 10 mm", "Plus de 10 mm")
d_normal[, pluie_classe := factor(fcase(pluie_mm == 0, classes_pluie[1], pluie_mm < 1, classes_pluie[2],
                                        pluie_mm < 5, classes_pluie[3], pluie_mm < 10, classes_pluie[4],
                                        default = classes_pluie[5]), levels = classes_pluie)]
eff_pluie <- d_normal[!is.na(relatif), .(mediane = median(relatif), n = .N), by = pluie_classe]
g3 <- ggplot(d_normal[!is.na(relatif)], aes(pluie_classe, relatif, fill = pluie_classe)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3) +
  geom_text(data = eff_pluie, aes(y = 1.45, label = sprintf("%+.0f %%\n(%d j)", 100 * (mediane - 1), n)), size = 3) +
  scale_fill_brewer(palette = "Blues", guide = "none") +
  coord_cartesian(ylim = c(0.3, 1.5)) +
  labs(title = "Trafic relatif selon la pluie du jour", subtitle = "Effet médian et nombre de jours",
       x = NULL, y = "Trafic relatif")
d_normal[, vent_classe := cut(vent_moyen_m_s, c(0, 2, 3, 4, 5, Inf),
                              labels = c("< 2 m/s", "2-3", "3-4", "4-5", "> 5 m/s"))]
eff_vent <- d_normal[!is.na(relatif), .(mediane = median(relatif), n = .N), by = vent_classe]
g4 <- ggplot(d_normal[!is.na(relatif)], aes(vent_classe, relatif, fill = vent_classe)) +
  geom_hline(yintercept = 1, linetype = 2) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3) +
  geom_text(data = eff_vent, aes(y = 1.45, label = sprintf("%+.0f %%\n(%d j)", 100 * (mediane - 1), n)), size = 3) +
  scale_fill_brewer(palette = "Greens", guide = "none") + coord_cartesian(ylim = c(0.3, 1.5)) +
  labs(title = "Trafic relatif selon le vent", x = "Vent moyen", y = "Trafic relatif")
g <- (g1 | g2) / (g3 | g4) +
  plot_annotation(title = "La météo : la pluie est le facteur le plus fort",
                  caption = "Trafic relatif = trafic du jour / médiane des mêmes jours de semaine ±3 semaines. Hors événements exceptionnels")
sauver(g, "12_meteo.png", l = 12, h = 9)

# Pluie et température ensemble : la pluie pèse-t-elle autant quand il fait chaud ?
d_normal[, temp_classe := cut(temperature_moyenne, c(-Inf, 5, 10, 15, 20, Inf),
                              labels = c("< 5 °C", "5-10 °C", "10-15 °C", "15-20 °C", "> 20 °C"))]
croise <- d_normal[!is.na(relatif), .(relatif = median(relatif), n = .N), by = .(pluie_classe, temp_classe)]
g <- ggplot(croise[n >= 10], aes(temp_classe, pluie_classe, fill = relatif)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%+.0f %%\n%d j", 100 * (relatif - 1), n)), size = 3.2) +
  scale_fill_gradient2(low = rouge, mid = "white", high = bleu, midpoint = 1, labels = percent,
                       name = "Trafic relatif médian") +
  labs(title = "Effet croisé de la température et de la pluie",
       subtitle = "Cases d'au moins 10 jours. Les deux effets s'additionnent (en log) : pas d'interaction marquée",
       x = "Température moyenne", y = "Pluie du jour") +
  theme(legend.key.width = unit(1.5, "cm"))
sauver(g, "13_pluie_temperature.png", h = 5.5)

# ===========================================================================
# 8. Calendrier : jours fériés, ponts, vacances
# ===========================================================================
feries <- d[est_ferie == TRUE & !est_week_end & evenement == "" & !is.na(relatif),
            .(relatif = median(relatif), n = .N), by = nom_jour_ferie]
feries <- rbind(feries, d[pont == TRUE & evenement == "" & !is.na(relatif),
                          .(nom_jour_ferie = "Ponts", relatif = median(relatif), n = .N)])
g1 <- ggplot(feries, aes(reorder(nom_jour_ferie, relatif), relatif - 1)) +
  geom_col(aes(fill = nom_jour_ferie == "Ponts"), width = 0.7) +
  geom_text(aes(label = sprintf("%+.0f %%", 100 * (relatif - 1))), hjust = 1.1, colour = "white", size = 3.2) +
  scale_fill_manual(values = c("FALSE" = rouge, "TRUE" = orange), guide = "none") +
  scale_y_continuous(labels = percent) + coord_flip() +
  labs(title = "Jours fériés tombant en semaine, et ponts",
       subtitle = "Trafic relatif médian : un férié vaut à peu près un dimanche", x = NULL, y = "Écart à un jour comparable")
# Vacances : comparées aux périodes scolaires (on ne peut pas utiliser le trafic
# relatif, ±3 semaines couvrant les vacances elles-mêmes). On utilise le rapport
# à la moyenne des jours ouvrés hors vacances du même mois et de la même année.
ouvres <- d_normal[!est_week_end & !est_ferie]
ref_mois <- ouvres[type_vacances == "Hors vacances", .(ref = mean(velos)), by = .(annee, mois)]
vac <- merge(ouvres, ref_mois, by = c("annee", "mois"))
vac[, rapport := velos / ref]
eff_vac <- vac[type_vacances != "Hors vacances", .(rapport = median(rapport), n = .N), by = nom_vacances]
g2 <- ggplot(eff_vac, aes(reorder(nom_vacances, rapport), rapport - 1)) +
  geom_col(fill = violet, width = 0.7) +
  geom_text(aes(label = sprintf("%+.0f %%", 100 * (rapport - 1))), hjust = 1.1, colour = "white", size = 3.2) +
  scale_y_continuous(labels = percent) + coord_flip() +
  labs(title = "Vacances scolaires (jours ouvrés)",
       subtitle = "Écart aux jours ouvrés hors vacances du même mois. Août est absent : pas de référence",
       x = NULL, y = "Écart")
# Noël en détail : moyenne des années, du 10 décembre au 15 janvier
noel <- d[(month(date) == 12 & mday(date) >= 10) | (month(date) == 1 & mday(date) <= 15)]
noel[, saison := ifelse(month(date) == 12, annee, annee - 1L)]
noel[, jour := as.numeric(date - as.IDate(paste0(saison, "-12-25")))]
noel[, indice := velos / mean(velos[jour < -7]), by = saison]
g3 <- ggplot(noel[saison >= 2019 & saison <= 2025], aes(jour, indice, colour = factor(saison), group = saison)) +
  geom_vline(xintercept = c(0, 7), linetype = 2, colour = "grey50") +
  geom_line(linewidth = 0.7) + scale_colour_viridis_d(name = NULL, end = 0.9) +
  annotate("text", x = c(0, 7), y = 1.3, label = c("Noël", "1er janvier"), size = 3) +
  scale_y_continuous(labels = percent) +
  labs(title = "Autour de Noël (base 100 = 10-17 décembre)",
       subtitle = "Les deux semaines de fêtes divisent le trafic par 2 à 3, chaque année",
       x = "Jours depuis le 25 décembre", y = NULL) + guides(colour = guide_legend(nrow = 1))
g <- (g1 | g2) / g3 + plot_annotation(title = "Le calendrier : fériés, ponts et vacances")
sauver(g, "14_calendrier.png", l = 12, h = 10)

# ===========================================================================
# 9. Événements exceptionnels : observé contre une année « normale »
# ===========================================================================
# Référence : même jour de la semaine, 52 semaines plus tard (ou plus tôt pour les JO).
ev <- rbindlist(lapply(seq_len(nrow(evenements)), function(i) {
  e <- evenements[i]
  fenetre <- d[date >= e$debut - 21L & date <= e$fin + 21L, .(date, velos)]
  decalage <- if (e$evenement == "Jeux olympiques") -364L else 364L
  fenetre[, reference := d$velos[match(date + decalage, d$date)]]
  fenetre[, `:=`(evenement = e$evenement, debut = e$debut, fin = e$fin,
                 annee_reference = if (decalage > 0) "52 semaines après" else "52 semaines avant")]
}))
ev[, evenement := factor(evenement, levels = evenements$evenement)]
effet_ev <- ev[date >= debut & date <= fin, .(effet = mean(velos) / mean(reference) - 1), by = evenement]
g <- ggplot(ev, aes(as.Date(date))) +
  geom_rect(data = unique(ev[, .(evenement, debut, fin)]), inherit.aes = FALSE,
            aes(xmin = as.Date(debut), xmax = as.Date(fin), ymin = -Inf, ymax = Inf), fill = "grey90") +
  geom_line(aes(y = reference, colour = "Même période un an à côté"), linewidth = 0.5) +
  geom_line(aes(y = velos, colour = "Observé"), linewidth = 0.7) +
  geom_text(data = effet_ev, aes(x = as.Date(-Inf), y = Inf, label = sprintf("Écart moyen : %+.0f %%", 100 * effet)),
            hjust = -0.1, vjust = 1.5, size = 3.4, fontface = "bold") +
  facet_wrap(~evenement, scales = "free", ncol = 2) +
  scale_colour_manual(values = c("Observé" = rouge, "Même période un an à côté" = "grey50"), name = NULL) +
  scale_y_continuous(labels = milliers) + scale_x_date(date_labels = "%d %b %y") +
  labs(title = "Les événements exceptionnels, comparés à la même période un an à côté",
       subtitle = "La grève de 2019 a fait monter le vélo, les confinements l'ont effondré ; pendant les JO, Paris a pédalé plus qu'un été normal",
       x = NULL, y = "Vélos par jour")
sauver(g, "15_evenements.png", l = 11, h = 9)

# ===========================================================================
# 10. Synthèse : corrélations entre la série et les variables explicatives
# ===========================================================================
vars <- d_normal[!is.na(relatif), .(
  `log(vélos)` = log_velos, `Trafic relatif` = relatif,
  `Week-end` = as.numeric(est_week_end), `Férié` = as.numeric(est_ferie),
  Vacances = as.numeric(est_vacances), `Température` = temperature_moyenne,
  `Pluie (log)` = log1p(pluie_mm), Vent = vent_moyen_m_s)]
m <- cor(vars)
cm <- as.data.table(as.table(m))
cm[, `:=`(V1 = factor(V1, levels = colnames(m)), V2 = factor(V2, levels = rev(colnames(m))))]
g <- ggplot(cm, aes(V1, V2, fill = N)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", N)), size = 3.3) +
  scale_fill_gradient2(low = rouge, mid = "white", high = bleu, limits = c(-1, 1), name = "Corrélation") +
  labs(title = "Corrélations entre le trafic et les variables explicatives",
       subtitle = "Le week-end pèse sur le niveau ; la pluie et le vent pèsent sur le trafic relatif",
       x = NULL, y = NULL) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.key.width = unit(1.5, "cm"))
sauver(g, "16_correlations.png", l = 8.5, h = 7)

# ===========================================================================
# Résultats chiffrés repris dans le rapport
# ===========================================================================
resume <- data.table(
  indicateur = c("Jours", "Vélos par jour (moyenne)", "Élasticité écart-type / moyenne (brut)",
                 "Lambda de Box-Cox (Guerrero)", "Effet médian d'une pluie > 10 mm",
                 "Effet médian d'un jour férié", "Trafic week-end / semaine (2019)",
                 "Trafic week-end / semaine (2025)"),
  valeur = c(nrow(d), round(mean(d$velos)), round(pente_brut, 2), round(lambda, 2),
             round(100 * (eff_pluie[pluie_classe == "Plus de 10 mm", mediane] - 1)),
             round(100 * (median(d[est_ferie & !est_week_end & evenement == "", relatif], na.rm = TRUE) - 1)),
             round(ratio[year(mois) == 2019, mean(week_end / semaine)], 2),
             round(ratio[year(mois) == 2025, mean(week_end / semaine)], 2)))
ecrire_resultat(resume, "03_exploration_resume.csv")
ecrire_resultat(effet_ev, "03_effet_evenements.csv")
print(resume); print(effet_ev)
