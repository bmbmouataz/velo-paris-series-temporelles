# 10 - Évaluation des prévisions sur l'année de test
#
# Entrées : data/3_final/velo_paris_par_jour.csv, resultats/modeles/*
# Sorties : figures/10_evaluation/*.png, resultats/10_*.csv
#
# Protocole : origine glissante (« validation croisée chronologique »).
#   - Les paramètres sont estimés une seule fois, sur l'apprentissage.
#   - Chaque jour o du test, on suppose connues les données jusqu'à o et on
#     prévoit o+1 ... o+14 (h = 1 à 14 jours). Les modèles dynamiques mettent à
#     jour leur état avec les nouvelles observations, sans ré-estimation.
#   - 352 origines × 14 horizons = 4 928 prévisions par modèle, 14 modèles.
# Mesures : REQM (en vélos), MAPE (en %), couverture des intervalles à 80 et 95 %,
# test de Diebold-Mariano pour comparer deux modèles.
#
# Question de la météo : pour prévoir demain, on ne connaît pas la météo
# de demain. Deux versions du SARIMAX encadrent la réalité :
#   - « météo observée » : comme si la prévision météo était parfaite (borne haute) ;
#   - « météo climatologique » : on ne sait rien, on met la météo moyenne de la saison.
# Une vraie prévision météo à quelques jours se situe entre les deux. Troisième
# version : la météo de la veille (Z_{t-1}), connue au moment de prévoir le lendemain.

source("R/00_outils.R")
source("R/00_previsions.R")
preparer_figures("10_evaluation")
d <- charger_serie()
d <- d[date >= debut_apprentissage]                  # les modèles commencent en juin 2021
app <- d[periode == "Apprentissage"]
n_app <- nrow(app); n <- nrow(d); H <- 14
origines <- n_app:(n - H)
if (Sys.getenv("TEST_RAPIDE") == "1") origines <- origines[c(1, 100, 300)]   # vérification rapide du script
cat("Origines :", length(origines), "\n")

sar <- readRDS("resultats/modeles/sarima.rds")
sx  <- readRDS("resultats/modeles/sarimax.rds")
sx_v <- readRDS("resultats/modeles/sarimax_meteo_veille.rds")   # météo de la veille
sar11 <- readRDS("resultats/modeles/sarima_d1D1.rds")            # (d, D) = (1, 1)
conc <- readRDS("resultats/modeles/concurrents.rds")

# Régresseurs : météo observée, et météo climatologique (moyennes de l'apprentissage)
X_obs  <- regresseurs(d, K = sx$K, colonnes = sx$colonnes)
d_clim <- climatologie(d, app)
X_clim <- regresseurs(d_clim, K = sx$K, colonnes = sx$colonnes)
X_nn   <- regresseurs(d, K = 2, colonnes = conc$colonnes_nnar)
# Météo de la veille : la ligne t contient la météo de t - 1. Pour prévoir o + 1,
# on connaît la météo de o (observée) ; au-delà, la météo de o + h - 1 est inconnue
# et remplacée par la climatologie.
veille <- function(dd) { out <- copy(dd)
  for (v in c("temperature_moyenne", "pluie_mm", "vent_moyen_m_s")) out[, (v) := shift(get(v), 1L)]
  out }
d_full <- charger_serie()
d_veille_obs  <- veille(d_full)[date >= debut_apprentissage]
d_veille_clim <- veille(climatologie(d_full, d_full[periode == "Apprentissage"]))[date >= debut_apprentissage]
X_v_obs  <- regresseurs(d_veille_obs,  K = sx_v$K, colonnes = sx_v$colonnes)
X_v_clim <- regresseurs(d_veille_clim, K = sx_v$K, colonnes = sx_v$colonnes)
y_tout <- d$velos

# Une ligne par (modèle, origine, horizon)
ligne <- function(modele, o, moy, b80 = NA, h80 = NA, b95 = NA, h95 = NA)
  data.table(modele = modele, origine = o, h = seq_len(H), cible = o + seq_len(H),
             prevision = as.numeric(moy), b80 = as.numeric(b80), h80 = as.numeric(h80),
             b95 = as.numeric(b95), h95 = as.numeric(h95))

# ===========================================================================
# 1. Prévisions glissantes
# ===========================================================================
t0 <- Sys.time()
z_tout <- d$log_velos - (conc$tendance[1] + conc$tendance[2] * d$t)
res <- rbindlist(lapply(origines, function(o) {
  y_o <- ts(y_tout[1:o], frequency = 7)
  # Naïf saisonnier : même jour de la semaine précédente
  naif <- y_tout[o + seq_len(H) - 7 * ceiling(seq_len(H) / 7)]
  # SARIMA : mise à jour de l'état avec les données jusqu'à o (paramètres fixes)
  m1 <- Arima(y_o, model = sar$modele)
  f1 <- forecast(m1, h = H, level = c(80, 95), biasadj = TRUE)
  # SARIMAX, météo observée puis climatologique pour les jours à prévoir
  m2 <- Arima(y_o, model = sx$modele, xreg = X_obs[1:o, , drop = FALSE])
  f2 <- forecast(m2, xreg = X_obs[o + seq_len(H), , drop = FALSE], level = c(80, 95), biasadj = TRUE)
  f3 <- forecast(m2, xreg = X_clim[o + seq_len(H), , drop = FALSE], level = c(80, 95), biasadj = TRUE)
  # SARIMAX avec la météo de la veille : seule la prévision à 1 jour utilise une météo observée
  m2v <- Arima(y_o, model = sx_v$modele, xreg = X_v_obs[1:o, , drop = FALSE])
  Xf <- X_v_clim[o + seq_len(H), , drop = FALSE]; Xf[1, ] <- X_v_obs[o + 1, ]
  f2v <- forecast(m2v, xreg = Xf, level = c(80, 95), biasadj = TRUE)
  # SARIMA doublement différencié (d = D = 1)
  m11 <- Arima(y_o, model = sar11$modele)
  f11 <- forecast(m11, h = H, level = c(80, 95), biasadj = TRUE)
  # Lissage exponentiel simple et Holt (paramètres fixés)
  les <- HoltWinters(log(y_o), alpha = unname(conc$lissage_simple[1]), beta = FALSE, gamma = FALSE)
  f_les <- predict(les, n.ahead = H)
  holt <- HoltWinters(log(y_o), alpha = unname(conc$holt[1]), beta = unname(conc$holt[2]), gamma = FALSE)
  f_holt <- predict(holt, n.ahead = H)
  # Holt-Winters, paramètres fixés
  hw <- HoltWinters(log(y_o), alpha = unname(conc$holt_winters[1]), beta = unname(conc$holt_winters[2]),
                    gamma = unname(conc$holt_winters[3]), seasonal = "additive")
  f4 <- predict(hw, n.ahead = H, prediction.interval = TRUE, level = 0.95)
  f4b <- predict(hw, n.ahead = H, prediction.interval = TRUE, level = 0.80)
  # TBATS, paramètres fixés
  m5 <- tbats(msts(log(y_tout[1:o]), seasonal.periods = c(7, 365.25)), model = conc$tbats)
  f5 <- forecast(m5, h = H, level = c(80, 95))
  # Blocs ressemblants
  f6 <- prevision_blocs_noyau(z_tout, o, H = H, b = conc$blocs$b, w = conc$blocs$w)$prevision +
    (conc$tendance[1] + conc$tendance[2] * d$t[o + seq_len(H)])
  # NNAR (sans intervalles : ils demandent des simulations coûteuses)
  m7 <- nnetar(log(y_o), model = conc$nnar, xreg = X_nn[1:o, , drop = FALSE])
  f7 <- forecast(m7, h = H, xreg = X_nn[o + seq_len(H), , drop = FALSE])
  rbind(
    ligne("Naïf saisonnier", o, naif),
    ligne("SARIMA", o, f1$mean, f1$lower[, 1], f1$upper[, 1], f1$lower[, 2], f1$upper[, 2]),
    ligne("SARIMAX (météo observée)", o, f2$mean, f2$lower[, 1], f2$upper[, 1], f2$lower[, 2], f2$upper[, 2]),
    ligne("SARIMAX (météo climatologique)", o, f3$mean, f3$lower[, 1], f3$upper[, 1], f3$lower[, 2], f3$upper[, 2]),
    ligne("SARIMAX (météo de la veille)", o, f2v$mean, f2v$lower[, 1], f2v$upper[, 1], f2v$lower[, 2], f2v$upper[, 2]),
    ligne("SARIMA (d = D = 1)", o, f11$mean, f11$lower[, 1], f11$upper[, 1], f11$lower[, 2], f11$upper[, 2]),
    ligne("Lissage simple", o, exp(f_les)),
    ligne("Holt", o, exp(f_holt)),
    ligne("Holt-Winters", o, exp(f4[, "fit"]), exp(f4b[, "lwr"]), exp(f4b[, "upr"]), exp(f4[, "lwr"]), exp(f4[, "upr"])),
    ligne("TBATS", o, exp(f5$mean), exp(f5$lower[, 1]), exp(f5$upper[, 1]), exp(f5$lower[, 2]), exp(f5$upper[, 2])),
    ligne("Blocs ressemblants", o, exp(f6)),
    ligne("NNAR", o, exp(f7$mean)))
}))
cat("Prévisions glissantes :", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")

# Prophet : ne dépend que de la date et des régresseurs, une seule prédiction par jour
suppressPackageStartupMessages(library(prophet))
tp <- data.table(ds = as.Date(d$date), regresseurs(d, K = 1, colonnes = c("vacances_noel", "vacances_ete",
                 "autres_vacances", "temperature", "temperature_2", "pluie", "vent")))
pp <- as.data.table(predict(conc$prophet, tp))
res_p <- rbindlist(lapply(origines, function(o) {
  k <- o + seq_len(H)
  ligne("Prophet", o, exp(pp$yhat[k]), exp(pp$yhat_lower[k]), exp(pp$yhat_upper[k]))
}))
# LSTM : toutes les origines d'un coup
reticulate::use_virtualenv("r-keras", required = TRUE)
suppressPackageStartupMessages(library(keras3))
lstm <- load_model("resultats/modeles/lstm.keras")
V <- variables_jour(d)
tab <- tableaux_lstm(z_tout / conc$lstm$echelle, V, origines, L = conc$lstm$L, H = H)
pl <- predict(lstm, list(tab$passe, tab$futur), verbose = 0)
res_l <- rbindlist(lapply(seq_along(origines), function(i) {
  o <- origines[i]
  ligne("LSTM", o, exp(pl[i, ] * conc$lstm$echelle + conc$tendance[1] + conc$tendance[2] * d$t[o + seq_len(H)]))
}))
res <- rbind(res, res_p, res_l)
res[, `:=`(date = d$date[cible], observe = y_tout[cible], coupure = d$coupure_generale[cible])]
res <- res[coupure == FALSE]           # jours entièrement estimés : exclus de l'évaluation
res[, erreur := observe - prevision]
fwrite(res, "resultats/10_previsions_glissantes.csv")

ordre_modeles <- c("Naïf saisonnier", "SARIMA", "SARIMA (d = D = 1)", "SARIMAX (météo observée)",
                   "SARIMAX (météo de la veille)", "SARIMAX (météo climatologique)",
                   "Lissage simple", "Holt", "Holt-Winters", "TBATS", "Prophet", "Blocs ressemblants", "NNAR", "LSTM")
couleurs_modeles <- setNames(c(gris, orange, "#FFB482", bleu, "#1F3F7A", "#4C9BB0",
                               "#C9C9C9", "#8C8C8C", "#A1C9F4", vert, violet, rouge, "#DA8BC3", "#937860"),
                             ordre_modeles)
res[, modele := factor(modele, levels = ordre_modeles)]

# ===========================================================================
# 2. Erreurs selon l'horizon
# ===========================================================================
scores <- res[, .(REQM = sqrt(mean(erreur^2)), EAM = mean(abs(erreur)), MAPE = 100 * mean(abs(erreur) / observe),
                  couverture_80 = 100 * mean(observe >= b80 & observe <= h80),
                  couverture_95 = 100 * mean(observe >= b95 & observe <= h95), n = .N), by = .(modele, h)]
ecrire_resultat(scores, "10_scores_par_horizon.csv")
global <- res[, .(REQM = sqrt(mean(erreur^2)), MAPE = 100 * mean(abs(erreur) / observe),
                  MAPE_h1 = 100 * mean((abs(erreur) / observe)[h == 1]),
                  MAPE_h7 = 100 * mean((abs(erreur) / observe)[h == 7]),
                  MAPE_h14 = 100 * mean((abs(erreur) / observe)[h == 14]),
                  couverture_95 = 100 * mean(observe >= b95 & observe <= h95)), by = modele][order(MAPE)]
ecrire_resultat(global, "10_scores_globaux.csv")
print(global)

g1 <- ggplot(scores, aes(h, MAPE, colour = modele)) + geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  scale_colour_manual(values = couleurs_modeles, name = NULL) + scale_x_continuous(breaks = 1:14) +
  labs(title = "MAPE selon l'horizon", x = "Horizon h (jours)", y = "Erreur absolue moyenne (%)")
g2 <- ggplot(scores, aes(h, REQM, colour = modele)) + geom_line(linewidth = 0.9) + geom_point(size = 1.5) +
  scale_colour_manual(values = couleurs_modeles, name = NULL) + scale_x_continuous(breaks = 1:14) +
  scale_y_continuous(labels = milliers) +
  labs(title = "REQM selon l'horizon", x = "Horizon h (jours)", y = "Vélos par jour")
g <- (g1 | g2) + plot_layout(guides = "collect") +
  plot_annotation(title = sprintf("Précision des %d modèles sur l'année de test (origine glissante, %d origines)", uniqueN(res$modele), length(origines)),
                  subtitle = "Plus bas = meilleur") & theme(legend.position = "bottom")
sauver(g, "01_erreur_selon_horizon.png", l = 13, h = 6.5)

tab_h <- scores[h %in% c(1, 2, 3, 7, 14)]
tab_h[, rang := rank(MAPE), by = h]
g <- ggplot(tab_h, aes(factor(h), modele, fill = MAPE)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.1f %%", MAPE), fontface = ifelse(rang == 1, "bold", "plain")), size = 3.5) +
  scale_fill_viridis_c(option = "rocket", direction = -1, name = "MAPE (%)") +
  scale_y_discrete(limits = rev(ordre_modeles)) +
  labs(title = "MAPE par modèle et par horizon", subtitle = "En gras : le meilleur modèle à chaque horizon",
       x = "Horizon (jours)", y = NULL) + theme(panel.grid = element_blank())
sauver(g, "02_tableau_mape.png", l = 9, h = 6)

# ===========================================================================
# 3. Intervalles de prévision : couverture réelle contre couverture annoncée
# ===========================================================================
couv <- melt(scores[!is.na(couverture_95) | !is.na(couverture_80), .(modele, h, `80 %` = couverture_80, `95 %` = couverture_95)],
             id.vars = c("modele", "h"), variable.name = "niveau", value.name = "couverture")[!is.na(couverture)]
g <- ggplot(couv, aes(h, couverture, colour = modele)) +
  geom_hline(data = data.table(niveau = c("80 %", "95 %"), v = c(80, 95)), aes(yintercept = v), linetype = 2) +
  geom_line(linewidth = 0.9) + facet_wrap(~niveau) +
  scale_colour_manual(values = couleurs_modeles, name = NULL) + scale_x_continuous(breaks = c(1, 7, 14)) +
  labs(title = "Les intervalles de prévision tiennent-ils leurs promesses ?",
       subtitle = "Part des jours où l'observé tombe dans l'intervalle. Pointillés : niveau annoncé. En dessous = intervalles trop étroits",
       x = "Horizon (jours)", y = "Couverture observée (%)")
sauver(g, "03_couverture_intervalles.png", l = 11, h = 5.5)

# ===========================================================================
# 4. La valeur de la météo
# ===========================================================================
meteo <- scores[modele %in% c("SARIMA", "SARIMAX (météo observée)", "SARIMAX (météo de la veille)",
                              "SARIMAX (météo climatologique)", "Naïf saisonnier")]
g1 <- ggplot(meteo, aes(h, MAPE, colour = modele)) + geom_line(linewidth = 1) + geom_point() +
  scale_colour_manual(values = couleurs_modeles, name = NULL) + scale_x_continuous(breaks = 1:14) +
  labs(title = "Gain apporté par la connaissance de la météo",
       subtitle = "Météo observée = prévision parfaite\nmétéo de la veille = information sûre à 1 jour ; climatologie = aucune",
       x = "Horizon (jours)", y = "MAPE (%)") + guides(colour = guide_legend(ncol = 2))
pluie_err <- res[modele %in% c("SARIMAX (météo observée)", "SARIMAX (météo de la veille)", "SARIMAX (météo climatologique)") & h == 1]
pluie_err <- merge(pluie_err, d[, .(date, pluie_mm)], by = "date")
pluie_err[, classe := cut(pluie_mm, c(-Inf, 0, 1, 5, Inf), labels = c("Sec", "< 1 mm", "1-5 mm", "> 5 mm"))]
pe <- pluie_err[, .(MAPE = 100 * mean(abs(erreur) / observe)), by = .(modele, classe)]
g2 <- ggplot(pe, aes(classe, MAPE, fill = modele)) + geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = couleurs_modeles, name = NULL) +
  labs(title = "Erreur à 1 jour selon la pluie du jour prévu",
       subtitle = "Sans prévision météo,\nles jours de pluie sont mal prévus", x = NULL, y = "MAPE (%)") +
  guides(fill = guide_legend(ncol = 1))
sauver((g1 | g2) + plot_annotation(title = "Faut-il connaître la météo pour prévoir le trafic ?"),
       "04_valeur_de_la_meteo.png", l = 13, h = 6)

# ===========================================================================
# 5. Prévisions à 1 jour sur toute l'année de test
# ===========================================================================
h1 <- res[h == 1 & modele %in% c("Naïf saisonnier", "SARIMA", "SARIMAX (météo observée)")]
g1 <- ggplot(h1, aes(as.Date(date))) +
  geom_line(data = unique(h1[, .(date, observe)]), aes(y = observe), colour = "grey30", linewidth = 0.4) +
  geom_line(aes(y = prevision, colour = modele), linewidth = 0.35, alpha = 0.85) +
  scale_colour_manual(values = couleurs_modeles, name = NULL) +
  scale_y_continuous(labels = milliers) + scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
  labs(title = "Prévision à 1 jour sur l'année de test", subtitle = "Gris foncé : observé", x = NULL, y = "Vélos par jour")
g2 <- ggplot(h1, aes(as.Date(date), erreur / observe, colour = modele)) +
  geom_hline(yintercept = 0) + geom_line(linewidth = 0.3) +
  facet_wrap(~modele, ncol = 1) + coord_cartesian(ylim = c(-1, 1)) +
  scale_colour_manual(values = couleurs_modeles, guide = "none") + scale_y_continuous(labels = percent) +
  scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
  labs(title = "Erreur relative (observé - prévu) / observé", x = NULL, y = NULL)
sauver(g1 / g2 + plot_layout(heights = c(1, 1.4)), "05_prevision_1_jour_annee_test.png", l = 12, h = 10)

# ===========================================================================
# 6. Où se trompe le meilleur modèle ?
# ===========================================================================
meilleur <- as.character(global$modele[1])
mb <- res[modele == meilleur & h == 1]
mb <- merge(mb, d[, .(date, jour_semaine, mois, type_vacances, est_ferie)], by = "date")
g1 <- ggplot(mb, aes(jour_semaine, abs(erreur) / observe)) + geom_boxplot(fill = "#A1C9F4", outlier.size = 0.5) +
  scale_y_continuous(labels = percent) + coord_cartesian(ylim = c(0, 0.6)) +
  labs(title = "Par jour de la semaine", x = NULL, y = "Erreur absolue relative")
g2 <- ggplot(mb, aes(mois, abs(erreur) / observe)) + geom_boxplot(fill = "#B5D99C", outlier.size = 0.5) +
  scale_y_continuous(labels = percent) + coord_cartesian(ylim = c(0, 0.6)) +
  labs(title = "Par mois", x = NULL, y = NULL)
mb[, ecart := observe / prevision - 1]
pires <- mb[order(-abs(log(observe / prevision)))][1:10]
g3 <- ggplot(pires, aes(reorder(format(date, "%a %d/%m/%y"), abs(log(observe / prevision))), ecart)) +
  geom_col(aes(fill = ecart > 0), width = 0.6, show.legend = FALSE) + coord_flip() +
  geom_text(aes(label = sprintf("%s prévus, %s observés", milliers(prevision), milliers(observe)),
                hjust = ifelse(ecart > 0, -0.05, 1.05)), size = 3) +
  scale_fill_manual(values = c("TRUE" = bleu, "FALSE" = rouge)) +
  scale_y_continuous(labels = percent, expand = expansion(mult = 0.6)) +
  labs(title = "Les 10 pires jours", subtitle = "Écart à la prévision (observé / prévu - 1)", x = NULL, y = NULL)
sauver((g1 | g2) / g3 + plot_annotation(title = paste("Erreurs à 1 jour du meilleur modèle :", meilleur)),
       "06_erreurs_meilleur_modele.png", l = 12, h = 9)
ecrire_resultat(pires[, .(date, observe, prevision, ecart_a_la_prevision = ecart)], "10_pires_jours.csv")

# ===========================================================================
# 7. Test de Diebold-Mariano : les écarts sont-ils significatifs ?
# ===========================================================================
# H0 : les deux modèles ont la même précision (perte quadratique). Pour h > 1,
# les erreurs se chevauchent : dm.test corrige la variance avec h - 1 retards.
dm <- rbindlist(lapply(c(1, 7), function(hh) {
  e <- dcast(res[h == hh], cible ~ modele, value.var = "erreur")
  e <- na.omit(e)
  rbindlist(lapply(setdiff(ordre_modeles, meilleur), function(m) {
    t <- tryCatch(forecast::dm.test(e[[meilleur]], e[[m]], h = hh, power = 2, alternative = "less"),
                  error = function(err) list(statistic = NA_real_, p.value = NA_real_))
    data.table(h = hh, modele = m, statistique = as.numeric(t$statistic), p_valeur = t$p.value)
  }))
}))
ecrire_resultat(dm, "10_diebold_mariano.csv")
g <- ggplot(dm, aes(factor(h), modele, fill = p_valeur < 0.05)) + geom_tile(colour = "white") +
  geom_text(aes(label = format.pval(p_valeur, digits = 2, eps = 1e-4)), size = 3.4) +
  scale_fill_manual(values = c("TRUE" = "#B5D99C", "FALSE" = "#F4A6A6"),
                    labels = c("TRUE" = "significativement meilleur", "FALSE" = "écart non significatif"), name = NULL) +
  scale_y_discrete(limits = rev(setdiff(ordre_modeles, meilleur))) +
  labs(title = paste("Diebold-Mariano :", meilleur, "contre chaque autre modèle"),
       subtitle = "H1 : le meilleur modèle a une erreur quadratique plus faible. p-valeurs à 1 et 7 jours",
       x = "Horizon (jours)", y = NULL) + theme(panel.grid = element_blank())
sauver(g, "07_diebold_mariano.png", l = 9, h = 5.5)

# ===========================================================================
# 8. Prévision réelle : les 14 jours après la fin des données
# ===========================================================================
# Le SARIMAX est ré-estimé sur toutes les données (juin 2021 - 11 novembre 2026)
# avec les mêmes ordres. La météo future est inconnue : on met la climatologie.
futur <- data.table(date = max(d$date) + seq_len(H))
futur[, `:=`(t = as.numeric(date - as.IDate("2019-01-01")), jour_annee = yday(date),
             jour_semaine = factor(jours_fr[(wday(date) + 5) %% 7 + 1], levels = jours_fr))]
feries <- fread("data/2_intermediaire/jours_feries.csv")
vacances <- fread("data/2_intermediaire/vacances_scolaires.csv")
futur[, est_week_end := jour_semaine %in% c("samedi", "dimanche")]
futur[, ferie_semaine := date %in% as.IDate(feries$date) & !est_week_end]
futur[, pont := FALSE]
futur[, nom_vacances := ""]
for (i in seq_len(nrow(vacances))) futur[date >= as.IDate(vacances$premier_jour[i]) & date < as.IDate(vacances$jour_de_rentree[i]),
                                         nom_vacances := vacances$nom_vacances[i]]
futur[, type_vacances := fcase(nom_vacances == "Vacances de Noël", "Noël", nom_vacances == "Vacances d'Été", "Été",
                               nom_vacances != "", "Autres vacances", default = "Hors vacances")]
futur[, evenement := ""]
futur <- climatologie(futur[, `:=`(temperature_moyenne = NA_real_, pluie_mm = NA_real_, vent_moyen_m_s = NA_real_)], d)
X_tout <- regresseurs(d, K = sx$K, colonnes = sx$colonnes)
final <- Arima(ts(d$velos, frequency = 7), order = sx$ordre, seasonal = sx$saison, xreg = X_tout,
               include.drift = TRUE, lambda = 0, method = "CSS-ML")
pf <- forecast(final, xreg = regresseurs(futur, K = sx$K, colonnes = sx$colonnes), level = c(80, 95), biasadj = TRUE)
fut <- data.table(date = futur$date, prevision = as.numeric(pf$mean), b80 = as.numeric(pf$lower[, 1]),
                  h80 = as.numeric(pf$upper[, 1]), b95 = as.numeric(pf$lower[, 2]), h95 = as.numeric(pf$upper[, 2]))
ecrire_resultat(fut, "10_prevision_futur.csv")
g <- ggplot() +
  geom_line(data = d[date > max(date) - 42L], aes(as.Date(date), velos), colour = "grey30") +
  geom_ribbon(data = fut, aes(as.Date(date), ymin = b95, ymax = h95), fill = bleu, alpha = 0.15) +
  geom_ribbon(data = fut, aes(as.Date(date), ymin = b80, ymax = h80), fill = bleu, alpha = 0.25) +
  geom_line(data = fut, aes(as.Date(date), prevision), colour = bleu, linewidth = 1) +
  geom_point(data = fut, aes(as.Date(date), prevision), colour = bleu, size = 1.5) +
  geom_vline(xintercept = as.Date(max(d$date)), linetype = 2) +
  scale_y_continuous(labels = milliers) + scale_x_date(date_breaks = "1 week", date_labels = "%d/%m") +
  labs(title = "Prévision du trafic du 12 au 25 novembre 2026",
       subtitle = paste(sx$nom, "ré-estimé sur toutes les données ; météo future = climatologie. Bandes : 80 % et 95 %"),
       x = NULL, y = "Vélos par jour")
sauver(g, "08_prevision_14_prochains_jours.png", l = 11, h = 5)
