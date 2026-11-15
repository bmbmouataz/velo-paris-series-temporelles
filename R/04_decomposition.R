# 04 - Décomposition : tendance + saisons + résidu
#
# Entrée  : data/3_final/velo_paris_par_jour.csv
# Sorties : figures/04_decomposition/*.png, resultats/04_*.csv
#
# Modèle multiplicatif X_t = m_t × s_t × z_t, donc additif en log :
#   log X_t = tendance + saison hebdomadaire + saison annuelle + résidu.
# Cinq méthodes, de la plus simple à la plus complète :
#   1. Filtrage linéaire : la moyenne mobile d'ordre 7 annule la semaine.
#   2. decompose() : décomposition classique à une seule saison, 7 jours.
#   3. MSTL : deux saisons (7 et 365,25 jours) estimées par loess.
#   4. Régression paramétrique : tendance + jour de semaine + série de Fourier annuelle.
#   5. Prophet : tendance par morceaux + saisons de Fourier + jours fériés.

source("R/00_outils.R")
preparer_figures("04_decomposition")
d <- charger_serie()
app <- d[periode == "Apprentissage"]

# ===========================================================================
# 1. Filtrage linéaire : moyennes mobiles
# ===========================================================================
# Fonction de transfert de la moyenne arithmétique d'ordre 2k+1 :
#   H(ω) = sin((2k+1)ω/2) / ((2k+1) sin(ω/2)), nulle en ω = 2πj/(2k+1).
# Pour l'ordre 7, elle vaut 0 aux fréquences 1/7, 2/7, 3/7 : elle supprime
# exactement la saison hebdomadaire et conserve les basses fréquences (tendance).
omega <- seq(0.001, pi, length.out = 1000)
gain <- rbindlist(lapply(c(7, 31, 365), function(o)
  data.table(frequence = omega / (2 * pi), gain = abs(sin(o * omega / 2) / (o * sin(omega / 2))),
             filtre = sprintf("Moyenne mobile d'ordre %d", o))))
gain[, filtre := factor(filtre, levels = unique(filtre))]
g1 <- ggplot(gain, aes(frequence, gain, colour = filtre)) +
  geom_vline(xintercept = c(1, 2, 3) / 7, linetype = 2, colour = "grey60") +
  geom_line(linewidth = 0.8) +
  annotate("text", x = c(1, 2, 3) / 7, y = 1.05, label = c("1/7", "2/7", "3/7"), size = 3.2) +
  scale_colour_manual(values = c(bleu, vert, rouge), name = NULL) +
  labs(title = "Gain des moyennes mobiles", subtitle = "L'ordre 7 vaut 0 aux fréquences de la semaine : il l'efface",
       x = "Fréquence (cycles par jour)", y = "|H(ω)|")
d[, `:=`(mm7 = frollmean(log_velos, 7, align = "center"),
         mm365 = frollmean(log_velos, 365, align = "center"))]
zoom <- d[date >= as.IDate("2023-01-01") & date <= as.IDate("2023-12-31")]
g2 <- ggplot(zoom, aes(as.Date(date))) +
  geom_line(aes(y = log_velos, colour = "log(vélos)"), linewidth = 0.3) +
  geom_line(aes(y = mm7, colour = "Moyenne mobile 7 j"), linewidth = 0.8) +
  geom_line(aes(y = mm365, colour = "Moyenne mobile 365 j"), linewidth = 1) +
  scale_colour_manual(values = c("log(vélos)" = "grey60", "Moyenne mobile 7 j" = bleu,
                                 "Moyenne mobile 365 j" = rouge), name = NULL) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  labs(title = "Application à 2023", subtitle = "MM 7 : tendance + saison annuelle ; MM 365 : tendance seule",
       x = NULL, y = "log(vélos)")
g3 <- ggplot(zoom, aes(as.Date(date), log_velos - mm7)) +
  geom_hline(yintercept = 0) + geom_line(colour = violet, linewidth = 0.4) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  labs(title = "log(vélos) - MM 7 = saison hebdomadaire + résidu", x = NULL, y = NULL)
g <- (g1 | g2) / g3 + plot_annotation(title = "Filtrer pour séparer les fréquences")
sauver(g, "01_moyennes_mobiles.png", l = 12, h = 8)

# ===========================================================================
# 2. Décomposition classique à une saison : decompose()
# ===========================================================================
y7 <- ts(app$log_velos, frequency = 7)
dec <- decompose(y7, type = "additive")
compo <- data.table(date = app$date, `Série log(vélos)` = as.numeric(y7),
                    `Tendance (MM 7)` = as.numeric(dec$trend),
                    `Saison hebdomadaire` = as.numeric(dec$seasonal),
                    `Résidu` = as.numeric(dec$random))
compo_long <- melt(compo, id.vars = "date")
g1 <- ggplot(compo_long, aes(as.Date(date), value)) +
  geom_line(colour = bleu, linewidth = 0.3) +
  facet_wrap(~variable, ncol = 1, scales = "free_y", strip.position = "left") +
  labs(title = "decompose() avec une période de 7 jours (apprentissage 2021-2025)",
       subtitle = "La « tendance » garde toute la saison annuelle : une seule saison ne suffit pas",
       x = NULL, y = NULL) + theme(strip.placement = "outside")
coefs_semaine <- data.table(jour = jours_fr[(as.numeric(app$jour_semaine[1:7]) - 1) %% 7 + 1],
                            effet = as.numeric(dec$figure))
coefs_semaine[, jour := factor(jour, levels = jours_fr)]
g2 <- ggplot(coefs_semaine, aes(jour, en_pourcent(effet) / 100, fill = effet > 0)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = sprintf("%+.0f %%", en_pourcent(effet)), vjust = ifelse(effet > 0, -0.4, 1.3)), size = 3.3) +
  scale_fill_manual(values = c("TRUE" = bleu, "FALSE" = orange), guide = "none") +
  scale_y_continuous(labels = percent, expand = expansion(mult = 0.15)) +
  labs(title = "Coefficients saisonniers ŝ_j", subtitle = "En % : exp(ŝ_j) - 1, somme des ŝ_j = 0",
       x = NULL, y = NULL)
g <- g1 + g2 + plot_layout(widths = c(2.2, 1))
sauver(g, "02_decompose_hebdomadaire.png", l = 12, h = 7)

# ===========================================================================
# 3. MSTL : deux saisons à la fois (7 et 365,25 jours)
# ===========================================================================
# Sur toute la série 2019-2026 pour voir les ruptures dans la tendance.
y_multi <- msts(d$log_velos, seasonal.periods = c(7, 365.25))
ms <- mstl(y_multi, s.window = 21, robust = TRUE)
compo_ms <- data.table(date = d$date, as.data.table(as.matrix(ms)))
setnames(compo_ms, c("date", "Série", "Tendance", "Saison 7 jours", "Saison 365 jours", "Résidu"))
long_ms <- melt(compo_ms, id.vars = "date")
g <- ggplot(long_ms, aes(as.Date(date), value)) +
  bandes_evenements(0.25) +
  geom_line(colour = "grey20", linewidth = 0.25) +
  facet_wrap(~variable, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_fill_manual(values = couleurs_evenement, name = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "MSTL : log(vélos) = tendance + saison 7 j + saison 365 j + résidu",
       subtitle = "Loess robuste. Les confinements et la grève se retrouvent dans la tendance et le résidu",
       x = NULL, y = NULL) +
  guides(fill = guide_legend(nrow = 1)) + theme(strip.placement = "outside")
sauver(g, "03_mstl_deux_saisons.png", l = 11, h = 10)

# Part de la variance portée par chaque composante (sur la période d'apprentissage)
dans_app <- compo_ms$date >= debut_apprentissage & compo_ms$date <= fin_apprentissage
variances <- compo_ms[dans_app, lapply(.SD, var), .SDcols = c("Tendance", "Saison 7 jours", "Saison 365 jours", "Résidu")]
variances <- melt(variances, measure.vars = names(variances), variable.name = "composante", value.name = "variance")
variances[, part := variance / sum(variance)]
g1 <- ggplot(variances, aes(reorder(composante, part), part)) +
  geom_col(aes(fill = composante), width = 0.6, show.legend = FALSE) +
  scale_fill_manual(values = c("Tendance" = rouge, "Saison 7 jours" = bleu, "Saison 365 jours" = vert, "Résidu" = gris)) +
  geom_text(aes(label = percent(part, 0.1)), hjust = -0.1) +
  scale_y_continuous(labels = percent, expand = expansion(mult = c(0, 0.15))) + coord_flip() +
  labs(title = "Part de la variance de chaque composante", subtitle = "Apprentissage 2021-2025 (MSTL)",
       x = NULL, y = NULL)
# Le profil hebdomadaire évolue-t-il ? Amplitude de la saison 7 j année par année
compo_ms[, annee := year(date)]
amplitude <- compo_ms[, .(amplitude = diff(range(`Saison 7 jours`))), by = .(mois = as.IDate(format(date, "%Y-%m-01")))]
g2 <- ggplot(amplitude, aes(as.Date(mois), exp(amplitude))) +
  geom_line(colour = bleu) + geom_point(colour = bleu, size = 1) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Amplitude de la saison hebdomadaire",
       subtitle = "Jour le plus fort / jour le plus faible, mois par mois : l'écart semaine / week-end est plus grand en hiver",
       x = NULL, y = "Jour max / jour min")
g <- g1 | g2
sauver(g, "04_mstl_variance_amplitude.png", l = 12, h = 4.8)

# ===========================================================================
# 4. Régression paramétrique : tendance + semaine + Fourier annuel
# ===========================================================================
# log X_t = a + b t + c t² + Σ_j δ_j 1{jour = j} + Σ_k (α_k sin(2πkt/365,25) + β_k cos(...)) + z_t
# Le nombre d'harmoniques K est choisi par AIC / BIC.
criteres <- rbindlist(lapply(1:25, function(K) {
  m <- lm(app$log_velos ~ poly(app$t, 2) + app$jour_semaine + fourier_annuel(app$t, K))
  data.table(K = K, AIC = AIC(m), BIC = BIC(m), R2 = summary(m)$r.squared)
}))
K_bic <- criteres[which.min(BIC), K]; K_aic <- criteres[which.min(AIC), K]
cat(sprintf("Fourier annuel : K retenu par BIC = %d, par AIC = %d\n", K_bic, K_aic))
crit_long <- melt(criteres, id.vars = "K", measure.vars = c("AIC", "BIC"))
g1 <- ggplot(crit_long, aes(K, value, colour = variable)) + geom_line() + geom_point() +
  geom_vline(xintercept = K_bic, linetype = 2) +
  scale_colour_manual(values = c(AIC = bleu, BIC = rouge), name = NULL) +
  scale_x_continuous(breaks = seq(0, 25, 5)) +
  labs(title = "Nombre d'harmoniques annuelles", subtitle = sprintf("Le BIC choisit K = %d : les hautes harmoniques servent à imiter les creux brutaux (Noël, août)", K_bic),
       x = "K (paires sin/cos)", y = NULL)
modele_fourier <- lm(log_velos ~ poly(t, 2, raw = TRUE) + jour_semaine + fourier_annuel(t, K_bic), data = app)
app[, ajuste_fourier := fitted(modele_fourier)]
app[, residu_fourier := residuals(modele_fourier)]
# Saison annuelle estimée (partie Fourier seule) comparée au profil empirique
cf <- coef(modele_fourier); idx <- grep("fourier", names(cf))
profil <- data.table(jour_annee = 1:365)
profil[, estime := as.numeric(fourier_annuel(jour_annee - 1 + as.numeric(as.IDate("2019-01-01") - as.IDate("2019-01-01")), K_bic) %*% cf[idx])]
# le temps t est compté depuis le 1er janvier 2019 : t mod 365,25 ≈ jour de l'année - 1
emp <- app[, .(residu = mean(log_velos - (ajuste_fourier - as.numeric(fourier_annuel(t, K_bic) %*% cf[idx])))),
           by = jour_annee][order(jour_annee)]
g2 <- ggplot() +
  geom_point(data = emp[jour_annee <= 365], aes(jour_annee, residu), colour = "grey60", size = 0.8) +
  geom_line(data = profil, aes(jour_annee, estime), colour = rouge, linewidth = 1.1) +
  scale_x_continuous(breaks = cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30)), labels = mois_fr) +
  labs(title = sprintf("Saison annuelle : série de Fourier, K = %d", K_bic),
       subtitle = "Points : moyenne par jour de l'année des données sans tendance ni semaine",
       x = NULL, y = "Effet sur log(vélos)")
g3 <- ggplot(app, aes(as.Date(date))) +
  geom_line(aes(y = log_velos, colour = "Observé"), linewidth = 0.3) +
  geom_line(aes(y = ajuste_fourier, colour = "Ajusté"), linewidth = 0.4) +
  scale_colour_manual(values = c("Observé" = "grey55", "Ajusté" = rouge), name = NULL) +
  labs(title = sprintf("Ajustement : R² = %.2f", summary(modele_fourier)$r.squared), x = NULL, y = "log(vélos)")
g4 <- gg_acf(app$residu_fourier, 42, "acf", "ACF des résidus de la régression") +
  labs(subtitle = "Encore très autocorrélés : place aux modèles ARMA") + theme(legend.position = "none")
g <- (g1 | g2) / (g3 | g4) + plot_annotation(title = "Décomposition paramétrique par régression")
sauver(g, "05_regression_fourier.png", l = 12, h = 9)
ecrire_resultat(criteres, "04_choix_K_fourier.csv")

# Les résidus de cette régression contiennent la météo et le calendrier
app[, pluie_classe := cut(pluie_mm, c(-Inf, 0, 1, 5, 10, Inf),
                          labels = c("Sec", "< 1 mm", "1-5 mm", "5-10 mm", "> 10 mm"))]
g1 <- ggplot(app, aes(pluie_classe, residu_fourier)) + geom_hline(yintercept = 0) +
  geom_boxplot(fill = "#A1C9F4", outlier.size = 0.4) +
  labs(title = "Résidu selon la pluie", x = NULL, y = "Résidu (log)")
g2 <- ggplot(app, aes(type_vacances, residu_fourier)) + geom_hline(yintercept = 0) +
  geom_boxplot(fill = "#D0BBFF", outlier.size = 0.4) +
  labs(title = "Résidu selon les vacances", x = NULL, y = NULL)
g3 <- ggplot(app, aes(ferie_semaine | pont, residu_fourier)) + geom_hline(yintercept = 0) +
  geom_boxplot(fill = "#FFB482", outlier.size = 0.4) +
  scale_x_discrete(labels = c("Jour normal", "Férié ou pont")) +
  labs(title = "Résidu les jours fériés", x = NULL, y = NULL)
g <- (g1 | g2 | g3) + plot_annotation(
  title = "Ce que la décomposition n'explique pas : météo et calendrier",
  subtitle = "Ces écarts systématiques justifient d'ajouter des régresseurs exogènes (SARIMAX, partie 7)")
sauver(g, "06_residus_contre_exogenes.png", l = 12, h = 4.8)

# ===========================================================================
# 5. Prophet
# ===========================================================================
# y(t) = g(t) tendance linéaire par morceaux + s(t) saisons de Fourier + h(t) fêtes + ε
library(prophet)
fetes <- d[est_ferie == TRUE, .(holiday = "jour férié", ds = as.Date(date), lower_window = 0, upper_window = 0)]
ponts <- d[pont == TRUE, .(holiday = "pont", ds = as.Date(date), lower_window = 0, upper_window = 0)]
m_prophet <- prophet(holidays = rbind(fetes, ponts), yearly.seasonality = 10, weekly.seasonality = TRUE,
                     daily.seasonality = FALSE, changepoint.prior.scale = 0.05)
m_prophet <- add_regressor(m_prophet, "vacances_noel")
m_prophet <- add_regressor(m_prophet, "vacances_ete")
app_p <- app[, .(ds = as.Date(date), y = log_velos,
                 vacances_noel = as.numeric(type_vacances == "Noël"),
                 vacances_ete = as.numeric(type_vacances == "Été"))]
m_prophet <- fit.prophet(m_prophet, app_p)
pr <- as.data.table(predict(m_prophet, app_p))
pr[, ds := as.IDate(ds)]
# Composantes extraites à la main pour des graphiques homogènes
semaine_p <- pr[1:7, .(jour = factor(jours_fr[(wday(ds) + 5) %% 7 + 1], levels = jours_fr), weekly)]
annee_p <- pr[year(ds) == 2023, .(ds, yearly)]
g1 <- ggplot(pr, aes(as.Date(ds), trend)) + geom_line(colour = rouge, linewidth = 1) +
  geom_vline(xintercept = as.Date(m_prophet$changepoints), colour = "grey80", linewidth = 0.3) +
  labs(title = "Tendance par morceaux", subtitle = "Traits gris : points de rupture possibles de la pente",
       x = NULL, y = "log(vélos)")
g2 <- ggplot(semaine_p, aes(jour, weekly)) + geom_col(fill = bleu, width = 0.6) +
  labs(title = "Saison hebdomadaire", x = NULL, y = NULL)
g3 <- ggplot(annee_p, aes(as.Date(ds), yearly)) + geom_line(colour = vert, linewidth = 1) +
  scale_x_date(date_labels = "%b", date_breaks = "1 month") +
  labs(title = "Saison annuelle (Fourier, 10 harmoniques)", x = NULL, y = NULL)
fetes_eff <- pr[, .(`jour férié` = `jour férié`, pont, vacances_noel, vacances_ete)]
fetes_eff <- melt(fetes_eff, measure.vars = names(fetes_eff))[value != 0, .(effet = mean(value)), by = variable]
fetes_eff[, variable := c("jour férié" = "Jour férié", pont = "Pont", vacances_noel = "Vacances de Noël", vacances_ete = "Vacances d'été")[as.character(variable)]]
g4 <- ggplot(fetes_eff, aes(reorder(variable, effet), en_pourcent(effet) / 100)) +
  geom_col(fill = orange, width = 0.6) + coord_flip() +
  geom_text(aes(label = sprintf("%+.0f %%", en_pourcent(effet))), hjust = 1.1, colour = "white") +
  scale_y_continuous(labels = percent) +
  labs(title = "Effets du calendrier", x = NULL, y = NULL)
g5 <- ggplot(pr[app_p, on = "ds"], aes(as.Date(ds))) +
  geom_line(aes(y = y, colour = "Observé"), linewidth = 0.25) +
  geom_line(aes(y = yhat, colour = "Prophet"), linewidth = 0.35) +
  scale_colour_manual(values = c("Observé" = "grey55", "Prophet" = violet), name = NULL) +
  labs(title = sprintf("Ajustement (R² = %.2f)", 1 - var(app_p$y - pr$yhat) / var(app_p$y)), x = NULL, y = "log(vélos)")
g <- (g1 | g2) / (g3 | g4) / g5 + plot_annotation(title = "Prophet : décomposition additive de log(vélos)",
                                                 subtitle = "Apprentissage 2021-2025, jours fériés, ponts et vacances en « fêtes »")
sauver(g, "07_prophet_composantes.png", l = 12, h = 11)

ecrire_resultat(variances, "04_part_variance_mstl.csv")
ecrire_resultat(coefs_semaine[, .(jour, effet_log = effet, effet_pourcent = en_pourcent(effet))],
                "04_coefficients_saisonniers.csv")
print(variances); print(coefs_semaine)

# ===========================================================================
# 6. Méthode par blocs (Buys-Ballot) pour la saison hebdomadaire
# ===========================================================================
# Chaque semaine complète (lundi -> dimanche) est un bloc j. On note X_{j,k} la
# valeur du jour k de la semaine j, m_j la moyenne de la semaine, puis
#   ŝ_k = (1/J) Σ_j (X_{j,k} - m_j)        (somme nulle automatiquement).
blocs <- app[, .(date, log_velos, jour_semaine, semaine = as.IDate(cut(as.Date(date), "week")))]
blocs <- blocs[, if (.N == 7) .SD, by = semaine]
blocs[, ecart := log_velos - mean(log_velos), by = semaine]
s_blocs <- blocs[, .(blocs = mean(ecart)), by = jour_semaine]
s_comp <- merge(s_blocs, coefs_semaine[, .(jour_semaine = jour, decompose = effet)], by = "jour_semaine")
s_long <- melt(s_comp, id.vars = "jour_semaine", variable.name = "methode", value.name = "effet")
s_long[, methode := fifelse(methode == "blocs", "Méthode par blocs", "decompose()")]
g1 <- ggplot(blocs, aes(jour_semaine, ecart)) +
  geom_hline(yintercept = 0, colour = "grey50") +
  geom_boxplot(fill = "#A1C9F4", outlier.size = 0.4, outlier.alpha = 0.4) +
  labs(title = sprintf("Écarts X_{j,k} - m_j des %d semaines complètes", uniqueN(blocs$semaine)),
       subtitle = "Un point par semaine et par jour ; la moyenne de chaque boîte est ŝ_k", x = NULL, y = "Écart (log)")
g2 <- ggplot(s_long, aes(jour_semaine, effet, fill = methode)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("Méthode par blocs" = vert, "decompose()" = bleu), name = NULL) +
  labs(title = "Coefficients saisonniers : deux méthodes",
       subtitle = sprintf("Écart maximal entre les deux : %.3f (en log)", max(abs(s_comp$blocs - s_comp$decompose))),
       x = NULL, y = "ŝ_k")
sauver((g1 | g2) + plot_annotation(title = "Estimer la saison par la méthode des blocs (Buys-Ballot)"),
       "08_methode_par_blocs.png", l = 12, h = 5)

# ===========================================================================
# 7. Filtres complémentaires sur la série hebdomadaire
# ===========================================================================
# Série des moyennes hebdomadaires de log(vélos) : la semaine disparaît, il reste
# tendance + saison annuelle, de période 52 semaines (paire).
#  a) Moyenne mobile d'ordre 5 contre filtre (-3, 12, 17, 12, -3)/35 : les deux ont
#     des poids de somme 1 (conservent les constantes et les droites) ; le second a
#     aussi Σ k² c_k = 0, il conserve les tendances quadratiques et aplatit moins
#     les pics et les creux.
#  b) Moyenne mobile modifiée d'ordre 52 (période paire) :
#     (1/52) [ Σ_{j=-25}^{25} X_{t-j} + (X_{t-26} + X_{t+26}) / 2 ] annule la saison
#     annuelle et ne garde que la tendance.
hebdo <- d[, .(x = mean(log_velos), n = .N), by = .(semaine = as.IDate(cut(as.Date(date), "week")))][n == 7]
filtrer <- function(x, poids) as.numeric(stats::filter(x, poids, sides = 2))
c_mm5 <- rep(1 / 5, 5); c_quad <- c(-3, 12, 17, 12, -3) / 35
c_mod52 <- c(0.5, rep(1, 51), 0.5) / 52
hebdo[, `:=`(mm5 = filtrer(x, c_mm5), quad = filtrer(x, c_quad), mod52 = filtrer(x, c_mod52))]
cat(sprintf("Filtre quadratique : somme %.2f, Σ k c_k = %.2f, Σ k² c_k = %.2f\n",
            sum(c_quad), sum((-2:2) * c_quad), sum((-2:2)^2 * c_quad)))
w <- seq(0.001, pi, length.out = 600)
transfert <- function(cc) { k <- seq_along(cc) - (length(cc) + 1) / 2; sapply(w, function(o) sum(cc * cos(k * o))) }
gains <- rbind(data.table(frequence = w / (2 * pi), H = transfert(c_mm5), filtre = "Moyenne mobile d'ordre 5"),
               data.table(frequence = w / (2 * pi), H = transfert(c_quad), filtre = "Filtre (-3, 12, 17, 12, -3)/35"))
g1 <- ggplot(gains, aes(frequence, H, colour = filtre)) + geom_hline(yintercept = c(0, 1), colour = "grey70") +
  geom_line(linewidth = 1) + scale_colour_manual(values = c("Moyenne mobile d'ordre 5" = bleu, "Filtre (-3, 12, 17, 12, -3)/35" = rouge), name = NULL) +
  labs(title = "Fonctions de transfert H(ω)", subtitle = "Le filtre quadratique reste proche de 1 plus longtemps : il déforme moins",
       x = "Fréquence (cycles par semaine)", y = "H(ω)") + guides(colour = guide_legend(ncol = 1))
z <- hebdo[semaine >= as.IDate("2022-06-01") & semaine <= as.IDate("2024-12-31")]
z_long <- melt(z[, .(semaine, `Série hebdomadaire` = x, `MM 5` = mm5, `Filtre quadratique` = quad)], id.vars = "semaine")
g2 <- ggplot(z_long, aes(as.Date(semaine), value, colour = variable, linewidth = variable)) + geom_line() +
  scale_colour_manual(values = c("Série hebdomadaire" = "grey60", "MM 5" = bleu, "Filtre quadratique" = rouge), name = NULL) +
  scale_linewidth_manual(values = c(0.4, 0.9, 0.9), guide = "none") +
  labs(title = "Lissage de la série hebdomadaire", subtitle = "Le filtre quadratique suit mieux les creux d'août et de Noël",
       x = NULL, y = "log(vélos), moyenne de la semaine")
g3 <- ggplot(hebdo, aes(as.Date(semaine))) +
  bandes_evenements(0.2) +
  geom_line(aes(y = x), colour = "grey60", linewidth = 0.4) +
  geom_line(aes(y = mod52), colour = vert, linewidth = 1.2) +
  scale_fill_manual(values = couleurs_evenement, guide = "none") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Moyenne mobile modifiée d'ordre 52 (période paire)",
       subtitle = "Poids 1/104 aux deux extrémités, 1/52 ailleurs : la saison annuelle est annulée, il reste la tendance",
       x = NULL, y = "log(vélos)")
sauver((g1 | g2) / g3 + plot_annotation(title = "Filtres linéaires sur la série hebdomadaire"),
       "09_filtres_hebdomadaires.png", l = 12, h = 9)
ecrire_resultat(s_comp, "04_saison_blocs_contre_decompose.csv")
