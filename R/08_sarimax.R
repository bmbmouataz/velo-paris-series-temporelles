# 08 - SARIMAX : SARIMA + régresseurs exogènes
#
# Entrée  : data/3_final/velo_paris_par_jour.csv, resultats/modeles/sarima.rds
# Sorties : figures/08_sarimax/*.png, resultats/08_*.csv, resultats/modeles/sarimax.rds
#
# Modèle : régression à erreurs ARMA saisonnières
#   log X_t = a + b t + β' x_t + u_t,   u_t ~ SARIMA(p,0,q)(P,0,Q)_7
# La partie 6 a montré que la différence saisonnière (I - B⁷) est compensée par
# un MA saisonnier de racine unité : la semaine est déterministe. On la met donc
# dans la régression (une indicatrice par jour) et on ne différencie pas (D = 0).
# x_t contient :
#   - le jour de la semaine (6 indicatrices, lundi = référence) ;
#   - la saison annuelle (K paires de Fourier de période 365,25 jours) ;
#   - le calendrier : férié en semaine, pont, vacances (Noël, été, autres), JO ;
#   - la météo du jour : température (et son carré), pluie en log(1 + mm), vent.
# Chaque coefficient se lit comme un effet en % : exp(β) - 1.
#
# On répond ici à la première moitié de la problématique : quelle part du
# trafic s'explique par le calendrier, la météo et les événements ?

source("R/00_outils.R")
preparer_figures("08_sarimax")
d <- charger_serie()
app <- d[periode == "Apprentissage"]
test <- d[periode == "Test"]
y <- ts(app$velos, frequency = 7)
sar <- readRDS("resultats/modeles/sarima.rds")

ajuster <- function(o, s, X) tryCatch(
  Arima(y, order = o, seasonal = s, xreg = X, include.drift = TRUE, lambda = 0, method = "CSS-ML"),
  error = function(e) NULL)

# ===========================================================================
# 1. Choix des ordres ARMA, puis du nombre d'harmoniques K
# ===========================================================================
X4 <- regresseurs(app, K = 4, semaine = TRUE)
grille <- CJ(p = 0:3, q = 0:2, P = 0:1, Q = 0:1)
cat("SARIMAX : estimation de", nrow(grille), "structures d'erreur...\n")
grille[, c("AIC", "BIC") := {
  m <- ajuster(c(p, 0, q), c(P, 0, Q), X4)
  if (is.null(m)) list(NA_real_, NA_real_) else list(m$aic, m$bic)
}, by = .(p, q, P, Q)]
setorder(grille, BIC)
print(head(grille, 5))
best <- grille[1]
ordre <- c(best$p, 0, best$q); saison <- c(best$P, 0, best$Q)

choix_K <- rbindlist(lapply(c(1, 2, 3, 4, 6, 8, 10), function(K) {
  m <- ajuster(ordre, saison, regresseurs(app, K = K, semaine = TRUE))
  data.table(K = K, AIC = m$aic, BIC = m$bic)
}))
K <- choix_K[which.min(BIC), K]
cat("K retenu :", K, "\n")
ecrire_resultat(grille, "08_grille_sarimax.csv")
ecrire_resultat(choix_K, "08_choix_K.csv")

X <- regresseurs(app, K = K, semaine = TRUE)
sarimax <- ajuster(ordre, saison, X)
nom_sarimax <- sprintf("SARIMAX(%d,0,%d)(%d,0,%d)_7 + jours, Fourier K = %d", ordre[1], ordre[3], saison[1], saison[3], K)
cat(nom_sarimax, "\n")

g1 <- ggplot(head(grille[!is.na(BIC)], 10), aes(reorder(sprintf("(%d,0,%d)(%d,0,%d)", p, q, P, Q), -BIC), BIC - min(BIC))) +
  geom_col(fill = bleu, width = 0.6) + coord_flip() +
  labs(title = "Structure des erreurs : 10 meilleures (BIC)", x = NULL, y = "Écart au meilleur BIC")
g2 <- ggplot(melt(choix_K, id.vars = "K"), aes(K, value, colour = variable)) + geom_line() + geom_point() +
  geom_vline(xintercept = K, linetype = 2) +
  scale_colour_manual(values = c(AIC = bleu, BIC = rouge), name = NULL) +
  labs(title = "Nombre d'harmoniques annuelles K", subtitle = sprintf("Le BIC choisit K = %d", K),
       x = "K", y = NULL)
sauver((g1 | g2) + plot_annotation(title = "Sélection du SARIMAX"), "01_selection.png", l = 12, h = 5)

# ===========================================================================
# 2. Coefficients : effets du calendrier et de la météo
# ===========================================================================
coefs <- table_coefficients(sarimax)
ecrire_resultat(coefs, "08_coefficients_sarimax.csv")
print(coefs)
effets <- coefs[parametre %in% names(noms_lisibles)]
effets[, `:=`(libelle = noms_lisibles[parametre], effet = en_pourcent(estimation),
              bas = en_pourcent(estimation - 1.96 * ecart_type), haut = en_pourcent(estimation + 1.96 * ecart_type),
              groupe = fcase(parametre %in% c("temperature", "temperature_2", "pluie", "vent"), "Météo",
                             parametre %in% jours_fr, "Semaine", default = "Calendrier"))]
g <- ggplot(effets, aes(effet, reorder(libelle, effet), colour = groupe)) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_errorbarh(aes(xmin = bas, xmax = haut), height = 0.3, linewidth = 0.8) +
  geom_point(size = 3) +
  geom_text(aes(label = sprintf("%+.1f %%", effet)), vjust = -1, size = 3.2, show.legend = FALSE) +
  facet_grid(groupe ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c("Calendrier" = violet, "Météo" = bleu, "Semaine" = vert), guide = "none") +
  labs(title = "Effet de chaque variable sur le trafic, toutes choses égales par ailleurs",
       subtitle = paste(nom_sarimax, "- intervalles à 95 %. Les p-valeurs sont dans resultats/08_coefficients_sarimax.csv"),
       x = "Effet sur le nombre de vélos (%)", y = NULL)
sauver(g, "02_effets_variables.png", h = 6.5)

# Effet de la température (courbe en cloche) et de la pluie, avec intervalle
cf <- coef(sarimax); V <- sarimax$var.coef
temp <- data.table(temperature = seq(-3, 34, by = 0.5))
temp[, `:=`(x1 = temperature - 15, x2 = (temperature - 15)^2 / 10)]
idx <- c("temperature", "temperature_2")
temp[, effet := x1 * cf[idx[1]] + x2 * cf[idx[2]]]
temp[, se := sqrt(x1^2 * V[idx[1], idx[1]] + x2^2 * V[idx[2], idx[2]] + 2 * x1 * x2 * V[idx[1], idx[2]])]
t_opt <- 15 - cf[idx[1]] / (2 * cf[idx[2]] / 10)
g1 <- ggplot(temp, aes(temperature, en_pourcent(effet) / 100)) +
  geom_ribbon(aes(ymin = en_pourcent(effet - 1.96 * se) / 100, ymax = en_pourcent(effet + 1.96 * se) / 100),
              fill = orange, alpha = 0.25) +
  geom_line(colour = orange, linewidth = 1.2) + geom_hline(yintercept = 0, linetype = 2) +
  {if (t_opt > -3 && t_opt < 34) list(geom_vline(xintercept = t_opt, colour = "grey50", linetype = 3),
                                       annotate("text", x = t_opt, y = 0.03, label = sprintf("optimum %.0f °C", t_opt), size = 3.3))} +
  scale_y_continuous(labels = percent) +
  labs(title = "Effet de la température (référence 15 °C)",
       subtitle = sprintf("Terme au carré non significatif (p = %.2f) : effet presque linéaire", coefs[parametre == "temperature_2", p_valeur]),
       x = "Température moyenne (°C)", y = NULL)
pluie <- data.table(mm = seq(0, 40, by = 0.5))
pluie[, `:=`(effet = log1p(mm) * cf["pluie"], se = log1p(mm) * sqrt(V["pluie", "pluie"]))]
g2 <- ggplot(pluie, aes(mm, en_pourcent(effet) / 100)) +
  geom_ribbon(aes(ymin = en_pourcent(effet - 1.96 * se) / 100, ymax = en_pourcent(effet + 1.96 * se) / 100),
              fill = bleu, alpha = 0.25) +
  geom_line(colour = bleu, linewidth = 1.2) +
  scale_y_continuous(labels = percent) +
  labs(title = "Effet de la pluie (référence : temps sec)", subtitle = "Les premiers millimètres pèsent le plus",
       x = "Pluie du jour (mm)", y = NULL)
# Saison annuelle estimée par Fourier
four <- data.table(jour = 1:365)
four[, effet := as.numeric(fourier_annuel(jour - 1, K) %*% cf[grep("^(sin|cos)_", names(cf))])]
g3 <- ggplot(four, aes(jour, en_pourcent(effet) / 100)) + geom_line(colour = vert, linewidth = 1.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_x_continuous(breaks = cumsum(c(1, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30)), labels = mois_fr) +
  scale_y_continuous(labels = percent) +
  labs(title = sprintf("Saison annuelle (Fourier, K = %d), hors vacances et météo", K), x = NULL, y = NULL)
g <- (g1 | g2) / g3 + plot_annotation(title = "Effets non linéaires et saison annuelle estimés par le SARIMAX")
sauver(g, "03_effets_meteo_saison.png", l = 12, h = 8)

# ===========================================================================
# 3. Chaque bloc de variables est-il utile ? Tests du rapport de vraisemblance
# ===========================================================================
# Modèles emboîtés (mêmes ordres ARMA, D = 0) : on retire un bloc et on compare
#   LR = 2 (log L complet - log L réduit) ~ χ²(nombre de variables retirées).
# Le SARIMA de la partie 6 n'est pas emboîté (D = 1) : on le compare par la
# variance des erreurs de prévision à 1 jour, qui est sur la même échelle.
blocs_test <- list("Météo" = c("temperature", "temperature_2", "pluie", "vent"),
                   "Calendrier" = c("ferie", "pont", "vacances_noel", "vacances_ete", "autres_vacances", "jeux_olympiques"),
                   "Saison annuelle (Fourier)" = grep("^(sin|cos)_", colnames(X), value = TRUE),
                   "Jours de la semaine" = jours_fr[-1])
lr <- rbindlist(lapply(names(blocs_test), function(b) {
  m0 <- ajuster(ordre, saison, X[, setdiff(colnames(X), blocs_test[[b]]), drop = FALSE])
  stat <- 2 * (sarimax$loglik - m0$loglik)
  data.table(bloc_retire = b, ddl = length(blocs_test[[b]]), LR = stat,
             p_valeur = pchisq(stat, length(blocs_test[[b]]), lower.tail = FALSE),
             BIC_sans = m0$bic, BIC_complet = sarimax$bic, ecart_type_sans = sqrt(m0$sigma2))
}))
print(lr)
ecrire_resultat(lr, "08_tests_rapport_vraisemblance.csv")
comparaison <- data.table(modele = c(sar$nom, nom_sarimax),
                          ecart_type_erreur_1_jour = sqrt(c(sar$modele$sigma2, sarimax$sigma2)),
                          AIC = c(sar$modele$aic, sarimax$aic), BIC = c(sar$modele$bic, sarimax$bic))
print(comparaison)
ecrire_resultat(comparaison, "08_sarima_contre_sarimax.csv")

# ===========================================================================
# 4. Diagnostic des résidus
# ===========================================================================
nb_arma <- sum(ordre[c(1, 3)]) + sum(saison[c(1, 3)])
diag_x <- diagnostic_residus(residuals(sarimax), app$date, paste("Diagnostic des résidus du", nom_sarimax),
                  "04_diagnostic_residus.png", fitdf = nb_arma)
print(diag_x)
ecrire_resultat(diag_x, "08_diagnostic_sarimax.csv")

# Erreurs du SARIMA et du SARIMAX côte à côte
err <- data.table(date = app$date, SARIMA = as.numeric(residuals(sar$modele)), SARIMAX = as.numeric(residuals(sarimax)))
err_long <- melt(err, id.vars = "date", variable.name = "modele", value.name = "residu")
g1 <- ggplot(err_long, aes(as.Date(date), en_pourcent(residu) / 100, colour = modele)) +
  geom_line(linewidth = 0.3) + facet_wrap(~modele, ncol = 1) +
  scale_colour_manual(values = c(SARIMA = gris, SARIMAX = bleu), guide = "none") +
  scale_y_continuous(labels = percent) + coord_cartesian(ylim = c(-0.7, 0.7)) +
  labs(title = "Erreurs d'ajustement à un jour", x = NULL, y = NULL)
g2 <- ggplot(err_long, aes(residu, fill = modele)) +
  geom_density(alpha = 0.5, colour = NA) +
  scale_fill_manual(values = c(SARIMA = gris, SARIMAX = bleu), name = NULL) + coord_cartesian(xlim = c(-0.8, 0.6)) +
  labs(title = "Distribution des erreurs", subtitle = sprintf("Écart-type : %.3f -> %.3f",
                                                              sd(err$SARIMA), sd(err$SARIMAX)), x = "Résidu (log)", y = NULL)
g <- (g1 | g2) + plot_layout(widths = c(2, 1)) +
  plot_annotation(title = "Ajouter calendrier et météo réduit fortement les erreurs",
                  subtitle = sprintf("Écart-type de l'erreur à 1 jour (log) : SARIMA %.3f, SARIMAX %.3f. Retirer la météo : LR = %.0f (%d ddl), p = %.1g",
                                     sqrt(sar$modele$sigma2), sqrt(sarimax$sigma2), lr$LR[1], lr$ddl[1], lr$p_valeur[1]))
sauver(g, "05_sarima_contre_sarimax_erreurs.png", l = 12, h = 6)

# Les erreurs restantes : jours que le modèle ne sait pas expliquer
res <- data.table(date = app$date, residu = as.numeric(residuals(sarimax)), app[, .(jour_semaine, nom_jour_ferie, nom_vacances,
                  pluie_mm, temperature_moyenne)])
gros <- res[order(-abs(residu))][1:15]
ecrire_resultat(gros, "08_plus_grosses_erreurs.csv")
g <- ggplot(res, aes(as.Date(date), en_pourcent(residu) / 100)) +
  geom_point(size = 0.6, alpha = 0.5, colour = "grey40") +
  geom_point(data = gros, colour = rouge, size = 2) +
  ggrepel::geom_label_repel(data = gros, aes(label = sprintf("%s %s", jour_semaine, format(date, "%d/%m/%y"))),
                            size = 2.7, max.overlaps = 30) +
  scale_y_continuous(labels = percent) +
  scale_x_date(date_breaks = "6 months", date_labels = "%m/%Y") +
  labs(title = "Les 15 jours que le SARIMAX explique le moins bien",
       subtitle = "Surtout les 24 et 25 décembre (le réveillon n'est pas férié), et des jours exceptionnels comme la journée sans voiture du 19/09/2021",
       x = NULL, y = "Résidu (en % du trafic)")
sauver(g, "06_erreurs_restantes.png", l = 12, h = 5.5)

# Effet ARCH restant ?
r <- as.numeric(residuals(sarimax))
ml <- data.table(h = 1:28, SARIMA = sapply(1:28, function(k) Box.test(as.numeric(residuals(sar$modele))^2, lag = k, type = "Ljung-Box")$p.value),
                 SARIMAX = sapply(1:28, function(k) Box.test(r^2, lag = k, type = "Ljung-Box")$p.value))
ecrire_resultat(ml, "08_mcleod_li.csv")

# ===========================================================================
# 5. Quelle part du trafic s'explique par quoi ? (réponse à la problématique)
# ===========================================================================
# Décomposition de la variance de log(vélos) sur la période d'apprentissage.
# Cinq blocs de variables. Comme la part d'un bloc dépend de l'ordre dans lequel
# on les ajoute, on fait la moyenne sur tous les ordres possibles (valeurs de
# Shapley). Puis la dynamique ARMA explique une partie du reste (R² à 1 jour).
blocs <- list(
  "Tendance"          = cbind(t = app$t),
  "Semaine"           = X[, jours_fr[-1]],
  "Saison annuelle"   = fourier_annuel(app$t, K),
  "Calendrier"        = X[, c("ferie", "pont", "vacances_noel", "vacances_ete", "autres_vacances", "jeux_olympiques")],
  "Météo"             = X[, c("temperature", "temperature_2", "pluie", "vent")]
)
ly <- app$log_velos
r2 <- function(ens) if (length(ens) == 0) 0 else summary(lm(ly ~ do.call(cbind, blocs[ens])))$r.squared
noms_b <- names(blocs)
sous_ensembles <- unlist(lapply(0:length(noms_b), function(k) combn(noms_b, k, simplify = FALSE)), recursive = FALSE)
cle <- function(ens) if (length(ens) == 0) "aucun" else paste(noms_b[noms_b %in% ens], collapse = "+")
r2_ens <- setNames(sapply(sous_ensembles, r2), sapply(sous_ensembles, cle))
shapley <- sapply(noms_b, function(b) {
  autres <- setdiff(noms_b, b); n <- length(noms_b)
  sum(sapply(unlist(lapply(0:length(autres), function(k) combn(autres, k, simplify = FALSE)), recursive = FALSE), function(S) {
    poids <- factorial(length(S)) * factorial(n - length(S) - 1) / factorial(n)
    poids * (r2_ens[[cle(c(S, b))]] - r2_ens[[cle(S)]])
  }))
})
r2_total <- r2_ens[[cle(noms_b)]]
r2_dyn <- 1 - var(residuals(sarimax)) / var(ly)
parts <- data.table(bloc = c(noms_b, "Dynamique ARMA (jours précédents)", "Non expliqué"),
                    part = c(shapley, r2_dyn - r2_total, 1 - r2_dyn))
parts[, bloc := factor(bloc, levels = bloc)]
parts[, fin := cumsum(part)][, debut := fin - part]
print(parts)
ecrire_resultat(parts, "08_part_expliquee.csv")
couleurs_blocs <- c("Tendance" = rouge, "Semaine" = bleu, "Saison annuelle" = vert, "Calendrier" = violet,
                    "Météo" = orange, "Dynamique ARMA (jours précédents)" = "#937860", "Non expliqué" = "grey75")
g <- ggplot(parts, aes(xmin = as.numeric(bloc) - 0.4, xmax = as.numeric(bloc) + 0.4, ymin = debut, ymax = fin, fill = bloc)) +
  geom_rect() +
  geom_text(aes(x = as.numeric(bloc), y = fin, label = percent(part, 0.1)), vjust = -0.5, size = 3.6) +
  scale_x_continuous(breaks = seq_along(parts$bloc), labels = label_wrap(14)(as.character(parts$bloc))) +
  scale_y_continuous(labels = percent, limits = c(0, 1.05)) +
  scale_fill_manual(values = couleurs_blocs, guide = "none") +
  labs(title = "Quelle part des variations du trafic s'explique par quoi ?",
       subtitle = sprintf("Variance de log(vélos), juin 2021 - novembre 2025. Variables : %.0f %% (parts de Shapley) ; avec la dynamique : %.0f %%",
                          100 * r2_total, 100 * r2_dyn),
       x = NULL, y = "Part cumulée de la variance")
sauver(g, "07_part_expliquee.png", l = 11, h = 6)

# ===========================================================================
# 6. Prévision avec la météo observée (borne haute : météo parfaitement prévue)
# ===========================================================================
Xt <- regresseurs(test, K = K, colonnes = colnames(X))
prev <- forecast(sarimax, xreg = Xt, level = c(80, 95), biasadj = TRUE)
pa <- data.table(date = test$date, observe = test$velos, prevision = as.numeric(prev$mean),
                 bas = as.numeric(prev$lower[, 2]), haut = as.numeric(prev$upper[, 2]))
sar_prev <- forecast(sar$modele, h = nrow(test), biasadj = TRUE)
pa[, sarima := as.numeric(sar_prev$mean)]
g <- ggplot(pa, aes(as.Date(date))) +
  geom_ribbon(aes(ymin = bas, ymax = haut), fill = bleu, alpha = 0.15) +
  geom_line(aes(y = observe, colour = "Observé"), linewidth = 0.35) +
  geom_line(aes(y = sarima, colour = "SARIMA"), linewidth = 0.35) +
  geom_line(aes(y = prevision, colour = "SARIMAX (météo observée)"), linewidth = 0.4) +
  scale_colour_manual(values = c("Observé" = "grey35", "SARIMA" = orange, "SARIMAX (météo observée)" = bleu), name = NULL) +
  scale_y_continuous(labels = milliers) + scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
  labs(title = "Toute l'année de test prévue d'un coup, depuis le 11 novembre 2025",
       subtitle = sprintf("Le SARIMAX suit Noël, l'été et les jours de pluie. Erreur moyenne absolue : SARIMA %.0f %%, SARIMAX %.0f %%",
                          100 * mean(abs(pa$sarima / pa$observe - 1)), 100 * mean(abs(pa$prevision / pa$observe - 1))),
       x = NULL, y = "Vélos par jour")
sauver(g, "08_prevision_annee_test.png", l = 12, h = 5.5)

saveRDS(list(modele = sarimax, ordre = ordre, saison = saison, K = K, colonnes = colnames(X), nom = nom_sarimax),
        "resultats/modeles/sarimax.rds")

# ===========================================================================
# 7. Corrélations croisées (CCF) entre la météo et le trafic
# ===========================================================================
# κ(h) = Corr(Z_t, X_{t+h}) : h > 0 = la météo d'aujourd'hui annonce le trafic
# de dans h jours. Pour ne pas confondre avec la saison, on croise :
#   - X : les innovations d'un SARIMAX sans météo (trafic « inexpliqué ») ;
#   - Z : la météo moins sa climatologie (écart à la normale de saison).
# Attention aux corrélations fallacieuses : la météo est elle-même
# autocorrélée (une journée pluvieuse suit souvent une journée pluvieuse).
sans_meteo <- ajuster(ordre, saison, X[, setdiff(colnames(X), blocs_test[["Météo"]])])
u <- as.numeric(residuals(sans_meteo))
clim_app <- climatologie(app, app)
anomalies <- data.table(`Température (écart à la normale)` = app$temperature_moyenne - clim_app$temperature_moyenne,
                        `Pluie, log(1 + mm) (écart à la normale)` = log1p(app$pluie_mm) - log1p(clim_app$pluie_mm),
                        `Vent (écart à la normale)` = app$vent_moyen_m_s - clim_app$vent_moyen_m_s)
ccfs <- rbindlist(lapply(names(anomalies), function(v) {
  cc <- ccf(u, anomalies[[v]], lag.max = 7, plot = FALSE)   # lag h : Corr(u_{t+h}, Z_t)
  data.table(variable = v, h = as.numeric(cc$lag), ccf = as.numeric(cc$acf))
}))
ecrire_resultat(ccfs, "08_ccf_meteo_trafic.csv")
g <- ggplot(ccfs, aes(h, ccf, fill = h == 0)) +
  geom_col(width = 0.6, show.legend = FALSE) + geom_hline(yintercept = 0) +
  geom_hline(yintercept = c(-1.96, 1.96) / sqrt(length(u)), linetype = 2, colour = bleu) +
  facet_wrap(~variable) + scale_fill_manual(values = c("TRUE" = rouge, "FALSE" = gris)) +
  scale_x_continuous(breaks = -7:7) +
  labs(title = "Corrélation croisée entre la météo (Z_t) et le trafic inexpliqué (X_{t+h})",
       subtitle = paste("La pluie agit le jour même (h = 0) : la pluie d'hier n'explique presque rien du trafic d'aujourd'hui.",
                        "
Température : corrélation surtout pour h < 0 (trafic des jours d'avant) : persistance des écarts de température, pas causalité"),
       x = "Décalage h (jours)", y = "CCF")
sauver(g, "09_ccf_meteo_trafic.png", l = 12, h = 5)

# ===========================================================================
# 8. Variante « météo de la veille » : X_t expliqué par Z_{t-1}
# ===========================================================================
# Pour prévoir demain, on connaît la météo d'aujourd'hui : décaler la météo d'un
# jour rend la prévision à 1 jour réalisable sans prévision météo. On estime le
# même SARIMAX avec la météo de la veille, et on compare l'ajustement.
decaler_meteo <- function(dd) {
  out <- copy(dd)
  for (v in c("temperature_moyenne", "pluie_mm", "vent_moyen_m_s")) out[, (v) := shift(get(v), 1L)]
  out
}
d_veille <- decaler_meteo(d)
app_veille <- d_veille[periode == "Apprentissage"]
X_veille <- regresseurs(app_veille, K = K, semaine = TRUE)     # 1er jour : météo du 31 mai 2021
sarimax_veille <- ajuster(ordre, saison, X_veille)
comp_veille <- data.table(modele = c("Météo du jour (Z_t)", "Météo de la veille (Z_{t-1})", "Sans météo"),
                          BIC = c(sarimax$bic, sarimax_veille$bic, sans_meteo$bic),
                          ecart_type_erreur = sqrt(c(sarimax$sigma2, sarimax_veille$sigma2, sans_meteo$sigma2)),
                          effet_pluie = c(eff_p <- en_pourcent(coef(sarimax)["pluie"]), en_pourcent(coef(sarimax_veille)["pluie"]), NA))
print(comp_veille)
ecrire_resultat(comp_veille, "08_meteo_du_jour_ou_de_la_veille.csv")
g <- ggplot(comp_veille, aes(reorder(modele, ecart_type_erreur), ecart_type_erreur, fill = modele)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("σ = %.3f\nBIC = %.0f", ecart_type_erreur, BIC)), vjust = -0.2, size = 3.5) +
  scale_fill_manual(values = c("Météo du jour (Z_t)" = bleu, "Météo de la veille (Z_{t-1})" = "#4C9BB0", "Sans météo" = gris)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
  labs(title = "Météo du jour ou météo de la veille ?",
       subtitle = "Écart-type de l'erreur à 1 jour (log). La météo de la veille n'apporte presque rien : il faut une vraie prévision météo",
       x = NULL, y = "σ")
sauver(g, "10_meteo_du_jour_ou_veille.png", l = 9, h = 5)
saveRDS(list(modele = sarimax_veille, ordre = ordre, saison = saison, K = K, colonnes = colnames(X_veille)),
        "resultats/modeles/sarimax_meteo_veille.rds")
