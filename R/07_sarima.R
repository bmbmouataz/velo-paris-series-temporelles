# 07 - Modèle SARIMA
#
# Entrée  : data/3_final/velo_paris_par_jour.csv (période d'apprentissage)
# Sorties : figures/07_sarima/*.png, resultats/07_*.csv, resultats/modeles/sarima.rds
#
# Démarche de Box-Jenkins :
#   1. Identification : ACF / PACF de (I - B⁷) log(vélos) (partie 6 : d = 0, D = 1).
#   2. Estimation d'une grille de modèles SARIMA(p,0,q)(P,1,Q)_7, choix par AIC / BIC.
#   3. Significativité des coefficients (|z| > 1,96), causalité et inversibilité (racines).
#   4. Diagnostic des résidus (Ljung-Box, normalité), effet ARCH.
#   5. Sur-ajustement : ajouter un terme ne doit rien apporter.
#   6. GARCH sur les résidus : la variance change-t-elle dans le temps ?
#   7. Prévision.
#   Compléments : Yule-Walker et PACF, poids de Wold, correction log -> exp,
#   GARCH complet et test BDS, alternative (d, D) = (1, 1), simulation.
#
# Le modèle est estimé sur log(vélos) (lambda = 0) avec une dérive (include.drift) :
# après la différence saisonnière, la tendance linéaire devient une constante.
# Les prévisions sont ramenées en vélos avec la correction de biais exp(m + σ²/2).

source("R/00_outils.R")
preparer_figures("07_sarima")
dir.create("resultats/modeles", showWarnings = FALSE)
d <- charger_serie()
app <- d[periode == "Apprentissage"]
test <- d[periode == "Test"]
y <- ts(app$velos, frequency = 7)
ly <- log(y)

# ===========================================================================
# 1. Identification
# ===========================================================================
dy <- diff(ly, 7)
g1 <- gg_acf(dy, 49, "acf", "ACF de (I - B⁷) log(vélos)") +
  labs(subtitle = "Décroissance rapide (partie AR) et pic négatif au retard 7 (MA saisonnier, Q = 1)")
g2 <- gg_acf(dy, 49, "pacf", "PACF de (I - B⁷) log(vélos)") +
  labs(subtitle = "Pics aux retards 1-2 (AR d'ordre 1 ou 2) puis pics décroissants à 7, 14, 21 (signature d'un MA saisonnier)")
g <- g1 / g2 + plot_annotation(title = "Identification : quels ordres pour le SARIMA ?",
                               subtitle = "Candidats : SARIMA(1 à 3, 0, 0 à 2)(0 ou 1, 1, 1)_7")
sauver(g, "01_identification_acf_pacf.png", h = 8)

# ===========================================================================
# 2. Grille de modèles, choix par AIC et BIC
# ===========================================================================
grille <- CJ(p = 0:4, q = 0:3, P = 0:1, Q = 0:2)
cat("Estimation de", nrow(grille), "modèles...\n")
t0 <- Sys.time()
grille[, c("AIC", "BIC", "sigma2") := {
  m <- tryCatch(Arima(y, order = c(p, 0, q), seasonal = c(P, 1, Q), include.drift = TRUE,
                      lambda = 0, method = "CSS-ML"),
                error = function(e) NULL)
  if (is.null(m)) list(NA_real_, NA_real_, NA_real_) else list(m$aic, m$bic, m$sigma2)
}, by = .(p, q, P, Q)]
cat("Durée :", round(difftime(Sys.time(), t0, units = "mins"), 1), "min\n")
grille[, nb_parametres := p + q + P + Q + 1]
setorder(grille, BIC)
ecrire_resultat(grille, "07_grille_sarima.csv")
print(head(grille, 8))

meilleur <- grille[1]
grille[, saison := sprintf("(P = %d, Q = %d)", P, Q)]
g <- ggplot(grille, aes(factor(q), factor(p), fill = BIC)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(BIC)), size = 2.6) +
  geom_tile(data = grille[1], colour = rouge, fill = NA, linewidth = 1.2) +
  facet_wrap(~saison, nrow = 2) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "BIC") +
  labs(title = "BIC des SARIMA(p, 0, q)(P, 1, Q)_7 sur log(vélos)",
       subtitle = sprintf("Plus c'est foncé, meilleur c'est. Encadré rouge : le minimum, SARIMA(%d,0,%d)(%d,1,%d)_7",
                          meilleur$p, meilleur$q, meilleur$P, meilleur$Q),
       x = "q (ordre MA)", y = "p (ordre AR)") + theme(legend.key.width = unit(1.5, "cm"))
sauver(g, "02_grille_bic.png", l = 12, h = 7)

top <- head(grille, 12)
top[, nom := sprintf("(%d,0,%d)(%d,1,%d)", p, q, P, Q)]
top_long <- melt(top, id.vars = "nom", measure.vars = c("AIC", "BIC"))
top_long[, ecart := value - min(value), by = variable]
g <- ggplot(top_long, aes(reorder(nom, -ecart), ecart, fill = variable)) +
  geom_col(position = "dodge", width = 0.7) + coord_flip() +
  scale_fill_manual(values = c(AIC = bleu, BIC = rouge), name = NULL) +
  labs(title = "Les 12 meilleurs modèles (BIC) : écart au meilleur critère",
       subtitle = "Le BIC pénalise davantage le nombre de paramètres que l'AIC", x = NULL, y = "Écart au minimum")
sauver(g, "03_meilleurs_modeles.png", h = 6)

# Comparaison avec la recherche automatique de forecast
auto <- auto.arima(y, d = 0, D = 1, lambda = 0, stepwise = TRUE, approximation = FALSE, allowdrift = TRUE)
cat("auto.arima propose :", as.character(auto), " BIC =", round(auto$bic, 1), "\n")

# ===========================================================================
# 3. Modèle retenu : significativité, racines
# ===========================================================================
ordre <- c(meilleur$p, 0, meilleur$q); saison <- c(meilleur$P, 1, meilleur$Q)
sarima <- Arima(y, order = ordre, seasonal = saison, include.drift = TRUE, lambda = 0, method = "CSS-ML")
nom_sarima <- sprintf("SARIMA(%d,0,%d)(%d,1,%d)_7", ordre[1], ordre[3], saison[1], saison[3])
coefs <- table_coefficients(sarima)
print(coefs)
ecrire_resultat(coefs, "07_coefficients_sarima.csv")
g1 <- ggplot(coefs, aes(estimation, reorder(parametre, estimation))) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_errorbarh(aes(xmin = estimation - 1.96 * ecart_type, xmax = estimation + 1.96 * ecart_type),
                 height = 0.2, colour = bleu) +
  geom_point(size = 3, colour = bleu) +
  geom_text(aes(label = sprintf("z = %.1f", z)), vjust = -1, size = 3.2) +
  labs(title = "Coefficients estimés et intervalles à 95 %",
       subtitle = "Tous significatifs si l'intervalle ne coupe pas 0 (|z| > 1,96)", x = NULL, y = NULL)
racine_sma <- min(Mod(polyroot(c(1, coef(sarima)[grep("^sma", names(coef(sarima)))]))))
g2 <- gg_racines(sarima, "Causalité et inversibilité") +
  labs(subtitle = sprintf("Inverses des racines. MA saisonnier : |racine| = %.3f, sur le cercle unité
=> (I - B⁷) est en partie compensé : la saison hebdomadaire est surtout déterministe", racine_sma))
g <- (g1 | g2) + plot_annotation(title = nom_sarima)
sauver(g, "04_coefficients_racines.png", l = 12, h = 5.5)

# ===========================================================================
# 4. Diagnostic des résidus
# ===========================================================================
nb_arma <- sum(ordre[c(1, 3)]) + sum(saison[c(1, 3)])
diag_sarima <- diagnostic_residus(residuals(sarima), app$date, paste("Diagnostic des résidus du", nom_sarima),
                       "05_diagnostic_residus.png", fitdf = nb_arma)
print(diag_sarima)
ecrire_resultat(diag_sarima, "07_diagnostic_sarima.csv")

# Les plus gros résidus correspondent-ils à des jours particuliers ?
res <- data.table(date = app$date, residu = as.numeric(residuals(sarima)), app[, .(jour_semaine, nom_jour_ferie,
                  nom_vacances, pluie_mm, temperature_moyenne)])
gros <- res[order(-abs(residu))][1:20]
gros[, explication := fcase(nom_jour_ferie != "", nom_jour_ferie,
                            pluie_mm >= 10, sprintf("pluie %.0f mm", pluie_mm),
                            nom_vacances != "", nom_vacances,
                            temperature_moyenne <= 2, sprintf("froid %.0f °C", temperature_moyenne),
                            default = "?")]
g <- ggplot(gros, aes(reorder(format(date, "%d/%m/%Y"), residu), en_pourcent(residu) / 100,
                      fill = explication != "?")) +
  geom_col(width = 0.7) + coord_flip() +
  geom_text(aes(label = explication, hjust = ifelse(residu < 0, 1.05, -0.05)), size = 3) +
  scale_fill_manual(values = c("TRUE" = bleu, "FALSE" = gris), labels = c("Non expliqué", "Férié, pluie, vacances ou froid"),
                    name = NULL) +
  scale_y_continuous(labels = percent, expand = expansion(mult = 0.35)) +
  labs(title = "Les 20 plus grosses erreurs du SARIMA",
       subtitle = "Presque toutes s'expliquent par le calendrier ou la météo : le SARIMA ne les voit pas venir",
       x = NULL, y = "Erreur d'ajustement (en % du trafic)")
sauver(g, "06_plus_grosses_erreurs.png", h = 7)

# ===========================================================================
# 5. Sur-ajustement : un terme de plus doit être non significatif
# ===========================================================================
sur <- rbindlist(lapply(list(c(1, 0), c(0, 1)), function(ajout) {
  o <- ordre; o[1] <- o[1] + ajout[1]; o[3] <- o[3] + ajout[2]
  m <- tryCatch(Arima(y, order = o, seasonal = saison, include.drift = TRUE, lambda = 0, method = "CSS-ML"),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  cf <- table_coefficients(m)
  terme <- if (ajout[1] == 1) paste0("ar", o[1]) else paste0("ma", o[3])
  data.table(modele = sprintf("SARIMA(%d,0,%d)(%d,1,%d)_7", o[1], o[3], saison[1], saison[3]),
             terme_ajoute = terme, estimation = cf[parametre == terme, estimation],
             p_valeur = cf[parametre == terme, p_valeur], BIC = m$bic, BIC_retenu = sarima$bic)
}))
print(sur)
ecrire_resultat(sur, "07_sur_ajustement.csv")

# ===========================================================================
# 6. Effet ARCH et GARCH sur les résidus
# ===========================================================================
r <- as.numeric(residuals(sarima))
ml <- sapply(1:28, function(k) Box.test(r^2, lag = k, type = "Ljung-Box")$p.value)   # = test de McLeod-Li
suppressPackageStartupMessages(library(rugarch))
spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                   mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                   distribution.model = "std")
garch <- ugarchfit(spec, r, solver = "hybrid")
cg <- as.data.table(garch@fit$matcoef, keep.rownames = "parametre")
setnames(cg, c("parametre", "estimation", "ecart_type", "t", "p_valeur"))
print(cg)
ecrire_resultat(cg, "07_garch_residus.csv")
vol <- data.table(date = app$date, sigma = as.numeric(sigma(garch)), residu = r)
pers <- cg[parametre == "alpha1", estimation] + cg[parametre == "beta1", estimation]
g1 <- ggplot(data.table(h = 1:28, p = ml), aes(h, p)) + geom_point(colour = bleu) +
  geom_hline(yintercept = 0.05, linetype = 2, colour = rouge) + scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Test de McLeod-Li", subtitle = "Ljung-Box sur les résidus² : p < 0,05 = effet ARCH",
       x = "Retard h", y = "p-valeur")
g2 <- ggplot(vol, aes(as.Date(date))) +
  geom_line(aes(y = abs(residu)), colour = "grey70", linewidth = 0.3) +
  geom_line(aes(y = sigma), colour = rouge, linewidth = 0.8) +
  scale_x_date(date_breaks = "6 months", date_labels = "%m/%Y") +
  labs(title = sprintf("Volatilité conditionnelle σ_t du GARCH(1,1) (α + β = %.2f)", pers),
       subtitle = "Gris : |résidu|. La volatilité monte à chaque période de fêtes et de vacances, puis retombe",
       x = NULL, y = NULL)
g <- (g1 | g2) + plot_layout(widths = c(1, 2.5)) +
  plot_annotation(title = "Les résidus du SARIMA sont hétéroscédastiques",
                  subtitle = "ε_t = σ_t η_t, σ²_t = ω + α ε²_{t-1} + β σ²_{t-1}, η de Student")
sauver(g, "07_garch_residus.png", l = 13, h = 5)

# ===========================================================================
# 7. Ajustement et prévision
# ===========================================================================
ajuste <- data.table(date = app$date, observe = app$velos, ajuste = as.numeric(fitted(sarima)))
fin <- ajuste[date > fin_apprentissage - 120L]
prev <- forecast(sarima, h = 28, level = c(80, 95), biasadj = TRUE)
pv <- data.table(date = fin_apprentissage + 1:28, moyenne = as.numeric(prev$mean),
                 b80 = as.numeric(prev$lower[, 1]), h80 = as.numeric(prev$upper[, 1]),
                 b95 = as.numeric(prev$lower[, 2]), h95 = as.numeric(prev$upper[, 2]))
pv <- merge(pv, test[, .(date, observe = velos)], by = "date")
# Valeurs ajustées ± 1,96 σ̂ : bande où devrait tomber 95 % des observations
fin[, `:=`(fb = exp(log(ajuste) - 1.96 * sqrt(sarima$sigma2)), fh = exp(log(ajuste) + 1.96 * sqrt(sarima$sigma2)))]
cat(sprintf("Observations dans fitted ± 1,96 σ̂ (apprentissage) : %.1f %%\n",
            100 * mean(abs(residuals(sarima)) <= 1.96 * sqrt(sarima$sigma2))))
g <- ggplot() +
  geom_ribbon(data = fin, aes(as.Date(date), ymin = fb, ymax = fh), fill = orange, alpha = 0.15) +
  geom_ribbon(data = pv, aes(as.Date(date), ymin = b95, ymax = h95), fill = bleu, alpha = 0.15) +
  geom_ribbon(data = pv, aes(as.Date(date), ymin = b80, ymax = h80), fill = bleu, alpha = 0.25) +
  geom_line(data = fin, aes(as.Date(date), observe, colour = "Observé"), linewidth = 0.5) +
  geom_line(data = fin, aes(as.Date(date), ajuste, colour = "Ajusté (prévision à 1 jour)"), linewidth = 0.5) +
  geom_line(data = pv, aes(as.Date(date), moyenne, colour = "Prévision"), linewidth = 0.9) +
  geom_point(data = pv, aes(as.Date(date), observe, colour = "Observé"), size = 1.3) +
  geom_line(data = pv, aes(as.Date(date), observe, colour = "Observé"), linewidth = 0.5) +
  geom_vline(xintercept = as.Date(fin_apprentissage), linetype = 2) +
  scale_colour_manual(values = c("Observé" = "grey25", "Ajusté (prévision à 1 jour)" = orange, "Prévision" = bleu),
                      name = NULL) +
  scale_y_continuous(labels = milliers) + scale_x_date(date_breaks = "2 weeks", date_labels = "%d/%m/%y") +
  labs(title = paste("Prévision à 28 jours du", nom_sarima, "depuis le 11 novembre 2025"),
       subtitle = "Orange : valeurs ajustées ± 1,96 σ̂. Bleu : intervalles de prévision à 80 % et 95 %. Points : ce qui s'est réellement passé",
       x = NULL, y = "Vélos par jour")
sauver(g, "08_prevision_28_jours.png", l = 12, h = 5.5)

# Sans exogènes, la prévision à long terme retombe sur le profil hebdomadaire moyen
prev_an <- forecast(sarima, h = nrow(test), level = 95, biasadj = TRUE)
pa <- data.table(date = test$date, observe = test$velos, prevision = as.numeric(prev_an$mean),
                 bas = as.numeric(prev_an$lower), haut = as.numeric(prev_an$upper))
g <- ggplot(pa, aes(as.Date(date))) +
  geom_ribbon(aes(ymin = bas, ymax = haut), fill = bleu, alpha = 0.15) +
  geom_line(aes(y = observe, colour = "Observé"), linewidth = 0.3) +
  geom_line(aes(y = prevision, colour = "Prévision SARIMA"), linewidth = 0.4) +
  scale_colour_manual(values = c("Observé" = "grey40", "Prévision SARIMA" = bleu), name = NULL) +
  scale_y_continuous(labels = milliers) + scale_x_date(date_breaks = "1 month", date_labels = "%b %y") +
  labs(title = "Prévision sur toute l'année de test, sans mise à jour",
       subtitle = "Le SARIMA prolonge la tendance et la semaine, mais ignore la saison annuelle, Noël et août : il faut des régresseurs",
       x = NULL, y = "Vélos par jour")
sauver(g, "09_prevision_annee_test.png", l = 12, h = 5)

saveRDS(list(modele = sarima, ordre = ordre, saison = saison, nom = nom_sarima), "resultats/modeles/sarima.rds")
ecrire_resultat(data.table(modele = c(nom_sarima, paste("auto.arima :", as.character(auto))),
                           AIC = c(sarima$aic, auto$aic), BIC = c(sarima$bic, auto$bic),
                           sigma2 = c(sarima$sigma2, auto$sigma2)), "07_modele_retenu.csv")

# ===========================================================================
# 8. PACF et Yule-Walker
# ===========================================================================
# α(h) = φ_hh, dernier coefficient de la projection de X_h sur (X_0, ..., X_{h-1}).
# On le calcule « à la main » en résolvant les équations de Yule-Walker
#   Γ_h (φ_h1, ..., φ_hh)' = (γ(1), ..., γ(h))'
# pour h = 1, 2, ... et on compare à pacf(). Les φ_hj changent avec h (double indice).
x_st <- as.numeric(dy) - mean(dy)
gam <- acf(x_st, lag.max = 30, type = "covariance", plot = FALSE)$acf[, 1, 1]
phi_hh <- sapply(1:28, function(h) {
  G <- toeplitz(gam[1:h])
  solve(G, gam[2:(h + 1)])[h]
})
yw <- data.table(h = 1:28, `Yule-Walker (φ_hh)` = phi_hh, `pacf()` = as.numeric(pacf(x_st, 28, plot = FALSE)$acf))
g1 <- ggplot(yw, aes(h)) +
  geom_col(aes(y = `pacf()`), fill = "grey75", width = 0.6) +
  geom_point(aes(y = `Yule-Walker (φ_hh)`), colour = rouge, size = 2) +
  geom_hline(yintercept = c(-1.96, 1.96) / sqrt(length(x_st)), linetype = 2, colour = bleu) +
  scale_x_continuous(breaks = seq(0, 28, 7)) +
  labs(title = "PACF = dernier coefficient de Yule-Walker",
       subtitle = sprintf("Barres : pacf() ; points : φ_hh par Γ_h φ = γ. Écart maximal %.1e", max(abs(phi_hh - yw$`pacf()`))),
       x = "h", y = "α(h)")
# Estimation d'un AR par Yule-Walker contre maximum de vraisemblance
ar_yw <- ar(x_st, method = "yule-walker", order.max = 14, aic = TRUE)
ar_mle <- ar(x_st, method = "mle", order.max = ar_yw$order, aic = FALSE)
comp_ar <- data.table(coefficient = paste0("φ", seq_len(ar_yw$order)),
                      `Yule-Walker` = ar_yw$ar, `Maximum de vraisemblance` = ar_mle$ar)
g2 <- ggplot(melt(comp_ar, id.vars = "coefficient"), aes(factor(coefficient, levels = comp_ar$coefficient), value, fill = variable)) +
  geom_col(position = "dodge", width = 0.7) + geom_hline(yintercept = 0) +
  scale_fill_manual(values = c(vert, violet), name = NULL) +
  labs(title = sprintf("AR(%d) sur (I - B⁷) y : deux estimateurs", ar_yw$order),
       subtitle = "Ordre choisi par AIC ; estimations très proches (n grand)", x = NULL, y = NULL)
sauver((g1 | g2) + plot_annotation(title = "Autocorrélations partielles et équations de Yule-Walker"),
       "10_yule_walker_pacf.png", l = 12, h = 5)

# ===========================================================================
# 9. Décomposition de Wold et intervalle de prévision
# ===========================================================================
# En développant tous les polynômes (y compris la différence I - B⁷), le SARIMA
# s'écrit X_t = Σ ψ_k ε_{t-k}. L'erreur de prévision à h pas vaut Σ_{k<h} ψ_k ε_{t+h-k},
# de variance σ² Σ_{k<h} ψ_k² : c'est elle qui donne la largeur des intervalles.
produit <- function(a, b) { r <- rep(0, length(a) + length(b) - 1)
  for (i in seq_along(a)) r[i:(i + length(b) - 1)] <- r[i:(i + length(b) - 1)] + a[i] * b; r }
mod <- sarima$model
poly_ar <- produit(c(1, -mod$phi), c(1, -mod$Delta))     # Φ(B) Φ_s(B⁷) (I - B⁷)
psi <- c(1, ARMAtoMA(ar = -poly_ar[-1], ma = mod$theta, lag.max = 59))
pf <- forecast(sarima, h = 28, level = 95)
demi_largeur_forecast <- (log(pf$upper[, 1]) - log(pf$lower[, 1])) / 2
demi_largeur_main <- 1.96 * sqrt(sarima$sigma2 * cumsum(psi[1:28]^2))
g1 <- ggplot(data.table(k = 0:59, psi = psi), aes(k, psi)) +
  geom_segment(aes(xend = k, yend = 0), colour = bleu) + geom_point(colour = bleu, size = 1) +
  geom_hline(yintercept = 0) + scale_x_continuous(breaks = seq(0, 56, 7)) +
  labs(title = "Poids ψ_k de la forme MA(∞)", subtitle = "Ils ne s'annulent pas : la différence saisonnière garde la mémoire des chocs",
       x = "k", y = "ψ_k")
g2 <- ggplot(data.table(h = 1:28, main = demi_largeur_main, forecast = as.numeric(demi_largeur_forecast)), aes(h)) +
  geom_line(aes(y = main, colour = "1,96 σ √(Σ ψ_k²), calculé à la main"), linewidth = 1) +
  geom_point(aes(y = forecast, colour = "forecast()"), size = 2) +
  scale_colour_manual(values = c(rouge, "grey20"), name = NULL) +
  labs(title = "Demi-largeur de l'intervalle à 95 % (échelle log)",
       subtitle = "Les deux calculs coïncident : l'intervalle s'élargit avec h", x = "Horizon h", y = NULL) +
  guides(colour = guide_legend(ncol = 1))
sauver((g1 | g2) + plot_annotation(title = "Décomposition de Wold du SARIMA et intervalles de prévision"),
       "11_wold_intervalles.png", l = 12, h = 5)

# ===========================================================================
# 10. Retour du log à l'échelle des vélos : correction de biais
# ===========================================================================
# Si log Z = X et l'erreur est N(0, σ²), E[Z] = exp(X̂ + σ²/2) : exp(X̂) sous-estime.
# On compare sur les prévisions à 1 jour de l'année de test :
#   naïf exp(X̂) ; correction log-normale exp(X̂ + σ̂²/2) ; facteur empirique (1/n) Σ exp(ε̃_t).
d_mod <- d[date >= debut_apprentissage]
m_tout <- Arima(ts(d_mod$velos, frequency = 7), model = sarima)
e_tout <- as.numeric(residuals(m_tout, type = "innovation"))       # sur l'échelle log
dans_test <- d_mod$date >= debut_test
obs <- d_mod$velos[dans_test]; naif_exp <- obs / exp(e_tout[dans_test])
facteurs <- c("exp(X̂) sans correction" = 1, "exp(X̂ + σ̂²/2)" = exp(sarima$sigma2 / 2),
              "Facteur empirique (1/n) Σ exp(ε̃)" = mean(exp(residuals(sarima))))
biais <- data.table(methode = factor(names(facteurs), levels = names(facteurs)), facteur = facteurs,
                    biais = sapply(facteurs, function(fc) mean(obs - naif_exp * fc) / mean(obs)))
print(biais)
ecrire_resultat(biais, "07_correction_log_exp.csv")
g <- ggplot(biais, aes(methode, biais, fill = methode)) + geom_col(width = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0) +
  geom_text(aes(label = sprintf("%+.1f %%\n(facteur %.3f)", 100 * biais, facteur), vjust = ifelse(biais > 0, -0.3, 1.2)), size = 3.5) +
  scale_fill_manual(values = c(gris, bleu, vert)) +
  scale_y_continuous(labels = percent, expand = expansion(mult = 0.4)) +
  labs(title = "Biais moyen des prévisions à 1 jour selon la correction du retour à l'échelle",
       subtitle = "Année de test. Positif : la prévision sous-estime le trafic. Les bornes des intervalles, elles, se transforment directement par exp",
       x = NULL, y = "(observé - prévu) / observé")
sauver(g, "12_correction_log_exp.png", l = 10, h = 5)

# ===========================================================================
# 11. GARCH complet et bruit blanc faible ou fort
# ===========================================================================
# Bruit blanc fort = i.i.d. ; faible = seulement non corrélé. Les résidus passent
# Ljung-Box (non corrélés) mais leurs carrés sont corrélés : ils ne sont PAS
# indépendants. Le test BDS (H0 : i.i.d.) le confirme.
bds <- bds.test(r, m = 3)
bds_dt <- data.table(dimension = rep(2:3, each = ncol(bds$p.value)), epsilon = rep(round(bds$parameter$eps, 3), 2),
                     statistique = as.numeric(t(bds$statistic)), p_valeur = as.numeric(t(bds$p.value)))
print(bds_dt)
ecrire_resultat(bds_dt, "07_test_bds.csv")

# Comparaison ARCH(1), ARCH(2), GARCH(1,1) avec tseries::garch (order = c(β, α))
rc <- r - mean(r)
gfits <- list("ARCH(1)" = c(0, 1), "ARCH(2)" = c(0, 2), "GARCH(1,1)" = c(1, 1))
comp_garch <- rbindlist(lapply(names(gfits), function(nm) {
  g <- suppressWarnings(garch(rc, order = gfits[[nm]], trace = FALSE))
  ll <- as.numeric(logLik(g)); k <- length(coef(g))
  cf <- coef(g); somme <- sum(cf[-1])
  data.table(modele = nm, log_vraisemblance = ll, BIC = -2 * ll + k * log(length(rc)),
             omega = cf[1], somme_alpha_beta = somme,
             sigma_stationnaire = if (somme < 1) sqrt(cf[1] / (1 - somme)) else NA_real_)
}))
print(comp_garch)
ecrire_resultat(comp_garch, "07_comparaison_garch.csv")
g11 <- suppressWarnings(garch(rc, order = c(1, 1), trace = FALSE))
sig_t <- fitted(g11)[, 1]
sig_stat <- comp_garch[modele == "GARCH(1,1)", sigma_stationnaire]
vol <- data.table(date = app$date, r = rc, sigma = as.numeric(sig_t))
g1 <- ggplot(vol, aes(as.Date(date))) +
  geom_line(aes(y = r), colour = "grey65", linewidth = 0.25) +
  geom_line(aes(y = sigma, colour = "± σ_t conditionnel"), linewidth = 0.5) +
  geom_line(aes(y = -sigma, colour = "± σ_t conditionnel"), linewidth = 0.5) +
  geom_hline(aes(yintercept = sig_stat, colour = "± σ stationnaire = √(ω / (1 - α - β))"), linetype = 2) +
  geom_hline(aes(yintercept = -sig_stat, colour = "± σ stationnaire = √(ω / (1 - α - β))"), linetype = 2) +
  scale_colour_manual(values = c(rouge, bleu), name = NULL) +
  scale_x_date(date_breaks = "6 months", date_labels = "%m/%Y") +
  labs(title = "Résidus du SARIMA et volatilité du GARCH(1,1) (tseries::garch)", x = NULL, y = NULL)
# Prévision de la volatilité avec rugarch : elle revient vers le niveau stationnaire
pv_garch <- ugarchforecast(garch, n.ahead = 60)
vol_prev <- data.table(h = 1:60, sigma = as.numeric(sigma(pv_garch)))
sig_stat_rug <- sqrt(cg[parametre == "omega", estimation] / (1 - pers))
g2 <- ggplot(vol_prev, aes(h, sigma)) + geom_line(colour = rouge, linewidth = 1) +
  geom_hline(yintercept = sig_stat_rug, linetype = 2, colour = bleu) +
  labs(title = "Prévision de la volatilité σ_{n+h} (rugarch)",
       subtitle = "On prévoit la volatilité, pas le signe du résidu (prévision 0)", x = "Horizon h (jours)", y = "σ prévu")
crit <- ggplot(comp_garch, aes(modele, BIC - min(BIC))) + geom_col(fill = violet, width = 0.6) +
  geom_text(aes(label = sprintf("%.0f", BIC - min(BIC))), vjust = -0.3) +
  labs(title = "BIC (écart au meilleur)", x = NULL, y = NULL)
bds_g <- ggplot(bds_dt, aes(factor(epsilon), statistique, fill = factor(dimension))) +
  geom_col(position = "dodge", width = 0.7) + geom_hline(yintercept = 1.96, linetype = 2) +
  scale_fill_manual(values = c(bleu, vert), name = "Dimension") +
  labs(title = "Test BDS (H0 : i.i.d.)", subtitle = "Statistique > 1,96 : rejet. Bruit blanc faible, pas fort",
       x = "Seuil ε", y = "Statistique")
sauver(g1 / (crit | g2 | bds_g) + plot_annotation(title = "Hétéroscédasticité conditionnelle des résidus du SARIMA"),
       "13_garch_complet.png", l = 13, h = 9)

# ===========================================================================
# 12. Alternative (d, D) = (1, 1)
# ===========================================================================
# Différencier aussi en (I - B) : pas de dérive possible (d + D = 2). Les tests de
# la partie 6 disaient d = 0 ; on vérifie ce que donne quand même ce choix.
grille_11 <- CJ(p = 0:2, q = 0:2, P = 0:1, Q = 1)
grille_11[, BIC := {
  m <- tryCatch(Arima(y, order = c(p, 1, q), seasonal = c(P, 1, Q), lambda = 0, method = "CSS-ML"), error = function(e) NULL)
  if (is.null(m)) NA_real_ else m$bic
}, by = .(p, q, P, Q)]
b11 <- grille_11[which.min(BIC)]
sarima_11 <- Arima(y, order = c(b11$p, 1, b11$q), seasonal = c(b11$P, 1, b11$Q), lambda = 0, method = "CSS-ML")
nom_11 <- sprintf("SARIMA(%d,1,%d)(%d,1,%d)_7", b11$p, b11$q, b11$P, b11$Q)
coefs_11 <- table_coefficients(sarima_11); print(nom_11); print(coefs_11)
ecrire_resultat(coefs_11, "07_coefficients_sarima_d1D1.csv")
ma_11 <- coef(sarima_11)[grep("^ma", names(coef(sarima_11)))]
racine_ma_11 <- if (length(ma_11)) min(Mod(polyroot(c(1, ma_11)))) else NA
g1 <- gg_racines(sarima_11, paste("Racines du", nom_11)) +
  labs(subtitle = sprintf("Plus petite racine du MA non saisonnier : |z| = %.3f\n(proche de 1 => la différence I - B est de trop)", racine_ma_11))
diagnostic_residus(residuals(sarima_11), app$date, paste("Diagnostic des résidus du", nom_11), "14_diagnostic_sarima_d1D1.png",
        fitdf = b11$p + b11$q + b11$P + b11$Q)
sauver(g1, "15_racines_sarima_d1D1.png", l = 7, h = 6)
saveRDS(list(modele = sarima_11, nom = nom_11), "resultats/modeles/sarima_d1D1.rds")
ecrire_resultat(data.table(modele = nom_11, BIC = sarima_11$bic, sigma = sqrt(sarima_11$sigma2), racine_ma_min = racine_ma_11),
                "07_sarima_d1D1.csv")

# ===========================================================================
# 13. Simuler le modèle estimé
# ===========================================================================
# Si le modèle est bon, des trajectoires simulées doivent « ressembler » à la série :
# même niveau, même semaine, même ACF. Ce qui manque saute aux yeux.
# On simule 2 ans de trajectoires qui démarrent à la fin de l'apprentissage
# (future = TRUE : le modèle part des dernières valeurs observées), et on les
# compare aux 2 dernières années réellement observées.
set.seed(2026)
n_sim <- 730
sims <- rbindlist(lapply(1:3, function(i)
  data.table(jour = 1:n_sim, valeur = as.numeric(simulate(sarima, nsim = n_sim, future = TRUE)),
             trajectoire = paste("Simulation", i))))
sims <- rbind(data.table(jour = 1:n_sim, valeur = tail(app$velos, n_sim), trajectoire = "Série observée (nov. 2023 - nov. 2025)"),
              sims)
sims[, trajectoire := factor(trajectoire, levels = unique(trajectoire))]
g1 <- ggplot(sims, aes(jour, valeur)) +
  geom_line(aes(colour = grepl("observée", trajectoire)), linewidth = 0.3, show.legend = FALSE) +
  facet_wrap(~trajectoire, ncol = 1) + scale_colour_manual(values = c("TRUE" = "grey20", "FALSE" = bleu)) +
  scale_y_continuous(labels = milliers) +
  labs(title = "Deux ans observés contre trois trajectoires simulées",
       subtitle = "Départ le 11 novembre 2025.\nLes simulations dérivent plus que la réalité : la différence saisonnière accumule les chocs",
       x = "Jour", y = "Vélos par jour")
acfs <- sims[, .(retard = 1:28, acf = as.numeric(acf(diff(log(valeur), 7), 28, plot = FALSE)$acf[-1])), by = trajectoire]
g2 <- ggplot(acfs, aes(retard, acf, colour = grepl("observée", trajectoire), group = trajectoire)) +
  geom_line(linewidth = 0.8) + geom_hline(yintercept = 0) +
  scale_colour_manual(values = c("TRUE" = "grey20", "FALSE" = bleu), labels = c("Simulations", "Observée"), name = NULL) +
  scale_x_continuous(breaks = seq(0, 28, 7)) +
  labs(title = "ACF de (I - B⁷) log", subtitle = "Même dépendance à court terme ;\nmais pas de creux d'août ni de Noël dans les simulations",
       x = "Retard", y = NULL)
sauver((g1 | g2) + plot_layout(widths = c(1.6, 1)) + plot_annotation(title = "Simulation à partir du modèle estimé"),
       "16_simulations.png", l = 13, h = 8)
