# 09 - Modèles concurrents
#
# Entrée  : data/3_final/velo_paris_par_jour.csv
# Sorties : figures/09_concurrents/*.png, resultats/modeles/*.rds, resultats/modeles/lstm.keras
#
# Les SARIMA / SARIMAX sont mis en concurrence avec six autres approches :
#   1. Naïf saisonnier : demain = même jour la semaine dernière (référence minimale).
#   2. Lissage de Holt-Winters (pondération par l'éloignement).
#   3. TBATS : lissage exponentiel à deux saisons (7 et 365,25 jours).
#   4. Prophet avec calendrier et météo.
#   5. Ressemblance entre blocs, noyau de Nadaraya-Watson calibré par validation croisée.
#   6. Réseaux de neurones : NNAR (perceptron sur les retards) et LSTM.
# Ici on estime chaque modèle sur l'apprentissage et on illustre son
# fonctionnement ; la comparaison chiffrée est faite dans R/10_evaluation.R.

source("R/00_outils.R")
source("R/00_previsions.R")
preparer_figures("09_concurrents")
d <- charger_serie()
app <- d[periode == "Apprentissage"]
test <- d[periode == "Test"]
n_app <- nrow(app)
ly <- ts(app$log_velos, frequency = 7)
H <- 28
fenetre_test <- test[1:H]
modeles <- list()

# ===========================================================================
# 1. Holt-Winters additif sur log(vélos), saison 7
# ===========================================================================
# c_t = a (X_t - s_{t-7}) + (1 - a)(c_{t-1} + l_{t-1})   niveau
# l_t = b (c_t - c_{t-1}) + (1 - b) l_{t-1}               pente
# s_t = g (X_t - c_t) + (1 - g) s_{t-7}                   saison
# a, b, g minimisent l'erreur de prévision à un pas.
# Progression : lissage simple (niveau), Holt (niveau + pente),
# Holt-Winters (niveau + pente + saison). Paramètres optimisés sur les erreurs à 1 pas.
les  <- HoltWinters(ly, beta = FALSE, gamma = FALSE)
holt <- HoltWinters(ly, gamma = FALSE)
hw <- HoltWinters(ly, seasonal = "additive")
modeles$lissage_simple <- c(alpha = unname(les$alpha))
modeles$holt <- c(alpha = unname(holt$alpha), beta = unname(holt$beta))
trois <- list("Lissage simple" = les, "Holt" = holt, "Holt-Winters" = hw)
prev3 <- rbindlist(lapply(names(trois), function(nm) {
  p <- predict(trois[[nm]], n.ahead = H)
  data.table(modele = nm, date = fenetre_test$date, prevision = exp(as.numeric(p)))
}))
prev3[, modele := factor(modele, levels = names(trois))]
params3 <- data.table(modele = factor(names(trois), levels = names(trois)),
                      parametres = c(sprintf("a = %.2f", les$alpha), sprintf("a = %.2f, b = %.3f", holt$alpha, holt$beta),
                                     sprintf("a = %.2f, b = %.3f, c = %.2f", hw$alpha, hw$beta, hw$gamma)),
                      SSE = c(les$SSE, holt$SSE, hw$SSE))
g3 <- ggplot() +
  geom_line(data = app[date > fin_apprentissage - 42L], aes(as.Date(date), velos), colour = "grey45") +
  geom_point(data = fenetre_test, aes(as.Date(date), velos), size = 1, colour = "grey30") +
  geom_line(data = prev3, aes(as.Date(date), prevision, colour = modele), linewidth = 1) +
  geom_text(data = params3, aes(x = as.Date(fin_apprentissage - 40L), y = 85000,
                                label = sprintf("%s\nSSE (log) = %.1f", parametres, SSE)), hjust = 0, size = 3.2) +
  facet_wrap(~modele) + scale_colour_manual(values = c(gris, orange, bleu), guide = "none") +
  scale_y_continuous(labels = milliers, limits = c(0, 92000)) + scale_x_date(date_labels = "%d/%m") +
  labs(title = "Trois lissages exponentiels de complexité croissante",
       subtitle = "Le lissage simple prolonge un niveau plat, Holt une droite ; seul Holt-Winters reproduit la semaine",
       x = NULL, y = "Vélos par jour")
sauver(g3, "07_lissages_simple_holt_hw.png", l = 13, h = 5)

cat(sprintf("Holt-Winters : alpha = %.3f, beta = %.3f, gamma = %.3f\n", hw$alpha, hw$beta, hw$gamma))
modeles$holt_winters <- c(alpha = hw$alpha, beta = hw$beta, gamma = hw$gamma)
comp <- data.table(date = app$date[-(1:7)], as.data.table(hw$fitted))
comp_long <- melt(comp[, .(date, `Niveau c_t` = level, `Pente l_t` = trend, `Saison s_t` = season)], id.vars = "date")
g1 <- ggplot(comp_long, aes(as.Date(date), value)) + geom_line(colour = bleu, linewidth = 0.3) +
  facet_wrap(~variable, ncol = 1, scales = "free_y") +
  labs(title = sprintf("Composantes mises à jour chaque jour (α = %.2f, β = %.3f, γ = %.2f)", hw$alpha, hw$beta, hw$gamma),
       subtitle = "α élevé : le niveau suit de près les derniers jours ; β ≈ 0 : pente presque fixe", x = NULL, y = NULL)
p_hw <- predict(hw, n.ahead = H, prediction.interval = TRUE, level = 0.95)
g2 <- ggplot() +
  geom_line(data = app[date > fin_apprentissage - 56L], aes(as.Date(date), velos), colour = "grey40") +
  geom_ribbon(data = data.table(date = fenetre_test$date, b = exp(p_hw[, "lwr"]), h = exp(p_hw[, "upr"])),
              aes(as.Date(date), ymin = b, ymax = h), fill = bleu, alpha = 0.2) +
  geom_line(data = data.table(date = fenetre_test$date, p = exp(p_hw[, "fit"])), aes(as.Date(date), p), colour = bleu, linewidth = 1) +
  geom_point(data = fenetre_test, aes(as.Date(date), velos), size = 1.2) +
  scale_y_continuous(labels = milliers) +
  labs(title = "Prévision à 28 jours", subtitle = "Points : observé. Niveau prolongé, sans saison annuelle ni météo", x = NULL, y = "Vélos par jour")
sauver((g1 | g2) + plot_annotation(title = "Lissage exponentiel de Holt-Winters (saison 7 jours, sur log)"),
       "01_holt_winters.png", l = 13, h = 6)

# ===========================================================================
# 2. TBATS : deux saisons trigonométriques + erreurs ARMA
# ===========================================================================
t0 <- Sys.time()
tb <- tbats(msts(app$log_velos, seasonal.periods = c(7, 365.25)), use.box.cox = FALSE,
            use.parallel = FALSE, use.arma.errors = TRUE)
cat("TBATS :", as.character(tb), "en", round(difftime(Sys.time(), t0, units = "secs")), "s\n")
modeles$tbats <- tb
comp_tb <- as.data.table(tbats.components(tb))
comp_tb[, date := app$date]
tb_long <- melt(comp_tb[, .(date, Niveau = level, `Saison 7 jours` = season1, `Saison 365 jours` = season2)], id.vars = "date")
g <- ggplot(tb_long, aes(as.Date(date), value)) + geom_line(colour = vert, linewidth = 0.3) +
  facet_wrap(~variable, ncol = 1, scales = "free_y") +
  labs(title = paste("TBATS :", as.character(tb)),
       subtitle = "Les deux saisons sont des sommes de sinusoïdes dont les coefficients évoluent dans le temps",
       x = NULL, y = NULL)
sauver(g, "02_tbats_composantes.png", h = 7)

# ===========================================================================
# 3. Prophet avec calendrier et météo
# ===========================================================================
suppressPackageStartupMessages(library(prophet))
vars_p <- c("vacances_noel", "vacances_ete", "autres_vacances", "temperature", "temperature_2", "pluie", "vent")
table_prophet <- function(dd) {
  X <- regresseurs(dd, K = 1, colonnes = c("vacances_noel", "vacances_ete", "autres_vacances",
                                           "temperature", "temperature_2", "pluie", "vent"))
  data.table(ds = as.Date(dd$date), y = dd$log_velos, X)
}
fetes <- d[(est_ferie == TRUE) | pont == TRUE,
           .(holiday = ifelse(pont, "pont", "jour férié"), ds = as.Date(date), lower_window = 0, upper_window = 0)]
mp <- prophet(holidays = fetes, yearly.seasonality = 10, weekly.seasonality = TRUE, daily.seasonality = FALSE)
for (v in vars_p) mp <- add_regressor(mp, v)
mp <- fit.prophet(mp, table_prophet(app))
modeles$prophet <- mp
pp <- as.data.table(predict(mp, table_prophet(fenetre_test)))
# Validation croisée propre à Prophet : estimations successives sur une
# fenêtre qui s'agrandit (au moins 2 ans), une nouvelle origine tous les 60 jours.
cv <- suppressMessages(cross_validation(mp, initial = 730, period = 60, horizon = 14, units = "days"))
perf <- as.data.table(performance_metrics(cv, rolling_window = 0))
perf[, h := as.numeric(horizon, units = "days")]
ecrire_resultat(perf[, .(h, mape, rmse, coverage)], "09_prophet_validation_croisee.csv")
g_cv <- ggplot(perf, aes(h, 100 * mape)) + geom_line(colour = violet, linewidth = 1) + geom_point(colour = violet) +
  scale_x_continuous(breaks = 1:14) +
  labs(title = "Validation croisée de Prophet (cross_validation)",
       subtitle = sprintf("%d origines dans l'apprentissage. Erreur mesurée sur log(vélos) : MAPE du log",
                          uniqueN(cv$cutoff)),
       x = "Horizon (jours)", y = "MAPE du log (%)")
sauver(g_cv, "08_prophet_validation_croisee.png", l = 9, h = 4.5)
eff <- data.table(variable = c("Température (effet du jour)", "Pluie", "Vent", "Vacances", "Jours fériés et ponts"),
                  ecart_type = c(sd(pp$temperature + pp$temperature_2), sd(pp$pluie), sd(pp$vent),
                                 sd(pp$vacances_noel + pp$vacances_ete + pp$autres_vacances), sd(pp$holidays)))
g1 <- ggplot(pp, aes(as.Date(ds))) +
  geom_ribbon(aes(ymin = exp(yhat_lower), ymax = exp(yhat_upper)), fill = violet, alpha = 0.2) +
  geom_line(aes(y = exp(yhat)), colour = violet, linewidth = 1) +
  geom_point(data = fenetre_test, aes(as.Date(date), velos), size = 1.2) +
  scale_y_continuous(labels = milliers) +
  labs(title = "Prophet : prévision à 28 jours (météo observée)",
       subtitle = "Prophet ne regarde pas les jours précédents : seulement la date et les régresseurs", x = NULL, y = "Vélos par jour")
comp_p <- melt(pp[, .(ds, Tendance = trend - mean(trend), Semaine = weekly, `Saison annuelle` = yearly,
                      `Fêtes` = holidays, `Météo` = temperature + temperature_2 + pluie + vent,
                      Vacances = vacances_noel + vacances_ete + autres_vacances)], id.vars = "ds")
g2 <- ggplot(comp_p, aes(as.Date(ds), value, colour = variable)) + geom_line(linewidth = 0.9) +
  geom_hline(yintercept = 0, colour = "grey50") +
  scale_colour_manual(values = c(rouge, bleu, vert, orange, "#4C9BB0", violet), name = NULL) +
  labs(title = "Ce qui compose la prévision (échelle log)", x = NULL, y = NULL) +
  guides(colour = guide_legend(nrow = 2))
sauver((g1 | g2) + plot_annotation(title = "Prophet avec jours fériés, vacances et météo"), "03_prophet.png", l = 13, h = 5.5)

# ===========================================================================
# 4. Ressemblance entre blocs
# ===========================================================================
tendance <- lm(log_velos ~ t, app)
modeles$tendance <- coef(tendance)
z_tout <- d$log_velos - predict(tendance, d)
o <- which(d$date == fin_apprentissage)
# Calibration de (b, w) par validation croisée chronologique : sur la
# dernière année d'apprentissage, une origine tous les 3 jours, h = 1 à 14.
z_app <- z_tout[seq_len(n_app)]
orig_val <- seq(n_app - 365L, n_app - 14L, by = 3L)
calib <- CJ(b = c(7L, 14L, 21L, 28L), w = c(0.03, 0.05, 0.08, 0.12, 0.2, 0.35))
calib[, MAPE := {
  e <- unlist(lapply(orig_val, function(oo) {
    p <- prevision_blocs_noyau(z_app, oo, H = 14, b = b, w = w)$prevision
    abs(1 - exp(p - z_app[oo + 1:14]))
  }))
  100 * mean(e)
}, by = .(b, w)]
best_b <- calib[which.min(MAPE)]
cat(sprintf("Blocs (noyau) : b = %d jours, w = %.2f, MAPE de validation %.1f %%\n", best_b$b, best_b$w, best_b$MAPE))
modeles$blocs <- list(b = best_b$b, w = best_b$w)
ecrire_resultat(calib, "09_calibration_blocs.csv")
g_cal <- ggplot(calib, aes(factor(w), factor(b), fill = MAPE)) + geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.1f", MAPE)), size = 3.3) +
  geom_tile(data = best_b, colour = rouge, fill = NA, linewidth = 1.2) +
  scale_fill_viridis_c(option = "mako", direction = -1, name = "MAPE (%)") +
  labs(title = "Calibration par validation croisée", subtitle = "Dernière année d'apprentissage, h = 1 à 14",
       x = "Fenêtre w du noyau", y = "Taille b du bloc (jours)")
pb <- prevision_blocs_noyau(z_tout, o, H = H, b = best_b$b, w = best_b$w)
top <- order(-pb$poids)[1:10]
pb$voisins <- pb$fins[top]; pb$poids <- pb$poids[top] / sum(pb$poids[top])
L <- best_b$b
vois <- rbindlist(lapply(seq_along(pb$voisins), function(i) {
  j <- pb$voisins[i]
  data.table(voisin = format(d$date[j], "%d/%m/%Y"), poids = pb$poids[i], jour = (-L + 1):H,
             z = z_tout[(j - L + 1):(j + H)] - mean(z_tout[(j - L + 1):j]))
}))
courant <- data.table(jour = (-L + 1):H, z = z_tout[(o - L + 1):(o + H)] - mean(z_tout[(o - L + 1):o]))
g1 <- ggplot() +
  geom_vline(xintercept = 0.5, linetype = 2) +
  geom_line(data = vois, aes(jour, z, group = voisin, alpha = poids), colour = vert, linewidth = 0.6) +
  geom_line(data = courant, aes(jour, z), colour = "grey20", linewidth = 1) +
  geom_line(data = data.table(jour = 1:H, z = pb$prevision - mean(z_tout[(o - L + 1):o])), aes(jour, z),
            colour = rouge, linewidth = 1.2) +
  annotate("text", x = c(-L / 2, 14), y = 0.9, label = c("Bloc de référence", "Suite à prévoir"), size = 3.4) +
  scale_alpha(range = c(0.25, 0.9), guide = "none") +
  labs(title = "Blocs les plus ressemblants et leur suite",
       subtitle = sprintf("Noir : %d derniers jours puis réalité.
Vert : voisins (%s...). Rouge : moyenne pondérée",
                          L, paste(head(unique(vois$voisin), 3), collapse = ", ")),
       x = "Jours depuis le 11 novembre 2025", y = "log(vélos) sans tendance, centré")
pb_dt <- data.table(date = fenetre_test$date, prevision = exp(pb$prevision + predict(tendance, fenetre_test)))
g2 <- ggplot() +
  geom_line(data = app[date > fin_apprentissage - 56L], aes(as.Date(date), velos), colour = "grey40") +
  geom_line(data = pb_dt, aes(as.Date(date), prevision), colour = rouge, linewidth = 1) +
  geom_point(data = fenetre_test, aes(as.Date(date), velos), size = 1.2) +
  scale_y_continuous(labels = milliers) +
  labs(title = "Prévision à 28 jours", subtitle = "Moyenne des suites des voisins", x = NULL, y = "Vélos par jour")
sauver((g_cal | g1 | g2) + plot_layout(widths = c(0.9, 1.4, 1)) +
         plot_annotation(title = sprintf("Ressemblance entre blocs, noyau gaussien (Nadaraya-Watson), b = %d jours, w = %.2f",
                                         best_b$b, best_b$w),
                         subtitle = "Au milieu : les 10 blocs passés de plus fort poids (tous les blocs passés sont utilisés)"),
       "04_ressemblance_blocs.png", l = 15, h = 5.5)

# ===========================================================================
# 5. NNAR : perceptron à une couche cachée sur les valeurs passées
# ===========================================================================
# Entrées : y_{t-1..t-p}, y_{t-7}, y_{t-14} et les régresseurs du jour ; une
# couche cachée ; 20 réseaux moyennés. Le réseau est non linéaire mais sans mémoire.
set.seed(2026)
X_nn <- regresseurs(app, K = 2)
nn <- nnetar(ly, P = 2, xreg = X_nn, repeats = 20, size = 10, MaxNWts = 5000)
cat("NNAR :", nn$method, "\n")
modeles$nnar <- nn
modeles$colonnes_nnar <- colnames(X_nn)
pn <- forecast(nn, h = H, xreg = regresseurs(fenetre_test, K = 2, colonnes = colnames(X_nn)))
nn_dt <- data.table(date = fenetre_test$date, prevision = exp(as.numeric(pn$mean)))

# ===========================================================================
# 6. LSTM
# ===========================================================================
# Réseau récurrent : la cellule c_t garde une mémoire, filtrée par trois portes
# (oubli f_t, entrée i_t, sortie o_t). Architecture :
#   séquence des 28 derniers jours -> LSTM(32) ─┐
#   variables des 14 jours à prévoir -> Dense(32) ┴-> Dense(32) -> 14 sorties
reticulate::use_virtualenv("r-keras", required = TRUE)
suppressPackageStartupMessages(library(keras3))
set_random_seed(2026)
echelle <- sd(residuals(tendance))
z_std <- z_tout / echelle
V <- variables_jour(d)
H_lstm <- 14; L_lstm <- 28
orig_app <- (L_lstm):(n_app - H_lstm)
tab <- tableaux_lstm(z_std, V, orig_app, L = L_lstm, H = H_lstm)
entree_passe <- keras_input(shape = c(L_lstm, ncol(V) + 1))
entree_futur <- keras_input(shape = H_lstm * ncol(V))
branche1 <- entree_passe |> layer_lstm(units = 32, dropout = 0.2)
branche2 <- entree_futur |> layer_dense(units = 32, activation = "relu")
sortie <- layer_concatenate(list(branche1, branche2)) |>
  layer_dense(units = 32, activation = "relu") |> layer_dense(units = H_lstm)
lstm <- keras_model(list(entree_passe, entree_futur), sortie)
lstm |> compile(optimizer = optimizer_adam(learning_rate = 5e-4), loss = "mse")
histo <- lstm |> fit(list(tab$passe, tab$futur), tab$cible, epochs = 300, batch_size = 64,
                     validation_split = 0.15, verbose = 0,
                     callbacks = list(callback_early_stopping(patience = 30, restore_best_weights = TRUE)))
save_model(lstm, "resultats/modeles/lstm.keras", overwrite = TRUE)
modeles$lstm <- list(echelle = echelle, L = L_lstm, H = H_lstm)
pertes <- data.table(epoque = seq_along(histo$metrics$loss), Apprentissage = histo$metrics$loss,
                     Validation = histo$metrics$val_loss)
cat("LSTM :", nrow(pertes), "époques, perte de validation minimale", round(min(pertes$Validation), 3), "\n")
g1 <- ggplot(melt(pertes, id.vars = "epoque"), aes(epoque, value, colour = variable)) + geom_line(linewidth = 0.8) +
  scale_colour_manual(values = c(Apprentissage = bleu, Validation = rouge), name = NULL) +
  labs(title = "Apprentissage du LSTM", subtitle = "Erreur quadratique (données standardisées). Arrêt quand la validation ne progresse plus",
       x = "Époque", y = "MSE")
tab_o <- tableaux_lstm(z_std, V, o, L = L_lstm, H = H_lstm)
pl <- predict(lstm, list(tab_o$passe, tab_o$futur), verbose = 0)
lstm_dt <- data.table(date = d$date[o + 1:H_lstm],
                      prevision = exp(as.numeric(pl) * echelle + predict(tendance, d[o + 1:H_lstm])))
g2 <- ggplot() +
  geom_line(data = app[date > fin_apprentissage - 56L], aes(as.Date(date), velos), colour = "grey40") +
  geom_line(data = lstm_dt, aes(as.Date(date), prevision), colour = "#937860", linewidth = 1) +
  geom_point(data = fenetre_test[1:H_lstm], aes(as.Date(date), velos), size = 1.2) +
  scale_y_continuous(labels = milliers) +
  labs(title = "Prévision à 14 jours", x = NULL, y = "Vélos par jour")
sauver((g1 | g2) + plot_annotation(title = "Réseau LSTM : 28 jours passés + calendrier et météo des 14 jours suivants"),
       "05_lstm.png", l = 12, h = 5)

# ===========================================================================
# 7. Tous les modèles, même départ : 28 premiers jours de test
# ===========================================================================
sar <- readRDS("resultats/modeles/sarima.rds"); sx <- readRDS("resultats/modeles/sarimax.rds")
p_sar <- forecast(sar$modele, h = H, biasadj = TRUE)
p_sx <- forecast(sx$modele, xreg = regresseurs(fenetre_test, K = sx$K, colonnes = sx$colonnes), biasadj = TRUE)
p_tb <- forecast(tb, h = H)
naif <- app$velos[n_app - 7 + ((0:(H - 1)) %% 7) + 1]
tous <- rbind(
  data.table(modele = "Naïf saisonnier", date = fenetre_test$date, prevision = as.numeric(naif)),
  data.table(modele = "SARIMA", date = fenetre_test$date, prevision = as.numeric(p_sar$mean)),
  data.table(modele = "SARIMAX", date = fenetre_test$date, prevision = as.numeric(p_sx$mean)),
  data.table(modele = "Holt-Winters", date = fenetre_test$date, prevision = as.numeric(exp(p_hw[, "fit"]))),
  data.table(modele = "TBATS", date = fenetre_test$date, prevision = exp(as.numeric(p_tb$mean))),
  data.table(modele = "Prophet", date = fenetre_test$date, prevision = as.numeric(exp(pp$yhat))),
  data.table(modele = "Blocs ressemblants", pb_dt),
  data.table(modele = "NNAR", nn_dt),
  data.table(modele = "LSTM", lstm_dt))
tous <- merge(tous, fenetre_test[, .(date, observe = velos)], by = "date")
tous[, modele := factor(modele, levels = unique(modele))]
erreurs <- tous[, .(mape = mean(abs(prevision / observe - 1))), by = modele]
g <- ggplot(tous, aes(as.Date(date))) +
  geom_line(aes(y = observe), colour = "grey40") + geom_point(aes(y = observe), colour = "grey40", size = 0.8) +
  geom_line(aes(y = prevision, colour = modele), linewidth = 0.9) +
  geom_text(data = erreurs, aes(x = as.Date(fin_apprentissage + 2L), y = 82000, label = sprintf("MAPE %.0f %%", 100 * mape)),
            hjust = 0, size = 3.3) +
  facet_wrap(~modele, ncol = 3) +
  scale_colour_manual(values = c(gris, orange, bleu, "#4C9BB0", vert, violet, rouge, "#DA8BC3", "#937860"), guide = "none") +
  scale_y_continuous(labels = milliers, limits = c(0, 90000)) + scale_x_date(date_labels = "%d/%m") +
  labs(title = "Neuf modèles, un même point de départ : prévoir les 28 jours après le 11 novembre 2025",
       subtitle = "Gris : observé. Une seule origine, donc une illustration : la comparaison rigoureuse est dans R/10_evaluation.R",
       x = NULL, y = "Vélos par jour")
sauver(g, "06_comparaison_28_jours.png", l = 12, h = 9)

saveRDS(modeles, "resultats/modeles/concurrents.rds")
print(erreurs)
