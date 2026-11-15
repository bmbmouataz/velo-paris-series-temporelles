# 05 - Ruptures et événements exceptionnels
#
# Entrée  : data/3_final/velo_paris_par_jour.csv
# Sorties : figures/05_ruptures/*.png, resultats/05_*.csv
#
# Trois questions :
#   1. Où la série change-t-elle de régime ? Détection par minimisation de
#      l'erreur quadratique, puis ruptures multiples (Bai-Perron).
#   2. La période d'apprentissage choisie (juin 2021 - novembre 2025) est-elle
#      stable ? Test CUSUM.
#   3. Combien chaque événement a-t-il coûté (ou rapporté) en vélos ?
#      Analyse d'intervention : une indicatrice par événement dans une régression.
#
# On travaille sur la série corrigée des variations saisonnières (CVS) :
# log(vélos) moins les deux saisons estimées par MSTL, moyennée par semaine.

source("R/00_outils.R")
suppressPackageStartupMessages(library(strucchange))
preparer_figures("05_ruptures")
d <- charger_serie()

ms <- mstl(msts(d$log_velos, seasonal.periods = c(7, 365.25)), s.window = 21, robust = TRUE)
d[, cvs := log_velos - as.numeric(ms[, "Seasonal7"]) - as.numeric(ms[, "Seasonal365.25"])]
hebdo <- d[, .(cvs = mean(cvs), n = .N), by = .(semaine = as.IDate(cut(as.Date(date), "week")))][n == 7]
hebdo[, i := .I]

# ===========================================================================
# 1. Une rupture : minimiser l'erreur quadratique
# ===========================================================================
# Pour chaque date candidate τ, on ajuste une moyenne avant et une moyenne après
# et on calcule la somme des carrés des résidus SCR(τ). La rupture estimée est
# le τ qui minimise SCR(τ). Variante : une droite avant et une droite après.
n <- nrow(hebdo); marge <- 26   # au moins 6 mois de chaque côté (sinon une droite colle à quelques points)
candidats <- (marge + 1):(n - marge)
scr <- rbindlist(lapply(candidats, function(tau) {
  avant <- hebdo$cvs[1:tau]; apres <- hebdo$cvs[(tau + 1):n]
  m_av <- lm(cvs ~ i, hebdo[1:tau]); m_ap <- lm(cvs ~ i, hebdo[(tau + 1):n])
  data.table(semaine = hebdo$semaine[tau],
             `Moyenne avant / après` = sum((avant - mean(avant))^2) + sum((apres - mean(apres))^2),
             `Droite avant / après` = sum(residuals(m_av)^2) + sum(residuals(m_ap)^2))
}))
scr_long <- melt(scr, id.vars = "semaine", variable.name = "modele", value.name = "scr")
meilleurs <- scr_long[, .SD[which.min(scr)], by = modele]
g1 <- ggplot(hebdo, aes(as.Date(semaine), cvs)) +
  bandes_evenements(0.25) + geom_line(colour = "grey30") +
  geom_vline(data = meilleurs, aes(xintercept = as.Date(semaine), colour = modele), linewidth = 1, linetype = 2) +
  scale_fill_manual(values = couleurs_evenement, name = NULL) +
  scale_colour_manual(values = c(bleu, rouge), name = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Série corrigée des saisons (moyenne par semaine)", x = NULL, y = "log(vélos) CVS") +
  guides(fill = "none")
g2 <- ggplot(scr_long, aes(as.Date(semaine), scr, colour = modele)) + geom_line(linewidth = 0.8) +
  geom_point(data = meilleurs, size = 3) +
  geom_label(data = meilleurs, aes(label = format(as.Date(semaine), "%d/%m/%Y")), vjust = -0.5, size = 3, show.legend = FALSE) +
  scale_colour_manual(values = c(bleu, rouge), name = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Somme des carrés des résidus selon la date de rupture τ",
       subtitle = "Le minimum désigne la rupture la plus nette", x = "Date de rupture candidate", y = "SCR(τ)")
g <- g1 / g2 + plot_annotation(title = "Détection d'une rupture par moindres carrés")
sauver(g, "01_rupture_unique_moindres_carres.png", h = 8)

# ===========================================================================
# 2. Plusieurs ruptures : Bai-Perron (programmation dynamique), choix par BIC
# ===========================================================================
bp <- breakpoints(cvs ~ 1, data = hebdo, h = 0.06)
nb <- summary(bp)$RSS["BIC", ]
m_opt <- as.integer(names(nb)[which.min(nb)])
bp_opt <- breakpoints(bp, breaks = m_opt)
ic <- tryCatch(confint(bp_opt)$confint, error = function(e) NULL)
dates_rupt <- hebdo$semaine[bp_opt$breakpoints]
hebdo[, segment := factor(breakfactor(bp_opt))]
hebdo[, moyenne_segment := mean(cvs), by = segment]
cat("Ruptures (Bai-Perron, BIC) :", format(dates_rupt), "\n")
g1 <- ggplot(data.table(m = as.integer(names(nb)), bic = nb), aes(m, bic)) +
  geom_line(colour = bleu) + geom_point(colour = bleu, size = 2) +
  geom_point(data = data.table(m = m_opt, bic = min(nb)), colour = rouge, size = 4) +
  labs(title = "BIC selon le nombre de ruptures", subtitle = sprintf("Minimum pour %d ruptures", m_opt),
       x = "Nombre de ruptures", y = "BIC")
intervalles <- if (!is.null(ic)) data.table(debut = hebdo$semaine[pmax(1, ic[, 1])], fin = hebdo$semaine[pmin(n, ic[, 3])]) else NULL
g2 <- ggplot(hebdo, aes(as.Date(semaine))) +
  bandes_evenements(0.25) +
  {if (!is.null(intervalles)) geom_rect(data = intervalles, inherit.aes = FALSE,
                                        aes(xmin = as.Date(debut), xmax = as.Date(fin), ymin = -Inf, ymax = Inf),
                                        fill = "grey50", alpha = 0.25)} +
  geom_line(aes(y = cvs), colour = "grey40") +
  geom_step(aes(y = moyenne_segment), colour = rouge, linewidth = 1.1) +
  geom_vline(xintercept = as.Date(dates_rupt), colour = rouge, linetype = 2) +
  annotate("label", x = as.Date(dates_rupt), y = max(hebdo$cvs) + 0.05, label = format(as.Date(dates_rupt), "%b %Y"),
           size = 2.8, colour = rouge) +
  scale_fill_manual(values = couleurs_evenement, name = NULL) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Ruptures de niveau détectées (Bai-Perron)",
       subtitle = "Rouge : moyenne de chaque régime. Gris foncé : intervalle de confiance à 95 % de la date",
       x = NULL, y = "log(vélos) CVS") + guides(fill = guide_legend(nrow = 1))
g <- (g1 | g2) + plot_layout(widths = c(1, 2.6))
sauver(g, "02_ruptures_multiples.png", l = 13, h = 5.5)

# ===========================================================================
# 3. Stabilité de la période d'apprentissage : CUSUM
# ===========================================================================
# Sous H0 (pas de rupture), le processus CUSUM des résidus récursifs reste dans
# les bandes. On teste un modèle « tendance linéaire » sur deux périodes.
cusum_periode <- function(debut, fin, nom) {
  h <- hebdo[semaine >= debut & semaine <= fin]
  e <- efp(cvs ~ i, data = h, type = "OLS-CUSUM")
  test <- sctest(e)
  bande <- boundary(e, alpha = 0.05)
  list(dt = data.table(semaine = h$semaine, cusum = as.numeric(e$process)[-1][seq_len(nrow(h))],
                       bande = as.numeric(bande)[seq_len(nrow(h))], periode = nom),
       p = test$p.value)
}
c1 <- cusum_periode(as.IDate("2019-01-01"), fin_apprentissage, "2019 - 2025 (avec la crise sanitaire)")
c2 <- cusum_periode(debut_apprentissage, fin_apprentissage, "Juin 2021 - novembre 2025 (apprentissage)")
cus <- rbind(c1$dt, c2$dt)
pv <- data.table(periode = c(c1$dt$periode[1], c2$dt$periode[1]), p = c(c1$p, c2$p))
g <- ggplot(cus, aes(as.Date(semaine))) +
  geom_ribbon(aes(ymin = -bande, ymax = bande), fill = "#A1C9F4", alpha = 0.5) +
  geom_hline(yintercept = 0, colour = "grey40") +
  geom_line(aes(y = cusum), colour = rouge, linewidth = 0.8) +
  geom_text(data = pv, aes(x = as.Date(-Inf), y = Inf, label = sprintf("Test OLS-CUSUM : p = %s", format.pval(p, 2))),
            hjust = -0.1, vjust = 1.5, fontface = "bold") +
  facet_wrap(~periode, ncol = 1, scales = "free_x") +
  scale_x_date(date_breaks = "6 months", date_labels = "%m/%Y") +
  labs(title = "Test de stabilité CUSUM : la série est-elle homogène sur la période ?",
       subtitle = "Modèle : tendance linéaire sur la série CVS hebdomadaire. Bleu : bande à 95 %. Sortir de la bande = rupture",
       x = NULL, y = "Processus CUSUM")
sauver(g, "03_cusum_stabilite.png", h = 7.5)

# ===========================================================================
# 4. Analyse d'intervention : effet de chaque événement
# ===========================================================================
# Régression sur toute la période 2019-2026, en log :
#   log X_t = tendance (spline) + jour de semaine + Fourier annuel + calendrier + météo
#             + Σ_e ω_e 1{t dans l'événement e} + z_t
# ω_e mesure l'effet de l'événement e, toutes choses égales par ailleurs ;
# exp(ω_e) - 1 est l'effet en %. On corrige l'écart-type de l'autocorrélation
# des résidus par Newey-West (sinon les intervalles seraient trop étroits).
X <- regresseurs(d, K = 6, meteo = TRUE, calendrier = TRUE)
X <- X[, colnames(X) != "jeux_olympiques"]
ev_mat <- sapply(evenements$evenement, function(e) as.numeric(d$evenement == e))
colnames(ev_mat) <- make.names(evenements$evenement)
tendance <- splines::ns(d$t, df = 8)
m_int <- lm(d$log_velos ~ tendance + d$jour_semaine + X + ev_mat)
vc <- sandwich::NeweyWest(m_int, lag = 14)
b <- coef(m_int)[grep("ev_mat", names(coef(m_int)))]
s <- sqrt(diag(vc))[grep("ev_mat", names(coef(m_int)))]
effets <- data.table(evenement = evenements$evenement, debut = evenements$debut, fin = evenements$fin,
                     effet = en_pourcent(b), bas = en_pourcent(b - 1.96 * s), haut = en_pourcent(b + 1.96 * s),
                     p_valeur = 2 * pnorm(-abs(b / s)))
effets[, evenement := factor(evenement, levels = rev(evenements$evenement))]
g1 <- ggplot(effets, aes(effet, evenement, colour = effet > 0)) +
  geom_vline(xintercept = 0, colour = "grey40") +
  geom_errorbarh(aes(xmin = bas, xmax = haut), height = 0.25, linewidth = 0.8) +
  geom_point(size = 3.5) +
  geom_text(aes(label = sprintf("%+.0f %%", effet)), vjust = -1.1, size = 3.5, show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = vert, "FALSE" = rouge), guide = "none") +
  labs(title = "Effet de chaque événement (à météo et calendrier égaux)",
       subtitle = "Intervalle à 95 % (écarts-types de Newey-West)", x = "Effet sur le trafic (%)", y = NULL)
# Contrefactuel : ce qu'aurait été le trafic sans les événements
d[, ajuste := fitted(m_int)]
d[, sans_evenement := ajuste - as.numeric(ev_mat %*% b)]
zoom <- d[date >= as.IDate("2019-10-01") & date <= as.IDate("2021-07-01")]
g2 <- ggplot(zoom, aes(as.Date(date))) +
  geom_rect(data = evenements[debut < as.IDate("2021-07-01")], inherit.aes = FALSE,
            aes(xmin = as.Date(debut), xmax = as.Date(fin), ymin = -Inf, ymax = Inf), fill = "grey88") +
  geom_line(aes(y = frollmean(velos, 7, align = "center"), colour = "Observé"), linewidth = 0.8) +
  geom_line(aes(y = frollmean(exp(sans_evenement), 7, align = "center"), colour = "Sans l'événement (modèle)"),
            linewidth = 0.8, linetype = 2) +
  scale_colour_manual(values = c("Observé" = rouge, "Sans l'événement (modèle)" = bleu), name = NULL) +
  scale_y_continuous(labels = milliers) + scale_x_date(date_breaks = "2 months", date_labels = "%b %y") +
  labs(title = "Contrefactuel : le trafic sans la grève ni les confinements (moyenne 7 jours)",
       subtitle = "L'écart entre les deux courbes est l'effet de l'événement", x = NULL, y = "Vélos par jour")
g <- g1 / g2 + plot_layout(heights = c(1, 1.3))
sauver(g, "04_analyse_intervention.png", h = 9)

# Vélos gagnés ou perdus pendant chaque événement
bilan <- d[evenement != "", .(jours = .N, observe = sum(velos), sans = sum(exp(sans_evenement))), by = evenement]
bilan[, ecart := observe - sans]
effets <- merge(effets, bilan, by = "evenement")

# ===========================================================================
# 5. Conséquence pour la modélisation
# ===========================================================================
g <- ggplot(d, aes(as.Date(date), velos)) +
  annotate("rect", xmin = as.Date(debut_apprentissage), xmax = as.Date(fin_apprentissage),
           ymin = -Inf, ymax = Inf, fill = bleu, alpha = 0.12) +
  annotate("rect", xmin = as.Date(debut_test), xmax = as.Date(max(d$date)),
           ymin = -Inf, ymax = Inf, fill = orange, alpha = 0.2) +
  geom_line(colour = "grey45", linewidth = 0.25) +
  geom_vline(xintercept = as.Date(dates_rupt), colour = rouge, linetype = 2) +
  annotate("text", x = as.Date(c("2023-08-15", "2026-05-15")), y = 95000,
           label = c("Apprentissage\njuin 2021 - nov. 2025", "Test\n365 jours"), size = 3.5, fontface = "bold") +
  annotate("text", x = as.Date("2020-03-01"), y = 95000, label = "Exclu :\nrégime de crise", size = 3.5, colour = rouge) +
  scale_y_continuous(labels = milliers, limits = c(0, 100000)) +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(title = "Découpage retenu pour la modélisation",
       subtitle = "Traits rouges : ruptures de niveau (Bai-Perron). Après juin 2021 la série suit une tendance stable (CUSUM, p = 0,65) ; la dernière année sert de test",
       x = NULL, y = "Vélos par jour")
sauver(g, "05_decoupage_apprentissage_test.png", h = 5)

ecrire_resultat(effets[order(debut), .(evenement, debut, fin, jours, effet_pourcent = round(effet, 1),
                                       ic_bas = round(bas, 1), ic_haut = round(haut, 1), p_valeur,
                                       velos_gagnes_ou_perdus = round(ecart))], "05_effets_evenements.csv")
ecrire_resultat(data.table(rupture = format(dates_rupt)), "05_dates_ruptures.csv")
ecrire_resultat(pv, "05_cusum.csv")
print(effets[order(debut), .(evenement, effet, bas, haut, ecart)]); print(pv)
