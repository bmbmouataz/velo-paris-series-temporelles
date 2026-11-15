# 00 - Outils communs à tous les scripts d'analyse (03 à 10)
#
# Ce fichier ne produit rien : il est chargé au début de chaque script par
#   source("R/00_outils.R")
# Il contient :
#   - les packages et le thème graphique commun ;
#   - les dates clés du projet (période d'apprentissage, période de test) ;
#   - la lecture de la série finale avec les variables dérivées ;
#   - les régresseurs exogènes du SARIMAX (calendrier, météo, Fourier) ;
#   - les graphiques ACF / PACF et la fonction de diagnostic des résidus.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(scales)
  library(forecast)
  library(tseries)
})

# ===========================================================================
# Thème et sauvegarde des figures
# ===========================================================================
theme_set(theme_minimal(base_size = 11) +
            theme(plot.title = element_text(face = "bold"),
                  plot.subtitle = element_text(colour = "grey35"),
                  panel.grid.minor = element_blank(),
                  legend.position = "bottom"))

bleu <- "#4C72B0"; rouge <- "#C44E52"; vert <- "#55A868"
orange <- "#DD8452"; violet <- "#8172B3"; gris <- "grey55"
couleurs_evenement <- c("Grève des transports" = "#DD8452", "Confinement 1" = "#C44E52",
                        "Confinement 2" = "#C44E52", "Confinement 3" = "#C44E52",
                        "Jeux olympiques" = "#8172B3")

# Chaque script appelle preparer_figures("03_exploration") puis sauver(g, "01_nom.png").
preparer_figures <- function(dossier) {
  dossier_figures <<- file.path("figures", dossier)
  dir.create(dossier_figures, recursive = TRUE, showWarnings = FALSE)
  dir.create("resultats", showWarnings = FALSE)
}
sauver <- function(g, nom, l = 10, h = 5.5) {
  ggsave(file.path(dossier_figures, nom), g, width = l, height = h, dpi = 130, bg = "white")
  cat("  figure :", nom, "\n")
}
# Tableaux de résultats repris dans le rapport (dossier resultats/)
ecrire_resultat <- function(x, nom) fwrite(x, file.path("resultats", nom))

milliers <- label_number(big.mark = " ", accuracy = 1)
en_pourcent <- function(b) 100 * (exp(b) - 1)   # coefficient d'un modèle en log -> effet en %

# ===========================================================================
# Dates clés
# ===========================================================================
# Apprentissage : de juin 2021 (fin du 3e confinement, retour à un trafic
# « normal », voir R/05_ruptures.R) au 11 novembre 2025.
# Test : la dernière année disponible, du 12 novembre 2025 au 11 novembre 2026.
debut_apprentissage <- as.IDate("2021-06-01")
fin_apprentissage   <- as.IDate("2025-11-11")
debut_test          <- fin_apprentissage + 1L

evenements <- data.table(
  evenement = c("Grève des transports", "Confinement 1", "Confinement 2",
                "Confinement 3", "Jeux olympiques"),
  debut = as.IDate(c("2019-12-05", "2020-03-17", "2020-10-30", "2021-04-03", "2024-07-26")),
  fin   = as.IDate(c("2020-01-19", "2020-05-10", "2020-12-14", "2021-05-02", "2024-08-11"))
)
# Bandes grisées des événements, à ajouter à un graphique dont l'axe x est une date
bandes_evenements <- function(alpha = 0.18) {
  geom_rect(data = evenements, inherit.aes = FALSE,
            aes(xmin = as.Date(debut), xmax = as.Date(fin), ymin = -Inf, ymax = Inf, fill = evenement),
            alpha = alpha)
}

# ===========================================================================
# Lecture de la série finale et variables dérivées
# ===========================================================================
jours_fr <- c("lundi", "mardi", "mercredi", "jeudi", "vendredi", "samedi", "dimanche")
mois_fr  <- c("janv.", "févr.", "mars", "avr.", "mai", "juin", "juil.", "août",
              "sept.", "oct.", "nov.", "déc.")

charger_serie <- function() {
  d <- fread("data/3_final/velo_paris_par_jour.csv")
  d[, date := as.IDate(date)]
  d[is.na(evenement), evenement := ""]
  d[is.na(nom_vacances), nom_vacances := ""]
  d[is.na(nom_jour_ferie), nom_jour_ferie := ""]
  d[, `:=`(
    log_velos   = log(velos),
    annee       = year(date),
    mois        = factor(mois_fr[month(date)], levels = mois_fr),
    jour_semaine = factor(jour_semaine, levels = jours_fr),
    jour_annee  = yday(date),
    t           = as.numeric(date - as.IDate("2019-01-01"))   # temps en jours
  )]
  # Jour férié tombant en semaine (le samedi et le dimanche sont déjà « bas »)
  d[, ferie_semaine := est_ferie & !est_week_end]
  # Pont : lundi précédant un mardi férié, ou vendredi suivant un jeudi férié
  d[, pont := (jour_semaine == "lundi"    & shift(ferie_semaine, -1L, fill = FALSE)) |
              (jour_semaine == "vendredi" & shift(ferie_semaine,  1L, fill = FALSE))]
  d[, type_vacances := fcase(nom_vacances == "Vacances de Noël", "Noël",
                             nom_vacances == "Vacances d'Été", "Été",
                             nom_vacances != "", "Autres vacances",
                             default = "Hors vacances")]
  d[, periode := fcase(date <  debut_apprentissage, "Avant (2019 - mai 2021)",
                       date <= fin_apprentissage,   "Apprentissage",
                       default = "Test")]
  d[]
}

# ===========================================================================
# Régresseurs exogènes (SARIMAX, modèle additif)
# ===========================================================================
# - Fourier : K paires sin/cos de période 365,25 jours (saison annuelle).
# - Météo : température et son carré (effet en cloche : trop chaud fait baisser),
#   pluie en log(1 + mm) (le premier millimètre compte plus que le dixième), vent.
# - Calendrier : férié en semaine, pont, trois types de vacances, Jeux olympiques.
# - Semaine (option) : indicatrices des jours, quand la saison 7 est traitée comme déterministe.
fourier_annuel <- function(t, K) {
  m <- do.call(cbind, lapply(seq_len(K), function(k)
    cbind(sin(2 * pi * k * t / 365.25), cos(2 * pi * k * t / 365.25))))
  colnames(m) <- paste0(rep(c("sin", "cos"), K), "_", rep(seq_len(K), each = 2))
  m
}

regresseurs <- function(d, K = 4, meteo = TRUE, calendrier = TRUE, semaine = FALSE, colonnes = NULL) {
  if (!is.null(colonnes)) semaine <- TRUE   # pour retrouver toutes les colonnes possibles
  X <- fourier_annuel(d$t, K)
  # Saison hebdomadaire déterministe : une indicatrice par jour (lundi = référence)
  if (semaine) {
    S <- sapply(jours_fr[-1], function(j) as.numeric(d$jour_semaine == j))
    colnames(S) <- jours_fr[-1]
    X <- cbind(X, S)
  }
  if (calendrier) {
    X <- cbind(X,
               ferie         = as.numeric(d$ferie_semaine),
               pont          = as.numeric(d$pont),
               vacances_noel = as.numeric(d$type_vacances == "Noël"),
               vacances_ete  = as.numeric(d$type_vacances == "Été"),
               autres_vacances = as.numeric(d$type_vacances == "Autres vacances"),
               jeux_olympiques = as.numeric(d$evenement == "Jeux olympiques"))
  }
  if (meteo) {
    X <- cbind(X,
               temperature   = d$temperature_moyenne - 15,
               temperature_2 = (d$temperature_moyenne - 15)^2 / 10,
               pluie         = log1p(d$pluie_mm),
               vent          = d$vent_moyen_m_s)
  }
  # Pour prévoir, on garde exactement les colonnes du modèle estimé (`colonnes`).
  # Sinon on retire les colonnes constantes (ex. JO absents de la période).
  if (!is.null(colonnes)) return(X[, colonnes, drop = FALSE])
  X[, apply(X, 2, function(x) length(unique(x)) > 1), drop = FALSE]
}

# Météo « climatologique » : ce qu'on sait de la météo sans prévision.
# Pour chaque jour de l'année, moyenne des années de référence, lissée sur ±7 jours.
# La pluie est moyennée en log(1 + mm), l'échelle où elle entre dans le modèle.
climatologie <- function(cible, reference) {
  ref <- reference[, .(temperature_moyenne = mean(temperature_moyenne), pluie_log = mean(log1p(pluie_mm)),
                       vent_moyen_m_s = mean(vent_moyen_m_s)), by = jour_annee][order(jour_annee)]
  lisser <- function(x) { n <- length(x); xx <- c(tail(x, 7), x, head(x, 7))
                          frollmean(xx, 15, align = "center")[8:(n + 7)] }
  cols <- c("temperature_moyenne", "pluie_log", "vent_moyen_m_s")
  ref[, (cols) := lapply(.SD, lisser), .SDcols = cols]
  out <- copy(cible)
  out[ref, on = "jour_annee", `:=`(temperature_moyenne = i.temperature_moyenne,
                                   pluie_mm = expm1(i.pluie_log), vent_moyen_m_s = i.vent_moyen_m_s)]
  out
}

noms_lisibles <- c(
  mardi = "Mardi (vs lundi)", mercredi = "Mercredi (vs lundi)", jeudi = "Jeudi (vs lundi)",
  vendredi = "Vendredi (vs lundi)", samedi = "Samedi (vs lundi)", dimanche = "Dimanche (vs lundi)",
  ferie = "Jour férié (en semaine)", pont = "Pont",
  vacances_noel = "Vacances de Noël", vacances_ete = "Vacances d'été",
  autres_vacances = "Autres vacances", jeux_olympiques = "Jeux olympiques",
  temperature = "Température (+1 °C autour de 15 °C)", temperature_2 = "Température² / 10",
  pluie = "Pluie (log(1 + mm), +1)", vent = "Vent (+1 m/s)"
)

# ===========================================================================
# ACF / PACF au format ggplot
# ===========================================================================
# Les retards multiples de la saison (7 jours) sont en rouge.
gg_acf <- function(x, lag.max = 56, type = c("acf", "pacf"), titre = NULL, saison = 7) {
  type <- match.arg(type)
  x <- as.numeric(na.omit(x)); n <- length(x)
  v <- if (type == "acf") acf(x, lag.max = lag.max, plot = FALSE)$acf[-1]
       else pacf(x, lag.max = lag.max, plot = FALSE)$acf
  d <- data.table(retard = seq_along(v), valeur = as.numeric(v))
  d[, saisonnier := ifelse(retard %% saison == 0, paste("multiple de", saison), "autre retard")]
  b <- 1.96 / sqrt(n)
  ggplot(d, aes(retard, valeur)) +
    geom_hline(yintercept = 0, colour = "grey40") +
    geom_hline(yintercept = c(-b, b), linetype = 2, colour = bleu) +
    geom_segment(aes(xend = retard, yend = 0, colour = saisonnier), linewidth = 0.7) +
    scale_colour_manual(values = c("autre retard" = "grey30", setNames(rouge, paste("multiple de", saison))),
                        name = NULL) +
    scale_x_continuous(breaks = seq(0, lag.max, by = saison)) +
    labs(title = titre, x = "Retard (jours)", y = if (type == "acf") "ACF" else "PACF")
}

# ===========================================================================
# Diagnostic des résidus : diagnostic complet des résidus
# ===========================================================================
# Produit une figure en 6 panneaux et renvoie les tests :
#   - Ljung-Box (blancheur), avec fitdf = nombre de paramètres ARMA ;
#   - Ljung-Box sur les résidus au carré (effet ARCH / GARCH) ;
#   - Shapiro-Wilk et Jarque-Bera (normalité).
diagnostic_residus <- function(res, dates, titre, fichier, fitdf = 0, lblag = 21, lagmax = 42) {
  ok <- !is.na(res); r <- as.numeric(res[ok]); dt <- as.Date(dates[ok])
  lb  <- Box.test(r,   lag = lblag, type = "Ljung-Box", fitdf = fitdf)
  lb2 <- Box.test(r^2, lag = lblag, type = "Ljung-Box")
  sw  <- shapiro.test(if (length(r) > 5000) r[seq_len(5000)] else r)
  jb  <- jarque.bera.test(r)
  tests <- data.table(
    test = c(sprintf("Ljung-Box (h = %d, fitdf = %d)", lblag, fitdf),
             sprintf("Ljung-Box résidus² (h = %d)", lblag), "Shapiro-Wilk", "Jarque-Bera"),
    hypothese_nulle = c("résidus non corrélés", "pas d'effet ARCH", "normalité", "normalité"),
    statistique = c(lb$statistic, lb2$statistic, sw$statistic, jb$statistic),
    p_valeur = c(lb$p.value, lb2$p.value, sw$p.value, jb$p.value)
  )
  # p-valeurs de Ljung-Box pour chaque horizon h (comme tsdiag)
  h <- (fitdf + 1):lagmax
  pv <- data.table(h = h, p = sapply(h, function(k) Box.test(r, lag = k, type = "Ljung-Box", fitdf = fitdf)$p.value))
  dr <- data.table(date = dt, r = r)
  sd_r <- sd(r)
  extremes <- dr[order(-abs(r))][1:5]

  # Résidus standardisés : hors de ±1,96, environ 5 % des points si gaussiens
  g1 <- ggplot(dr, aes(date, r / sd_r)) + geom_line(colour = "grey40", linewidth = 0.3) +
    geom_hline(yintercept = c(-1.96, 1.96), linetype = 2, colour = rouge) +
    geom_point(data = extremes, colour = rouge) +
    ggrepel::geom_text_repel(data = extremes, aes(label = format(date, "%d/%m/%y")), size = 2.8) +
    labs(title = sprintf("Résidus standardisés (%.1f %% hors de ±1,96)", 100 * mean(abs(r / sd_r) > 1.96)),
         x = NULL, y = NULL)
  # Nuage Res[i] contre Res[i-1] : aucune structure si non corrélés
  g8 <- ggplot(data.table(avant = head(r, -1), apres = r[-1]), aes(avant, apres)) +
    geom_point(size = 0.5, alpha = 0.4, colour = bleu) +
    geom_smooth(method = "lm", formula = y ~ x, colour = rouge, se = FALSE, linewidth = 0.7) +
    labs(title = sprintf("Res[i] contre Res[i-1] (r = %.2f)", cor(head(r, -1), r[-1])), x = "Res[i-1]", y = "Res[i]")
  g2 <- gg_acf(r, lagmax, "acf", "ACF des résidus") + theme(legend.position = "none")
  g3 <- gg_acf(r, lagmax, "pacf", "PACF des résidus") + theme(legend.position = "none")
  g4 <- ggplot(dr, aes(r)) + geom_histogram(aes(y = after_stat(density)), bins = 60, fill = "grey70") +
    stat_function(fun = dnorm, args = list(mean = mean(r), sd = sd_r), colour = rouge, linewidth = 0.8) +
    labs(title = "Histogramme et loi normale", x = NULL, y = NULL)
  g5 <- ggplot(dr, aes(sample = r)) + stat_qq(size = 0.6, colour = "grey30") +
    stat_qq_line(colour = rouge) + labs(title = "QQ-plot normal", x = "Quantiles théoriques", y = NULL)
  g6 <- ggplot(pv, aes(h, p)) + geom_point(colour = bleu) +
    geom_hline(yintercept = 0.05, linetype = 2, colour = rouge) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "p-valeurs de Ljung-Box selon h", x = "h", y = NULL)
  g7 <- gg_acf(r^2, lagmax, "acf", "ACF des résidus² (effet ARCH ?)") + theme(legend.position = "none")

  sous_titre <- sprintf("Ljung-Box p = %s | Ljung-Box sur résidus² p = %s | Shapiro p = %s | Jarque-Bera p = %s",
                        format.pval(lb$p.value, 2), format.pval(lb2$p.value, 2),
                        format.pval(sw$p.value, 2), format.pval(jb$p.value, 2))
  g <- (g1 | g6) / (g2 | g3 | g8) / (g7 | g4 | g5) +
    plot_annotation(title = titre, subtitle = sous_titre,
                    theme = theme(plot.title = element_text(face = "bold", size = 14)))
  sauver(g, fichier, l = 12, h = 11)
  invisible(tests)
}

# Tableau des coefficients avec test de significativité (z = coef / écart-type)
table_coefficients <- function(modele) {
  b <- coef(modele); s <- sqrt(diag(modele$var.coef))
  data.table(parametre = names(b), estimation = b, ecart_type = s,
             z = b / s, p_valeur = 2 * pnorm(-abs(b / s)))
}

# Inverses des racines des polynômes AR et MA, avec le cercle unité
gg_racines <- function(modele, titre) {
  # Polynôme AR : 1 - a1 z - ... ; polynôme MA : 1 + b1 z + ...
  r <- function(coefs, signe, nom) {
    if (length(coefs) == 0 || all(coefs == 0)) return(NULL)
    z <- polyroot(c(1, signe * coefs))
    data.table(re = Re(1 / z), im = Im(1 / z), polynome = nom)
  }
  s  <- modele$arma[5]   # période saisonnière
  cf <- coef(modele)
  ar  <- cf[grep("^ar", names(cf))];  ma  <- cf[grep("^ma", names(cf))]
  sar <- cf[grep("^sar", names(cf))]; sma <- cf[grep("^sma", names(cf))]
  # Un polynôme saisonnier en B^s : coefficients placés aux puissances s, 2s, ...
  en_saison <- function(x) if (length(x) == 0) numeric(0) else {
    v <- numeric(s * length(x)); v[s * seq_along(x)] <- x; v }
  d <- rbind(r(ar, -1, "AR"), r(ma, 1, "MA"),
             r(en_saison(sar), -1, "AR saisonnier"), r(en_saison(sma), 1, "MA saisonnier"))
  cercle <- data.table(a = seq(0, 2 * pi, length.out = 200))
  ggplot(d, aes(re, im, colour = polynome)) +
    geom_path(data = cercle, aes(cos(a), sin(a)), inherit.aes = FALSE, colour = "grey50") +
    geom_hline(yintercept = 0, colour = "grey85") + geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(size = 2.5) + coord_equal() +
    scale_colour_manual(values = c("AR" = bleu, "MA" = rouge, "AR saisonnier" = vert, "MA saisonnier" = orange),
                        name = NULL) +
    labs(title = titre, subtitle = "Inverses des racines : à l'intérieur du cercle = causal et inversible",
         x = "Partie réelle", y = "Partie imaginaire")
}
