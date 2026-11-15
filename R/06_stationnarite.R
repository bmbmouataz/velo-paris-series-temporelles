# 06 - Stationnarité et différenciations
#
# Entrée  : data/3_final/velo_paris_par_jour.csv (période d'apprentissage)
# Sorties : figures/06_stationnarite/*.png, resultats/06_tests_stationnarite.csv
#
# Une série doit être (faiblement) stationnaire pour être modélisée par un ARMA :
# moyenne et variance constantes, autocovariance qui ne dépend que de l'écart.
# On teste cinq versions de log(vélos) :
#   y                      la série
#   (I - B) y              différence d'ordre 1        -> enlève une tendance
#   (I - B^7) y            différence saisonnière      -> enlève la semaine
#   (I - B)(I - B^7) y     les deux
#   résidu de régression   y moins tendance + Fourier + calendrier + météo (base du SARIMAX)
#
# Tests (attention au sens des hypothèses) :
#   ADF  : H0 = racine unitaire (non stationnaire). p petite -> stationnaire.
#   PP   : idem, correction non paramétrique.
#   KPSS : H0 = stationnaire. p petite -> non stationnaire.
# Conclusion fiable quand ADF rejette ET KPSS ne rejette pas.

source("R/00_outils.R")
preparer_figures("06_stationnarite")
d <- charger_serie()
app <- d[periode == "Apprentissage"]

y <- app$log_velos
reg <- lm(y ~ app$t + regresseurs(app, K = 6))
versions <- list(
  "y = log(vélos)"              = y,
  "(I - B) y"                    = diff(y, 1),
  "(I - B⁷) y"                   = diff(y, 7),
  "(I - B)(I - B⁷) y"            = diff(diff(y, 7), 1),
  "Résidu de régression (tendance + Fourier + calendrier + météo)" = residuals(reg)
)

# ===========================================================================
# 1. Tests
# ===========================================================================
tests <- rbindlist(lapply(names(versions), function(nom) {
  x <- versions[[nom]]
  adf  <- suppressWarnings(adf.test(x, k = 14))
  pp   <- suppressWarnings(pp.test(x))
  kpl  <- suppressWarnings(kpss.test(x, null = "Level"))
  kpt  <- suppressWarnings(kpss.test(x, null = "Trend"))
  data.table(serie = nom,
             test = c("ADF (H0 : racine unitaire)", "Phillips-Perron (H0 : racine unitaire)",
                      "KPSS niveau (H0 : stationnaire)", "KPSS tendance (H0 : stationnaire)"),
             statistique = c(adf$statistic, pp$statistic, kpl$statistic, kpt$statistic),
             p_valeur = c(adf$p.value, pp$p.value, kpl$p.value, kpt$p.value),
             variance = var(x))
}))
# tseries borne les p-valeurs : 0,01 signifie « < 0,01 » et 0,1 « > 0,1 »
tests[, lecture := fifelse(grepl("KPSS", test),
                           fifelse(p_valeur < 0.05, "non stationnaire", "stationnaire"),
                           fifelse(p_valeur < 0.05, "stationnaire", "non stationnaire"))]
tests[, serie := factor(serie, levels = names(versions))]
print(tests[, .(serie, test, p_valeur, lecture)])
ecrire_resultat(tests, "06_tests_stationnarite.csv")

g <- ggplot(tests, aes(test, serie, fill = lecture)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = paste0("p ", ifelse(p_valeur <= 0.01, "< 0,01", ifelse(p_valeur >= 0.1, "> 0,10",
                                                                               formatC(p_valeur, format = "f", digits = 3))),
                               "\n", lecture)), size = 3.3) +
  scale_fill_manual(values = c("stationnaire" = "#B5D99C", "non stationnaire" = "#F4A6A6"), name = NULL) +
  scale_y_discrete(limits = rev(names(versions))) +
  scale_x_discrete(labels = function(x) gsub(" \\(", "\n(", x)) +
  labs(title = "Tests de stationnarité sur la période d'apprentissage",
       subtitle = "ADF et PP : H0 = racine unitaire. KPSS : H0 = stationnarité. Vert = conclusion « stationnaire » au seuil de 5 %",
       x = NULL, y = NULL) +
  theme(panel.grid = element_blank(), legend.position = "none")
sauver(g, "01_tableau_tests.png", l = 11, h = 5.5)

# ===========================================================================
# 2. Chaque version : chronogramme, ACF, PACF
# ===========================================================================
lignes <- lapply(names(versions)[1:4], function(nom) {
  x <- versions[[nom]]
  dd <- data.table(date = tail(app$date, length(x)), x = x)
  g1 <- ggplot(dd, aes(as.Date(date), x)) + geom_line(colour = bleu, linewidth = 0.25) +
    labs(title = nom, x = NULL, y = NULL)
  g2 <- gg_acf(x, 42, "acf", "ACF") + theme(legend.position = "none")
  g3 <- gg_acf(x, 42, "pacf", "PACF") + theme(legend.position = "none")
  g1 | g2 | g3
})
g <- wrap_plots(lignes, ncol = 1) +
  plot_annotation(title = "Différencier pour rendre la série stationnaire",
                  subtitle = paste("y : ACF lente avec pics tous les 7 jours (non stationnaire). (I - B) y : il reste la semaine.",
                                   "(I - B⁷) y : la semaine disparaît, ACF qui décroît vite + pic négatif au retard 7 (MA saisonnier)"),
                  caption = "En rouge : retards multiples de 7")
sauver(g, "02_differenciations_acf_pacf.png", l = 13, h = 12)

# ===========================================================================
# 3. Choix automatique du nombre de différences (forecast)
# ===========================================================================
y7 <- ts(y, frequency = 7)
choix <- data.table(
  question = c("d (KPSS)", "d (ADF)", "d (PP)", "D, saison 7 (test OCSB)", "D, saison 7 (test de saisonnalité)",
               "d après (I - B⁷) (KPSS)"),
  reponse = c(ndiffs(y, test = "kpss"), ndiffs(y, test = "adf"), ndiffs(y, test = "pp"),
              nsdiffs(y7, test = "ocsb"), nsdiffs(y7, test = "seas"), ndiffs(diff(y, 7), test = "kpss")))
print(choix)
ecrire_resultat(choix, "06_nombre_differences.csv")

# ===========================================================================
# 4. Ne pas sur-différencier : la variance comme garde-fou
# ===========================================================================
# Différencier une série déjà stationnaire augmente sa variance et crée une
# racine unitaire dans la partie MA (non inversible).
var_dt <- unique(tests[, .(serie, variance)])
var_dt[, serie := factor(serie, levels = names(versions))]
g1 <- ggplot(var_dt, aes(serie, variance, fill = serie)) + geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = sprintf("%.3f", variance)), vjust = -0.4) +
  scale_fill_manual(values = c(gris, orange, vert, rouge, bleu)) +
  scale_x_discrete(labels = label_wrap(14)) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = "Variance de chaque version", x = NULL, y = "Variance")
# Moyenne et écart-type glissants sur 8 semaines : constants si stationnaire
gl <- rbindlist(lapply(names(versions)[c(1, 3, 5)], function(nom) {
  x <- versions[[nom]]
  data.table(date = tail(app$date, length(x)), serie = sub(" [(].*", "", nom),
             moyenne = frollmean(x, 56, align = "center"),
             ecart_type = frollapply(x, 56, sd, align = "center"))
}))
gl_long <- melt(gl, id.vars = c("date", "serie"))
gl_long[, variable := fifelse(variable == "moyenne", "Moyenne sur 8 semaines", "Écart-type sur 8 semaines")]
g2 <- ggplot(gl_long, aes(as.Date(date), value, colour = serie)) + geom_line(linewidth = 0.6) +
  facet_wrap(variable ~ serie, scales = "free_y", ncol = 3, labeller = labeller(serie = label_wrap_gen(30), variable = label_wrap_gen(30), .multi_line = FALSE)) +
  scale_colour_manual(values = c(gris, vert, bleu), name = NULL) +
  labs(title = "Moments glissants (fenêtre de 8 semaines)",
       x = NULL, y = NULL) + theme(legend.position = "none")
v <- setNames(var_dt$variance, as.character(var_dt$serie))
g <- g1 + g2 + plot_layout(widths = c(1, 1.8)) +
  plot_annotation(title = "Variance et moments glissants des versions de la série",
                  subtitle = sprintf(paste("Variance : y %.3f ; (I - B⁷) y %.3f ; (I - B)(I - B⁷) y %.3f. Les tests déclarent déjà (I - B⁷) y stationnaire (d = 0 suffit).",
                                           "
La moyenne de y suit la saison ; (I - B⁷) y et le résidu restent centrés. L'écart-type monte à chaque Noël et chaque mois d'août"),
                                     v[1], v[3], v[4]))
sauver(g, "03_variance_et_moments.png", l = 13, h = 6)

# ===========================================================================
# 5. ADF avec aTSA : les trois types de régression
# ===========================================================================
# aTSA::adf.test donne trois tests, à lire selon la forme de la série :
#   type 1 : ΔX_t = ρ X_{t-1} + Σ φ_i ΔX_{t-i} + ε_t              (ni constante ni tendance)
#   type 2 : + constante α                                         (« drift »)
#   type 3 : + constante α + tendance β t                          (« drift + trend »)
# pour chaque nombre de retards augmentés (lag 0 à 7). Le bon type est celui qui
# correspond au graphique : y a une tendance -> type 3 ; (I - B⁷) y est centrée -> type 1 ou 2.
atsa <- rbindlist(lapply(names(versions)[c(1, 3)], function(nom) {
  a <- aTSA::adf.test(versions[[nom]], nlag = 8, output = FALSE)
  rbindlist(lapply(1:3, function(k) data.table(serie = nom, type = c("Type 1 : sans constante", "Type 2 : constante",
                                                                         "Type 3 : constante + tendance")[k],
                                                retards = a[[k]][, "lag"], ADF = a[[k]][, "ADF"], p_valeur = a[[k]][, "p.value"])))
}))
ecrire_resultat(atsa, "06_adf_atsa_trois_types.csv")
g <- ggplot(atsa, aes(retards, p_valeur, colour = type)) +
  geom_hline(yintercept = 0.05, linetype = 2) + geom_line(linewidth = 0.9) + geom_point() +
  facet_wrap(~serie) + scale_colour_manual(values = c(gris, bleu, rouge), name = NULL) +
  scale_x_continuous(breaks = 0:7) +
  labs(title = "Test ADF (aTSA) : p-valeurs selon le type de régression et le nombre de retards",
       subtitle = paste("aTSA borne les p-valeurs à 0,01. Pour y, le type 3 (avec tendance) est le bon : il rejette la racine unitaire.",
                        "\nLe type 1 ne rejette pas pour y car il impose une moyenne nulle, absurde pour un log autour de 10,7 : il faut lire le bon type"),
       x = "Nombre de retards augmentés", y = "p-valeur") + guides(colour = guide_legend(nrow = 1))
sauver(g, "04_adf_trois_types_atsa.png", l = 12, h = 5.5)
print(atsa[retards %in% c(0, 7)])
