# 00 - Fonctions de prévision partagées par R/09_concurrents.R et R/10_evaluation.R
#
#   prevision_blocs()     méthode de ressemblance entre blocs
#   tableaux_lstm()       mise en forme des données pour le réseau LSTM
#
# Convention : `serie` est le log du trafic sur toute la période (apprentissage +
# test) et `o` l'indice du dernier jour connu au moment de prévoir.

# ===========================================================================
# Ressemblance entre blocs
# ===========================================================================
# Idée : chercher dans le passé les périodes de L jours qui ressemblent le plus
# aux L derniers jours, et prolonger ce qui a suivi ces périodes.
#   1. On retire la tendance (droite estimée sur l'apprentissage) : z_t.
#   2. Bloc courant B_o = (z_{o-L+1}, ..., z_o).
#   3. Candidats : blocs B_j finissant le même jour de la semaine que o (pour
#      aligner la semaine) et dont la suite (j+1 ... j+H) est connue.
#   4. Distance euclidienne entre blocs centrés ; on garde les k plus proches,
#      pondérés par 1 / distance.
#   5. Prévision : niveau du bloc courant + moyenne pondérée des suites centrées.
prevision_blocs <- function(z, o, H = 14, L = 14, k = 10) {
  courant <- z[(o - L + 1):o]
  m0 <- mean(courant)
  fins <- seq(o - 7L, L, by = -7L)          # même jour de la semaine que o
  fins <- fins[fins + H <= o]               # suite entièrement connue
  dist <- vapply(fins, function(j) {
    b <- z[(j - L + 1):j]
    sqrt(sum(((b - mean(b)) - (courant - m0))^2))
  }, numeric(1))
  ordre <- order(dist)[seq_len(min(k, length(fins)))]
  w <- 1 / pmax(dist[ordre], 1e-6); w <- w / sum(w)
  suites <- vapply(fins[ordre], function(j) z[j + seq_len(H)] - mean(z[(j - L + 1):j]), numeric(H))
  list(prevision = m0 + as.numeric(suites %*% w), voisins = fins[ordre], poids = w)
}

# ===========================================================================
# LSTM
# ===========================================================================
# Entrées pour prévoir les jours o+1 ... o+H :
#   - séquence des 28 derniers jours : z (trafic sans tendance, standardisé)
#     et les variables du jour (semaine, saison, calendrier, météo) ;
#   - variables des H jours à prévoir (calendrier connu, météo observée).
# Sortie : z des H jours suivants (une sortie par horizon, prévision directe).
variables_jour <- function(d) {
  ws <- (as.numeric(d$jour_semaine) - 1) / 7
  cbind(sin_semaine = sin(2 * pi * ws), cos_semaine = cos(2 * pi * ws),
        week_end = as.numeric(d$est_week_end),
        sin_an = sin(2 * pi * d$t / 365.25), cos_an = cos(2 * pi * d$t / 365.25),
        ferie = as.numeric(d$ferie_semaine), pont = as.numeric(d$pont),
        noel = as.numeric(d$type_vacances == "Noël"), ete = as.numeric(d$type_vacances == "Été"),
        vacances = as.numeric(d$type_vacances == "Autres vacances"),
        temperature = (d$temperature_moyenne - 14) / 7, pluie = log1p(d$pluie_mm) - 0.5,
        vent = (d$vent_moyen_m_s - 3) / 1.2)
}

tableaux_lstm <- function(z, V, origines, L = 28, H = 14) {
  n <- length(origines); f <- ncol(V)
  passe <- array(0, c(n, L, f + 1)); futur <- array(0, c(n, H * f)); cible <- array(NA_real_, c(n, H))
  for (i in seq_len(n)) {
    o <- origines[i]
    passe[i, , ] <- cbind(z[(o - L + 1):o], V[(o - L + 1):o, ])
    futur[i, ] <- as.numeric(t(V[o + seq_len(H), ]))
    if (o + H <= length(z)) cible[i, ] <- z[o + seq_len(H)]
  }
  list(passe = passe, futur = futur, cible = cible)
}

# ===========================================================================
# Ressemblance entre blocs : noyau de Nadaraya-Watson
# ===========================================================================
#   X̃ = Σ_i W_i X_i,   W_i ∝ k((B_courant - B_i) / w),   k noyau gaussien.
# Produit de noyaux gaussiens univariés = noyau gaussien de la distance : on
# utilise la distance quadratique moyenne entre blocs centrés (comparable quelle
# que soit la taille b). Tous les blocs passés comptent, pondérés par leur
# ressemblance ; w (fenêtre) et b (taille du bloc) se calibrent par validation croisée.
prevision_blocs_noyau <- function(z, o, H = 14, b = 14, w = 0.1) {
  courant <- z[(o - b + 1):o]; m0 <- mean(courant)
  fins <- seq(o - 7L, b, by = -7L)
  fins <- fins[fins + H <= o]
  dist <- vapply(fins, function(j) { bj <- z[(j - b + 1):j]
                                     sqrt(mean(((bj - mean(bj)) - (courant - m0))^2)) }, numeric(1))
  # exp(-d²/2w²) normalisé = exp(-(d² - d_min²)/2w²) normalisé : on retire d_min²
  # dans l'exposant pour éviter que tous les poids soient nuls en machine (0/0)
  W <- exp(-0.5 * (dist^2 - min(dist)^2) / w^2)
  W <- W / sum(W)
  suites <- vapply(fins, function(j) z[j + seq_len(H)] - mean(z[(j - b + 1):j]), numeric(H))
  list(prevision = m0 + as.numeric(suites %*% W), fins = fins, poids = W, distances = dist)
}
