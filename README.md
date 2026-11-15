# Projet - Le vélo à Paris : comprendre et prévoir le trafic cyclable

Projet de séries chronologiques réalisé en R.

## Problématique

> Quelle part du trafic cyclable parisien s'explique par le calendrier, la météo et les événements,
> et peut-on le prévoir à quelques jours ?

## Données

| Source | Jeu de données | Contenu | Rôle |
|---|---|---|---|
| opendata.paris.fr | `comptage-velo-historique-donnees-compteurs` | Archives annuelles des comptages horaires (depuis 2018) | Série principale |
| opendata.paris.fr | `comptage-velo-donnees-compteurs` | Comptages horaires des 13 derniers mois | Série principale (fin) |
| opendata.paris.fr | `comptage-velo-compteurs` | Liste et localisation des compteurs | Choix du panel |
| meteo.data.gouv.fr | Station Paris-Montsouris | Température, pluie, vent (journaliers) | Exogènes |
| data.gouv.fr | Jours fériés, vacances scolaires (zone C) | Calendrier | Exogènes |

Points de vigilance :
- Garder un **panel de compteurs présents sur toute la période** (sinon fausse tendance).
- Repérer les **pannes** (longues suites de zéros, trous) avant d'agréger.
- Événements à dater : grèves de décembre 2019, confinements 2020.

## Plan et scripts

Tous les scripts se lancent depuis `projet_velo/`, dans l'ordre. `R/00_outils.R` (thème, régresseurs, diagnostic des résidus, ACF) et `R/00_previsions.R` (blocs, LSTM) sont chargés par les autres.

| Partie | Script | Contenu |
|---|---|---|
| 1. Données | `R/01_import.R`, `R/02_nettoyage.R` | Import, panel de 27 compteurs, pannes, météo, calendrier |
| 2. Exploration | `R/03_exploration.R` | Série, log, saisons semaine et année, ACF, périodogramme, météo, calendrier, événements |
| 3. Décomposition | `R/04_decomposition.R` | Moyennes mobiles (gain), `decompose`, méthode par blocs, MSTL, Fourier, Prophet, filtre quadratique, moyenne mobile modifiée (52 semaines) |
| 4. Ruptures | `R/05_ruptures.R` | Rupture par moindres carrés, Bai-Perron, CUSUM, analyse d'intervention |
| 5. Stationnarité | `R/06_stationnarite.R` | ADF (dont les 3 types d'aTSA), PP, KPSS ; (I − B), (I − B⁷) |
| 6. SARIMA | `R/07_sarima.R` | ACF/PACF, Yule-Walker, grille AIC/BIC, z-tests, racines, Wold (ψ_k) et intervalles, diagnostic des résidus, sur-ajustement, correction log → exp, ARCH/GARCH (tseries, rugarch), test BDS, (d, D) = (1, 1), simulation |
| 7. SARIMAX | `R/08_sarimax.R` | Semaine + Fourier + calendrier + météo, effets en %, rapports de vraisemblance, CCF, météo de la veille, part expliquée (Shapley) |
| 8. Concurrents | `R/09_concurrents.R` | Naïf, lissage simple, Holt, Holt-Winters, TBATS, Prophet (+ validation croisée), blocs à noyau calibrés, NNAR, LSTM (keras) |
| 9. Évaluation | `R/10_evaluation.R` | 14 modèles, origine glissante sur un an, h = 1 à 14, REQM/MAPE, couverture, Diebold-Mariano, météo observée / de la veille / climatologique, prévision des 14 prochains jours |
| 10. Rendu | `rapport/rapport.qmd`, `rapport/soutenance.qmd` | Rapport et diaporama Quarto |

Chaque script écrit ses figures dans `figures/<numéro>_<partie>/` et ses tableaux dans `resultats/`.

**Environnement.** R 4.4 avec `data.table`, `ggplot2`, `patchwork`, `forecast`, `tseries`, `aTSA`, `strucchange`, `sandwich`, `rugarch`, `prophet`, `keras3`. Le LSTM utilise l'environnement Python `r-keras` (TensorFlow), créé par `reticulate::virtualenv_create("r-keras")` puis `reticulate::virtualenv_install("r-keras", c("tensorflow-cpu", "keras"))`.

## Reproduire les résultats

1. Installer R (≥ 4.4) et les packages : `data.table`, `ggplot2`, `patchwork`, `scales`, `forecast`, `tseries`, `aTSA`, `strucchange`, `sandwich`, `rugarch`, `prophet`, `keras3`, `ggrepel`, `ggridges`, `zoo`, `stringi`, `jsonlite`.
2. Pour le LSTM uniquement : créer l'environnement Python `r-keras` (voir « Environnement » plus haut).
3. Depuis le dossier du projet, lancer les scripts dans l'ordre :
   ```r
   for (s in sort(list.files("R", pattern = "^(0[1-9]|10)_.*\\.R$", full.names = TRUE))) source(s, echo = FALSE)
   ```
   `R/01_import.R` télécharge les données publiques (environ 140 Mo, quelques minutes) dans `data/1_brut/`, qui n'est pas versionné. La série finale `data/3_final/` est fournie : les scripts 03 à 10 peuvent tourner sans téléchargement.
4. Générer le rapport : `quarto render rapport/rapport.qmd` et `quarto render rapport/soutenance.qmd`.

## Sources et licences des données

| Données | Producteur | Licence |
|---|---|---|
| Comptages vélo (historique et 13 derniers mois) | Ville de Paris - opendata.paris.fr | [ODbL](https://opendatacommons.org/licenses/odbl/) |
| Météo quotidienne, station Paris-Montsouris | Météo-France - meteo.data.gouv.fr | [Licence Ouverte 2.0](https://www.etalab.gouv.fr/licence-ouverte-open-licence/) |
| Jours fériés | calendrier.api.gouv.fr | Licence Ouverte 2.0 |
| Vacances scolaires | Ministère de l'Éducation nationale - data.education.gouv.fr | Licence Ouverte 2.0 |

Les fichiers de `data/3_final/` sont dérivés des comptages de la Ville de Paris : ils sont diffusés sous la même licence **ODbL**. Toute réutilisation doit citer les producteurs ci-dessus.

## Questions à discuter

- Météo future inconnue au moment de prévoir : météo observée vs météo de la veille.
- Comptages = données de type Poisson : justification du log et limites de l'hypothèse gaussienne.

## Arborescence

```
projet_velo/
├── R/          scripts numérotés, à exécuter dans l'ordre depuis projet_velo/
├── data/       voir data/LISEZMOI.md
│   ├── 1_brut/            fichiers téléchargés, jamais modifiés
│   ├── 2_intermediaire/   1 fichier par source, lisible, noms de colonnes en français
│   └── 3_final/           série journalière prête pour la modélisation
├── figures/        un sous-dossier par script
├── resultats/      tableaux (csv) et modèles estimés (modeles/)
└── rapport/        rapport et soutenance (Quarto)
```

## Avancement

- [x] 1. Données - `R/01_import.R`, `R/02_nettoyage.R` (16 figures)
- [x] 2. Exploration - `R/03_exploration.R` (16 figures)
- [x] 3. Décomposition - `R/04_decomposition.R` (9 figures)
- [x] 4. Ruptures - `R/05_ruptures.R` (5 figures)
- [x] 5. Stationnarité - `R/06_stationnarite.R` (4 figures)
- [x] 6. SARIMA - `R/07_sarima.R` (16 figures)
- [x] 7. SARIMAX - `R/08_sarimax.R` (10 figures)
- [x] 8. Concurrents - `R/09_concurrents.R` (8 figures)
- [x] 9. Évaluation - `R/10_evaluation.R` (8 figures)
- [x] 10. Rapport (`rapport/rapport.html`) et soutenance (`rapport/soutenance.html`)

Total : 92 figures. 

## Résultats principaux

- Part de la variance de log(vélos) expliquée (Shapley) : semaine 28 %, saison annuelle 22 %, calendrier 20 %, météo 13 %, tendance 3 %, dynamique ARMA 2 % ; non expliqué 12 %.
- SARIMA(3,0,3)(1,1,2)₇ : racine unité du MA saisonnier, donc semaine déterministe, d'où un SARIMAX(1,0,1)(1,0,1)₇ avec indicatrices de jour, Fourier K = 8, calendrier et météo.
- Évaluation sur un an (14 modèles, 352 origines, h = 1 à 14) : SARIMAX meilleur à 1 jour (MAPE 9,8 %), Prophet meilleur au-delà (≈ 11 %), naïf saisonnier 23 à 30 %. Sans prévision météo, le SARIMAX monte à environ 15 % ; avec la météo de la veille, 12,4 % à 1 jour.
- Résidus du SARIMA : bruit blanc faible mais pas fort (test BDS), effet ARCH ; (d, D) = (1, 1) donne des racines unité dans le MA (double sur-différenciation).
