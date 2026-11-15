# Le dossier `data/`

Les données passent par trois niveaux. Chaque niveau est produit par un script, on ne modifie jamais un fichier à la main.

```
data/
├── 1_brut/            fichiers téléchargés tels quels          (créés par R/01_import.R)
│   ├── velo/
│   ├── meteo/
│   └── calendrier/
├── 2_intermediaire/   fichiers lisibles, même format, 1 ligne = 1 jour  (créés par R/01_import.R)
└── 3_final/           série prête pour la modélisation          (créée par R/02_nettoyage.R)
```

Pour tout reconstruire : supprimer `2_intermediaire/` et `3_final/`, puis relancer les scripts depuis `projet_velo/`.
Supprimer `1_brut/` force un nouveau téléchargement (environ 30 Mo, 3 minutes).

---

## 1. `1_brut/` - les sources, sans aucune modification

| Fichier | Contenu | Source |
|---|---|---|
| `velo/velo_archive_2018.zip` … `velo_archive_2025.zip` | Comptages **horaires** de tous les compteurs vélo, une archive par année | Ville de Paris, opendata.paris.fr, jeu « Comptage vélo - Historique » |
| `velo/velo_13_derniers_mois_par_jour.csv` | Comptages des 13 derniers mois, **déjà additionnés par jour** | Ville de Paris, jeu « Comptage vélo - Données compteurs » |
| `meteo/meteo_paris_1950-2024.csv.gz` | Météo quotidienne de toutes les stations de Paris, 1950-2024 | Météo-France, meteo.data.gouv.fr |
| `meteo/meteo_paris_2025-2026.csv.gz` | Idem, 2025-2026 | Météo-France |
| `calendrier/jours_feries.json` | Jours fériés en métropole | calendrier.api.gouv.fr |
| `calendrier/vacances_scolaires_paris.csv` | Vacances scolaires de l'académie de Paris | Ministère de l'Éducation, data.education.gouv.fr |

Pièges rencontrés dans ces fichiers (corrigés à l'import) :
- Le format change d'une année à l'autre : séparateur `;` ou `,`, comptage écrit `12` ou `12.0`, deux formats de date.
- En 2025, les comptages ≥ 1000 sont écrits avec une espace (« 1 022 ») : 6 694 valeurs corrigées.
- L'archive 2025 fait 2 Go une fois décompressée, à cause des liens vers les photos.

---

## 2. `2_intermediaire/` - les fichiers de travail

### `velo_par_compteur_et_par_jour.csv`
Une ligne = **un compteur, un jour**. 375 688 lignes, 2018-01-01 → 2026-11-13.

| Colonne | Signification |
|---|---|
| `date` | Jour (heure de Paris) |
| `id_compteur` | Identifiant du compteur. Attention : peut changer d'une année à l'autre pour un même compteur |
| `nom_compteur` | Adresse et sens de circulation, ex. « 97 avenue Denfert Rochereau SO-NE » |
| `id_site` | Identifiant du lieu de comptage (un site peut avoir plusieurs compteurs, un par sens) |
| `nb_velos` | Nombre de vélos comptés dans la journée |
| `nb_heures_mesurees` | Nombre d'heures avec une mesure (24 = journée complète) |

**Attention : ce fichier ne doit pas être additionné tel quel.** Problèmes à traiter dans `R/02_nettoyage.R` :
- **Le nombre de compteurs augmente** : 14 en 2018, environ 45 en 2019, 90 en 2020, 170 ensuite. La somme brute monterait mécaniquement.
- **Doublons** : en 2021-2022, il y a 167 identifiants pour 102 noms, ce qui laisse penser qu'un même compteur apparaît deux fois.
- **Ligne sans identifiant** en 2024-2025 : une série complète sans identifiant, sans nom et sans coordonnées (environ 30 % du trafic). Elle sera écartée.
- Le dernier jour (2026-11-13) est incomplet.

### `meteo_paris_par_jour.csv`
Une ligne = un jour, station **Paris-Montsouris**. 3 237 jours, 2018-01-01 → 2026-11-11.

| Colonne | Signification |
|---|---|
| `date` | Jour |
| `pluie_mm` | Hauteur de pluie (mm) |
| `temperature_min`, `temperature_max`, `temperature_moyenne` | Températures du jour (°C) |
| `vent_moyen_m_s` | Vitesse moyenne du vent (m/s) |

### `jours_feries.csv`
| Colonne | Signification |
|---|---|
| `date` | Jour férié |
| `nom_jour_ferie` | Ex. « Lundi de Pâques » |

### `vacances_scolaires.csv`
Une ligne = une période de vacances, académie de Paris (zone C).

| Colonne | Signification |
|---|---|
| `nom_vacances` | Ex. « Vacances de la Toussaint » |
| `annee_scolaire` | Ex. « 2017-2018 » |
| `premier_jour` | Premier jour de vacances |
| `jour_de_rentree` | Jour de la reprise des cours (n'est **pas** un jour de vacances) |

---

## 3. `3_final/` - la série à modéliser (créée par `R/02_nettoyage.R`)

### `velo_paris_par_jour.csv` ← **le fichier à utiliser pour la suite**
Une ligne = un jour. 2 872 jours, 2019-01-01 → 2026-11-11.
Trafic = somme des **27 compteurs du panel**, présents et fiables sur toute la période.

| Colonne | Signification |
|---|---|
| `date` | Jour |
| `velos` | **Nombre de vélos comptés dans la journée par les 27 compteurs du panel** (la série à modéliser) |
| `jour_semaine` | lundi … dimanche |
| `est_week_end` | VRAI le samedi et le dimanche |
| `est_ferie`, `nom_jour_ferie` | Jour férié et son nom |
| `est_vacances`, `nom_vacances` | Vacances scolaires de Paris (zone C) et leur nom |
| `evenement` | Grève des transports, Confinement 1/2/3, Jeux olympiques (vide sinon) |
| `temperature_moyenne`, `temperature_min`, `temperature_max` | Températures à Paris-Montsouris (°C) |
| `pluie_mm` | Pluie du jour (mm) |
| `vent_moyen_m_s` | Vent moyen (m/s) |
| `meteo_interpolee` | VRAI si la météo du jour a été reconstituée (manquante à la source) |
| `part_estimee` | Part du trafic du jour estimée car un compteur était en panne (moyenne 1,9 %) |
| `compteurs_valides` | Nombre de compteurs du panel qui fonctionnaient ce jour-là (sur 27) |
| `coupure_generale` | VRAI les 6 jours où tout le système de comptage était coupé (trafic entièrement estimé) |

Les règles de passage des données brutes à cette série sont expliquées et illustrées par les figures de `figures/02_nettoyage/`.

### `panel_compteurs.csv`
Une ligne par compteur (161) : retenu ou non, et pourquoi.

| Colonne | Signification |
|---|---|
| `libelle` | Nom lisible du compteur |
| `compteur` | Nom normalisé (minuscules, sans accents), clé utilisée par le script |
| `retenu` | VRAI s'il fait partie du panel |
| `raison` | « retenu », « installé trop tard », « arrêté avant la fin », « x % de jours valides », « niveau qui dérive » |
| `premier_jour` | Premier jour de données |
| `dernier_valide` | Dernier jour où il fonctionnait correctement |
| `part_valide` | Part des jours valides depuis 2019 |
| `derive_max` | Écart maximal de son niveau annuel face aux autres compteurs (1 = aucun écart) |
| `niveau_moyen` | Vélos par jour en moyenne |
| `dans_panel_elargi` | VRAI s'il fait aussi partie du panel de contrôle (depuis 2020) |

### `jours_corriges.csv`
Une ligne par valeur remplacée dans le panel (2 072 jours compteur, 2,7 %).

| Colonne | Signification |
|---|---|
| `date`, `libelle`, `compteur` | Jour et compteur concernés |
| `motif` | jour absent, journée incomplète, zéro, anormalement bas, anormalement haut |
| `valeur_observee` | Ce que le compteur a enregistré (vide si absent) |
| `valeur_retenue` | La valeur estimée qui la remplace |
