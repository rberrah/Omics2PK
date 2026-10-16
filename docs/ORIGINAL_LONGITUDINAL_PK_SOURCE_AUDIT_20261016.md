# Qualification des premières données PK longitudinales — 16 octobre 2026

**Statut : données de concentration–temps réellement observées ingérées et contrôlées, PK_LONG_V1 validé ; aucune PopPK ajustée ni omique appariée à ces observations pour le moment.**

## Theoph — 12 sujets, observations réelles

Source de référence : R datasets::Theoph, d’après l’étude pharmacocinétique orale de la théophylline ; exemple compilé par Boeckmann, Sheiner et Beal (1994). Source CSV publiquement accessible dans Rdatasets, commit épinglé cc03c29690889dbe83089f5206e2422db8c3f71f, SHA-256 original f9b77033090644cf38c83056c25c8137bca71cf629d7c0159a72ef18095726b0. Données brutes non commitées dans le dépôt.

Le pipeline scientifique a analysé sur les données originales **12 vrais profils individuels** comprenant chacun **une dose orale et 11 concentrations sériques**, soit 12 doses + 132 mesures PK, et non 132 patients. Les doses initialement en mg/kg sont converties en mg à partir des poids réellement décrits. Les 12 prélèvements t=0 sont marqués predose ; chez 3 sujets, la concentration originale t=0 est non nulle, **conservée**, non remplacée par zéro.

Calculs descriptifs d’exposition par la règle des trapèzes linéaires, pour les concentrations observées jusqu’au dernier point (entre 23 et 26 h) :

| Résultat agrégé | Médiane originale observée |
|---|---:|
| Cmax | 8,465 mg/L |
| Tmax | 1,135 h |
| AUC(0,dernier temps observé) | 95,40665 mg·h/L |

Cela ne détermine PAS AUC(0,infini), clairance/F, demi-vie terminale fiable ni performance PopPK. Pas d’intervalle de confiance déclaré sans bootstrap patient rigoureux. Aucun RNA/protéome/métabolome véritable dans cette série.

Reproduction : commande Python -m scripts.qualify_public_theoph_pk. La source est réaccédée par URL figée, son SHA vérifié, le contrat contrôlé puis seuls effectifs et résultats agrégés imprimés ; absence d’export des observations individuelles.

## Inventaire indépendant des jeux originaux nlmixr2data

Sources d’origine : dépôt nlmixr2/nlmixr2data, commit f2cfb01ade88d4d30e66f56b19efce2a2c6e26bd, quatre objets RDA originaux lus de manière temporaire ; jamais de lignes PK originales publiées dans les logs. Les empreintes SHA-256 exactes sont contrôlables dans la sortie CI.

| Original RDA | Sujets | Lignes | Événements EVID | Interprétation |
|---|---:|---:|---|---|
| theo_sd | 12 | 144 | 12 doses code 101 ; 132 DV code 0 | Observations PK originales au jour 1 |
| theo_md | 12 | 348 | 84 doses code 101 ; 264 DV code 0 | Après phase initiale, **simulations** de schémas quotidiens et concentrations |
| warfarin | 32 | **515** | 32 doses code 1 ; 483 lignes code 0 | DV combinent concentrations PK et activité prothrombinique PD (DVID distinct) |
| mavoglurant | **120** | 2 678 | 198 doses code 1 ; 2 480 DV code 0 | Deux occasions OCC, 118 temps distincts, données PK cliniques, aucun RNA |

**Anomalie documentaire conservée :** la page de référence de nlmixr2data donne 519 lignes pour warfarin, mais l’objet RDA réellement épinglé en possède 515. On ne reconstruit pas les quatre lignes absentes. L’explication/version doit être contrôlée avant validation d’estimation PK/PD.

Pour mavoglurant, **une occasion OCC n’est pas automatiquement un épisode de doses répétées en temps global** ; il faut qualifier l’origine du temps, l’ordre des visites, RATE et les prélèvements par sujet/occasion. Les 198 doses sur 120 individus ne représentent pas, à elles seules, un schéma longitudinal continu complet sur 48 heures.

## Couplage futur omiques → PK

Le nouveau contrat vérifie séparément les événements de dose, concentrations, BLQ et historique de traitement et empêche toute assimilation d’une simulation à une observation. Les échantillons omiques horodatés (RNA/protéine/métabolite) seront disponibles comme covariables dynamiques uniquement **avant le temps de prédiction (landmark)**, avec fenêtre de validité définie avant d’analyser les résultats.

Un état omique mesuré 48 h après une dose de paracétamol ne peut pas servir à « prédire » sa Cmax à 1 h. L’étude Fannin et al. (2010), DOI 10.1002/hep.23330, documente un protocole avec PK sérique précoce et RNA/métabolomique répétés ultérieurement chez les mêmes volontaires ; son format individuel téléchargeable et l’appariement PK–RNA restent à qualifier. Candidat méthodologique important pour étudier le lien PK→réponse biologique, puis potentiellement une prédiction prospective ultérieure, mais pas démonstration actuelle.

Le prochain niveau **scientifique** sera un modèle PK identifiable sur vrais C(t) et doses, avec diagnostic residuals/FOCEi/SAEM dans le pilier B, puis comparaison clinico-PGx vs omiques sur des échantillons authentiquement appariés et à la date adéquate.

Aucune de ces données ne justifie encore un calculateur de dose patient ni une publication affirmant une amélioration de la PopPK par les omiques longitudinales.


## Actualisation cruciale mavoglurant : occasions multiples ≠ doses répétées continues

L'objet original contient exactement **120 sujets**, **198 couples sujet×OCC**, **198 événements EVID=1** et **2 480 lignes de concentrations EVID=0**, chaque couple sujet×OCC n'ayant **qu'une seule dose enregistrée**. **78 sujets** ont deux occasions. Le protocole original A2121 utilisait des perfusions IV de dix minutes, d'après Wendling et la vignette nlmixr2 ; concentrations publiées en ng/mL (https://nlmixr2.org/articles/mavoglurant.html). Ce jeu fournit donc une **véritable PK humaine longitudinale multi-occasion**, utile pour quantifier l'évolution entre périodes et valider les effets de covariables dans une PopPK, mais PAS un historique de prises répétées en continu. Les temps sont à interpréter dans leur occasion ; l'intervalle réel entre les deux occasions n'est pas présent comme calendrier absolu.
