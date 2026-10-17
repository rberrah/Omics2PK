# Rosuvastatine n40 : les cinq CSV originaux traités des auteurs sont référencés, accès automatisé bloqué

**Date : 17 octobre 2026.** Cohorte : Silveira et al., Frontiers in Pharmacology (2021), DOI https://doi.org/10.3389/fphar.2021.752960, n40 volontaires, authentique métabolomique prédose et rosuvastatine AUC/Cmax quantitatives.

## Nouveau chaînon original vérifié

L'article source cite explicitement le code original de l'auteur [GustavoHBDuarte/Metabolomics_PK_proj](https://github.com/GustavoHBDuarte/Metabolomics_PK_proj). Ce dépôt contient le notebook natif **`Metabolomics_Elastic_net.ipynb`**, à révision git authentique `dbf9182b0f72f8b290f623f7cc6560c5a5f99d73`.

Son README, vérifié dans le dépôt d'origine, précise que des fichiers CSV prétraités sont **liés à l'intérieur du notebook** : données de métabolomique ESI+/ESI−, examens biologiques et analyse combinée, issus initialement de ProgenesisQI puis normalisation standard interne/QC, racine cubique et sélection RSD QC. Les cinq références dans le code ne sont donc **pas** des données inventées, et leur origine est liée au modèle de l'article.

**Audit exécuté sur le fichier original de 21 Mo depuis GitHub Actions** : [run original scientifique](https://github.com/rberrah/Omics2PK/actions/runs/39024522995), module `scripts/audit_rosuvastatin_author_notebook_data_urls.py` ; SHA256 de la version du notebook source `25e361ef4e3329a8cc6e99336b177d614e5d8238d9741f5033d8e336dc1abecd`, **128 cellules dont 104 cellules code et 69 cellules avec sorties de calcul**. L'audit n'a **pas** lu/affiché ces sorties d'origine et n'a jamais journalisé les mesures biologiques individuelles.

- Exactement **5 appels `pd.read_csv()`**, utilisant **5 URL sur `drive.google.com`** dans le notebook d'origine. Les URL cloud contiennent des identifiants de fichiers opaques ; ils ne sont **pas** recopiés dans les rapports ou logs GitHub.
- Chaque tentative d'accès HTTPS par le script depuis GitHub Actions a retourné **une page d'accès/HTML Google Drive**, pas un CSV. Cela ne signifie **pas que les fichiers sont inexistants** ou que le partage est privé : un navigateur connecté ou l'accord de l'auteur peut suffire, mais cet audit ne les a **pas** obtenus.
- Par suite, **zéro ligne et zéro colonne de ces cinq CSV sources n'a pu être qualifiée** ; les identités exactes `volunteer_id ↔ metabolomics matrix ↔ AUC/Cmax` restent non démontrées par les fichiers accessibles à l'audit. Même un CSV prétraité récupéré devra être inspecté pour statut prédose, rôle ESI, présence de molécules post-dose, unités et correspondance origine AUC/Cmax.
- Les MAF officiels MetaboLights contiennent **69 et 73 colonnes d'échantillons originaux reconnues à l'identique**, mais **0 valeur numérique renseignée** dans les colonnes d'abondance de leurs 15 lignes d'annotation chacun ([rapport MAF](ROSUVASTATIN40_ORIGINAL_MAF_NO_NUMERIC_INTENSITY_20261017.md)). Les nombreux fichiers mzML du dépôt MTBLS3128 restent exploitables pour reconstruction par prétraitement MS.
- Aucun réentraînement ni performance E3 réelle n'est annoncé ici. Les R² très élevés du manuscrit restent des résultats des auteurs avec risque de fuite de sélection (plusieurs milliers de configurations), pas une preuve indépendante.

## Reproduction / demande scientifique minimale

1. Pour voir les cinq liens authentiques, ouvrir le [notebook original dans GitHub](https://github.com/GustavoHBDuarte/Metabolomics_PK_proj/blob/dbf9182b0f72f8b290f623f7cc6560c5a5f99d73/Metabolomics_Elastic_net.ipynb) ; les copier depuis les cellules `pd.read_csv` dans un navigateur autorisé, plutôt que depuis le rapport ou un miroir tiers. Les fichiers CSV doivent être vérifiés quant à leur partage et leurs droits de réutilisation.
2. Priorité : obtenir **les 5 CSV originaux** avec leurs colonnes de provenance et la clé appariant les 40 personnes à leurs AUC/Cmax mesurées. Si les fichiers sont indisponibles, demander à l'auteur le seul `40 participants x [metabolites pré-dose, AUC, Cmax, covariables]` dé-identifié, description des dates, unités et types de normalisation, avec une autorisation claire de recherche locale.
3. N'importer aucun profil patient ou ligne originale dans ce dépôt public ; local uniquement. Pour une validation honnête, nombre de sujets réels n40, folds groupés par individu, feature selection + normalisation + tuning dans chaque *train* et baseline clinique à mêmes 40 sujets.

**Statut : `5_AUTHOR_LINKED_ORIGINAL_PROCESSED_DATA_CSV_REFERENCES_FOUND; 5_GDRIVE_INTERSTITIAL; ORIGINAL_40_HUMAN_NUMERIC_PK_OMICS_CROSSWALK_NOT_YET_ACQUIRED; E3_ABSTAIN`.**

Si ces fichiers deviennent accessibles, cette cohorte est le **candidat P0** pour la première vraie validation E3, devant les n17 metformine/n14 Jarmusch puisque la source d'acquisition moléculaire brute et les scripts de traitement sont identifiés.
