# Rosuvastatine n40 — validation privée des 5 CSV originaux des auteurs

## Pourquoi

La cohorte indépendante Silveira et al. (*Front. Pharmacol.* 2021), DOI [10.3389/fphar.2021.752960](https://doi.org/10.3389/fphar.2021.752960), [MetaboLights MTBLS3128](https://www.ebi.ac.uk/metabolights/MTBLS3128) contient authentiquement 40 sujets, leurs métabolomes sériques avant rosuvastatine et les AUC/Cmax de rosuvastatine plasmatique mesurées. Le [notebook original de l'auteur](https://github.com/GustavoHBDuarte/Metabolomics_PK_proj/blob/main/Metabolomics_Elastic_net.ipynb) référence **5 fichiers CSV originaux** via Google Drive. Tous les téléchargements directs puis avec `gdown` ont été tentés de façon réelle dans GitHub Actions **sans succès (5/5)** ; cela indique une indisponibilité d'accès automatisé, pas une absence scientifique.

Dans les originaux publics MetaboLights, **142 lignes échantillons LC-MS**, réparties en assays ESI+ **69** et ESI− **73**, correspondent exactement aux ISA. Les tables d'annotation MAF ont **zéro cellule de signal numérique** dans leurs colonnes échantillons ; les vrais signaux sont dans les `mzML` déposés ou dans les cinq CSV traités des auteurs. AUC/Cmax individuels n'ont **pas** été récupérés comme table numérique authentique reliée par identifiants. Résultats présents dans [l'audit MTBLS3128](ROSUVASTATIN40_MTBLS3128_REAL_PREDOSE_OMICS_PK_SOURCE_AUDIT_20261017.md), [audit MAF](ROSUVASTATIN40_ORIGINAL_MAF_NO_NUMERIC_INTENSITY_20261017.md) et [audit gdown](ROSUVASTATIN40_AUTHOR_GDOWN_DOWNLOAD_ABSTENTION_AND_UPLOAD_20261017.md).

## Fichiers nécessaires — upload privé, pas dans GitHub

Télécharger dans le navigateur **les cinq CSV originaux**, à partir des cinq appels `pd.read_csv()` du [notebook original](https://github.com/GustavoHBDuarte/Metabolomics_PK_proj/blob/main/Metabolomics_Elastic_net.ipynb). Le nom des fichiers originaux reste libre, ne pas réordonner ni réécrire les IDs des lignes. **Au minimum :**

1. Un CSV `ESI+` de **métabolomique pré-dose traitée**, dont une ligne (ou une clé de jointure originalement fournie par l'auteur) correspond à un volontaire biologique. Idéalement récupérer également la matrice ESI−.
2. Le CSV source contenant des **AUC et Cmax quantitatives mesurées** pour les 40 mêmes volontaires, avec unités de l'article.
3. Si les IDs diffèrent entre données moléculaires et PK, la **table de correspondance originale de l'auteur**. **Ne jamais** relier les deux tableaux par position de ligne, nombre d'échantillons ou nom de fichier supposé.
4. Les autres CSV de biologie clinique/covariables et combinés pour comparer des modèles sur les mêmes volontaires ; la somme de cinq fichiers est donc préférable à un seul.

Les fichiers originaux ne doivent pas être commités dans ce dépôt public, et les données individuelles ne doivent pas apparaître dans les logs CI ou les pages du site Omics2PK. Dans ce projet, une analyse de cohorte originelle n40 reste **non vérifiable** tant que la jointure des vrais AUC/Cmax avec les vraies omiques pré-dose n'a pas été démontrée.

## Auditer localement les fichiers originaux sans afficher les mesures personnelles

Depuis la racine du dépôt :

```bash
python -m pip install -e .
python scripts/qualify_local_original_rosuvastatin_csvs.py \
  --source-file /CHEMIN/ORIGINAL_1.csv \
  --source-file /CHEMIN/ORIGINAL_2.csv \
  --source-file /CHEMIN/ORIGINAL_3.csv \
  --source-file /CHEMIN/ORIGINAL_4.csv \
  --source-file /CHEMIN/ORIGINAL_5.csv
```

**Cette première exécution ne prétend pas apparier les sujets**. Elle retourne uniquement SHA256 des cinq CSV, dimensions et présence apparente de colonnes AUC et Cmax, sans imprimer les noms de sujets, les en-têtes complets ou la moindre valeur de concentration.

**Certains CSV des auteurs peuvent déjà réunir leurs variables métabolomiques avec AUC/Cmax sur les mêmes 40 lignes.** Dans ce cas, aucune fusion de fichiers par ordre de ligne n'est nécessaire : inspecter les noms réels des colonnes et utiliser uniquement le mode **combined** ci-dessous pour une première qualification structurelle :

```bash
python scripts/qualify_local_original_rosuvastatin_csvs.py \\
  --combined /CHEMIN/CSV_ORIGINAL_AUTEUR.csv \\
  --combined-id-column NOM_EXACT_ID_SUJET \\
  --auc-column NOM_EXACT_AUC \\
  --cmax-column NOM_EXACT_CMAX
```

Le validateur écarte explicitement AUC et Cmax des variables moléculaires candidates, rejette les colonnes PK *post-dose* supplémentaires, confirme 40 lignes personnelles distinctes et **n'entraîne aucun modèle**. Cela ne garantit ni la provenance d'une vraie mesure clinique ni le caractère pré-dose de tous les signaux : attestation documentaire requise.

Ensuite, après contrôle *humain* des noms des colonnes (ne pas inventer les noms) et que les deux matrices contiennent une vraie clé auteur, exécuter :

```bash
python scripts/qualify_local_original_rosuvastatin_csvs.py \
  --omics /CHEMIN/ORIGINAL_METABOLOME_PREDOSE.csv \
  --pk /CHEMIN/ORIGINAL_PK_AUC_CMAX.csv \
  --omics-id-column NOM_EXACT_ID_ECHANTILLON_OU_SUJET \
  --pk-id-column NOM_EXACT_ID_SUJET_ORIGINAL \
  --auc-column NOM_EXACT_AUC \
  --cmax-column NOM_EXACT_CMAX \
  --expected-original-subjects 40
```

Si les IDs échantillon et participant diffèrent, ne pas inventer une correspondance. On peut ajouter `--original-crosswalk /CHEMIN/CORRESPONDANCE_ORIGINALE_AUTEUR.csv` contenant les colonnes `omics_sample_id`, `pk_subject_id`, après vérification qu'elle est authentique et relationnelle 1:1. Le CLI rejette le moindre ID dupliqué, omis ou inconnu.

**Le programme refuse** : (a) matrices annotées mais sans intensités numériques comme les MAF officiels, (b) cible PK AUC/Cmax/clairance réutilisée comme prédicteur, (c) comptabilisation des 142 acquisitions comme 142 donneurs, (d) imputations non documentées, (e) dossier de connexion Google Drive HTML traité comme CSV, (f) fichiers PK de taille autre que les 40 volontaires lorsque l'essai décrit n40.

Si une jointure structurelle exacte de 40 personnes est trouvée, **aucun score de modèle n'est encore produit**. Le logiciel imprime un résultat au statut `STRICT_ORIGINAL_40_PERSON_STRUCTURE_JOIN_ONLY` et force quatre validations manuelles toujours non vérifiées : provenance originale, temps biologique *pré-dose*, unité/étalonnage AUC et Cmax, droits d'utilisation. Ensuite seulement, à données réelles accessibles, effectuer ré-estimation nested LODO, prétraitement dans les plis et comparateurs non-omiques. Il ne faut pas reproduire la recherche de 4800 combinaisons du manuscrit en choisissant les meilleures sur les plis test.

## Etat actuel

- [x] Publication site Omics2PK et intégration Ecosystem PMx / MultiOmics, CI verte.
- [x] MAF authentiques vérifiés ; pas d'intensité numérique traitée dans ces tables.
- [x] Notebook original et cinq liens Drive identifiés ; essais urllib et gdown infructueux.
- [x] **Validateur local privé et sans fuite** permettant l'import de vrais CSV originaux et le contrôle exact d'identifiants.
- [ ] **Réception par upload privé des cinq CSV** ou accès auteur à leurs liens.
- [ ] Jointure indépendante 40 personnes × intensités omiques pré-dose × AUC/Cmax quantitatives mesurées.
- [ ] Validation E3 et comparaison avec baseline clinique. **Aucun résultat clinique E3 indépendant n'est encore établi.**
