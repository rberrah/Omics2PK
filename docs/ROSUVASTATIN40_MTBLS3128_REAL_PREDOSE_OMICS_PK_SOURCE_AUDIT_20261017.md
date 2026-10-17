# P0 E3 — rosuvastatine n40, véritable métabolomique prédose et PK quantitative : qualification des sources originales

**Date 17 octobre 2026. Étude de plus forte pertinence trouvée pour Omics2PK.**

**Article primaire** : Silveira AMR et al. (2021), *Frontiers in Pharmacology*, DOI https://doi.org/10.3389/fphar.2021.752960 ; [PMC8633954](https://pmc.ncbi.nlm.nih.gov/articles/PMC8633954/). Quarante volontaires sains ont eu leur profil **métabolomique sérique pré-dose** mesuré par LC-MS (+/−), avant quantification ciblée de la rosuvastatine plasmatique et calcul des **AUC et Cmax**. Article avec origine séparée de la mesure PK, pas une intensité m/z relative confondue avec une concentration ng/mL. Vrai essai humain, différent des n17 metformine et n14 Jarmusch.

**Données métabolomiques originales identifiées sans ambiguïté** : [MetaboLights MTBLS3128](https://www.ebi.ac.uk/metabolights/MTBLS3128), également référencé par [OmicsDI MTBLS3128](https://www.omicsdi.org/dataset/metabolights_dataset/MTBLS3128), qui recense environ **1 833 fichiers**, incluant de nombreux `mzML` annotés prédose. C'est le dépôt expressément déclaré par les auteurs.

## Accès aux fichiers sources réel (CI GitHub)

**Validation exécutée** : [qualification source d'article/MetaboLights](https://github.com/rberrah/Omics2PK/actions/runs/39023462155) et [correspondance des tables ISA originales](https://github.com/rberrah/Omics2PK/actions/runs/39023667343). Script `scripts/audit_rosuvastatin40_original_mtbls3128.py` et `scripts/audit_rosuvastatin40_isa_patient_crosswalk.py`. Aucun nom de volontaire, matrice de valeurs ou C(t) originel(le) n'a été exporté en CI ou sur GitHub.

| Origine véritable | Résultat confirmé | Conclusion |
|---|---|---|
| MetaboLights API originale `/studies/MTBLS3128` | Métadonnée officielle récupérée et identité MTBLS3128 contrôlée. SHA256 **`1b9d453d9ff668e5024292f06fede36518613fbd5a9ba7cceaa7ca7301bfaaa3`** | 2 assays LC-MS déclarés |
| MetaboLights `/files?include_raw_data=false` et téléchargement des ISA `download?file=` | **6 fichiers ISA originaux** inventoriés et inspectés, 1 source des échantillons, 2 fichiers d'assay + fichiers MAF/investigation. **142 lignes distinctes** dans la table des échantillons et deux assays de **69** et **73** lignes, dont les identifiants d'échantillon correspondent **exactement à 69 et 73** des lignes de la table mère | **Appariement original échantillon↔technique LC-MS vérifié**, mais 142 analyses/échantillons ne signifient PAS 142 personnes ; article n40 volontaires, avec contrôles QC, modes ionisation et possibles répétitions à qualifier |
| Métadonnées ISA réelles des échantillons, des assays et des MAF | Les colonnes inspectées ne comportent pas de champ explicitement quantifié **AUC, Cmax, concentration du médicament** ni d'identifiant clinique de participant associé au résultat PK | **Original PK outcome matrix NOT FOUND dans ces ISA** |
| Europe PMC `PMC8633954/supplementaryFiles` | Original article package d'empreinte **`643cb780de4213f0d6cdc29176ce31b53acbe7804b5041294d18fc56c008362b`**, 844 751 octets. Trois XLSX originaux `Table1.XLSX` (2 feuilles PT_BR/ENG ; 46 et 50 lignes), `Table2.XLSX` (72 lignes), `Table3.XLSX` (52 lignes), plus 1 PDF `DataSheet1.PDF` | Ces feuilles sont du matériel de publication, leurs en-têtes initiaux examinés n'affichent pas d'AUC ou Cmax individuelle ; ce constat n'exclut pas des valeurs PK enfouies plus loin, mais **aucune matrice PK×personne validée** n'est démontrée |

**Lacune binaire unique pour E3 :** récupérer les valeurs **AUC et Cmax (ou dose-temps-C(t) d'origine) pour les 40 vrais volontaires**, puis valider les correspondances exactes avec le préfixe/ID prédose des 142 échantillons MetaboLights. Une présence de n40 dans les statistiques de l'article ou des points sur Figures 2/3 ne suffit pas à reconstruire un vrai `sample_id↔subject_id↔AUC/Cmax`.

## Évaluation scientifique à faire seulement après acquisition PK originale

- Contrôles de provenance : données PK quantitatives, mesures pré-dose (ne jamais utiliser les caractéristiques post-dose), affectation positive/négative et contrôles blancs, doublons/QC, unités et statut BLQ.
- Définir **40 sujets biologiques** réellement indépendants ; une validation leave-one-*sample* n'est pas une validation leave-one-*patient*. Tous les échantillons/analyses du même individu restent dans le même pli.
- Comparateurs cliniques indispensables : moyenne/âge/sexe/poids/IMC/biologie, avec/sans métabolome et mêmes sujets complets. Sélection des métabolites, correction de batch et transformations ajustées **uniquement dans chaque entraînement**.
- **Éviter l'optimisme** : l'article annonce des R² LOOCV d'AUC 0,93 en ESI+ et 0,86 en ESI−, mais décrit une recherche de 4 800 configurations et diverses étapes de normalisation en amont ; ces chiffres sont **rapportés par les auteurs**, pas reproduits indépendamment, et doivent être requalifiés dans une CV intégralement imbriquée si l'optimisation a utilisé les données de validation.
- Sorties primaires `AUC/Cmax` en vraie unité PK et MAE/RMSE/biais, calibration, erreur de prédiction relative, bootstrap apparié des 40 personnes et intervalles exploratoires ; aucune recommandation de dose.
- Ce n'est pas directement un jumeau PBPK ni un SAEM PKpop dynamique si seules les AUC/Cmax sont disponibles, mais c'est un **E3 clinique métabolome→PK** crédible sur une cohorte indépendante.

**Statut actuel : `ORIGINAL_PREDOSE_OMICS_AND_ASSAY_IDS_VERIFIED; ORIGINAL_PATIENT_QUANTITATIVE_PK_CROSSWALK_NOT_VERIFIED`.** Il serait trompeur de comparer M1bis/M1/M2/M3/M4 avant le rapprochement des vraies sorties quantitatives.
