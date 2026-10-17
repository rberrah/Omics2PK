# Jarmusch et al. 2020 — audit des matrices individuelles PK ciblée × métabolomique (17 octobre 2026)

**Décision : E3 — UNVERIFIED PAIRING / ABSTAIN.** Les observations individuelles de concentrations étalonnées ne sont pas suffisamment accessibles pour joindre la matrice d'intensités métabolomiques originales. **Ne pas entraîner, afficher ni revendiquer un score Omics→PK clinique sur cette cohorte avant de pouvoir joindre les vraies mesures aux mêmes sujets.**

## Étude authentique et licences

- Étude primaire : Jarmusch AK et al., *Clinical and Translational Science* 2020, DOI [10.1111/cts.12785](https://doi.org/10.1111/cts.12785), PMCID [PMC7485946](https://pmc.ncbi.nlm.nih.gov/articles/PMC7485946/). Le protocole associe chez les mêmes **14 volontaires** cocktails de sondes CYP, pharmacocinétique **ciblée LC-MS/MS** (caféine, midazolam, oméprazole et ses métabolites), métabolomique MS non ciblée, microbiote 16S et génotypage CYP.
- **Licence du texte source** lue dans la métadonnée officielle PMC OA : **CC BY-NC** ; ne pas republier le contenu des tableaux/suppléments sans vérification de la licence adéquate.
- **Archive omique MassIVE** [MSV000082493](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=b8ffc66bab8d435aa863a0fed8fff6fa) indiquée comme **Partial Public** avec licence **CC0 1.0** du dépôt. Le dépôt affiche « urine, fecal and skin » en résumé ; ne pas en déduire sans inventaire que les échantillons **plasmatiques** nécessaires à l'analyse originale sont tous accessibles et joignables.
- Code original : https://github.com/alan-jarmusch/UntargetedMS-DrugMetabolism , son script `Code/20191030_UntargetedMS_DrugMetabolism.R` lit `GNPS_quant.csv` et `Metadata.txt` ; la consultation du code n'a pas établi la présence d'une matrice clinique source calibrée avec même identifiant. Un fichier omique relatif `quant.csv` ≠ des concentrations ciblées ng/mL vérifiées.

## Accès réellement contrôlés depuis GitHub Actions

**Exécution réelle réussie sur la branche PR #46** : [action source originale](https://github.com/rberrah/Omics2PK/actions/runs/39021336921).

Le script `scripts/audit_jarmusch14_original_omics_targeted_pk.py` a testé **les vraies sources**, ne journalisant que les métadonnées et le statut d'accès ; les tests synthétiques servent à prévenir les fuites de lignes individuelles.

| Origine | Résultat réel du 17/10/2026 | Conséquence |
|---|---|---|
| PMC OA officiel `PMC7485946.1` | Métadonnées d'article obtenues, licence **CC BY-NC**, **0** médias DOCX/XLSX originaux listés | Le package OA testé ne fournit pas le tableau PK clinique individuel |
| Wiley `cts12785-sup-0001-Supinfo.docx` | HTTP refusé (erreur HTTP) | Pas de récupération des méthodes supplémentaires |
| Wiley `cts12785-sup-0002-TableS2.xlsx` | HTTP refusé (erreur HTTP) | Impossible de confirmer si S2 contient des valeurs individuelles ou des données agrégées |
| MassIVE `MSV000082493` | Existence du dépôt, statut **Partial Public**, licence CC0 et format de données documentés | Inventaire précis plasmatiques × volontaires × temps, **pas encore numériquement qualifié** |
| GitHub original `alan-jarmusch/UntargetedMS-DrugMetabolism` | Script de traitement présent, références aux fichiers `GNPS_quant.csv`/`Metadata.txt` | Code analytique, **pas** matrice ciblée PK patient vérifiée |

**En particulier, HTTP 403 n'est pas une démonstration de non-existence ou d'interdiction juridique :** c'est une indisponibilité du document original dans les accès automatisés essayés. La même étude reste une très bonne cible pour une qualification par téléchargement autorisé ou correspondance avec l'auteur.

## Critères d'une véritable matrice E3 jointe

Pour chaque volontaire, exiger :
`volunteer_id, PK_analyte, dose_mg, dose_time, PK_sample_time, quantified_concentration, concentration_unit, quantified_is_BLQ`, et un second objet
`volunteer_id, untargeted_sample_file, sample_time_relative_to_dose, plasma_or_other_biofluid, feature_mz_rt, intensity, batch, normalization`. Exiger en plus la cohérence des noms d'échantillon avec les fichiers originaux, les références et contrôles qualité, et une preuve que l'omique candidate a été prélevée **avant l'administration** si elle sert à prédire l'exposition future.

Le simple fait qu'il existe 14 participants dans l'article et 14 participants dans MassIVE ne démontre **jamais** qu'une ligne omique et une concentration concernent le même volontaire. Il est interdit de joindre par ordre de lignes ou par proximité de nombre d'échantillons. Les ratios d'activité CYP publiés au niveau du groupe ne remplacent pas les profils individuels.

## Priorité relative

1. **Jarmusch n14** : original ciblé PK et MS non ciblée existent dans protocole, MassIVE public partiel ; il manque la clé sujet→quantité PK ciblée. Seule source nouvelle potentiellement réellement riche en 16S+PGx+PK dans le même protocole.
2. **Tee metformine n17** : les **255 concentrations** LC-MS et paramètres individuels sont décrits en Table S5.1/S5.2 (article DOI 10.3390/pharmaceutics14061268), mais la matrice `volunteer_id × LCMS-intensity` pré-dose n'a pas été vérifiée accessible ; [audit source](METFORMIN17_CLINICAL_OMICS_PK_OPEN_SOURCE_AUDIT_20261017.md).
3. **Degraeve tacrolimus n93** : PK (C0/AUC estimée) et 16S même cohorte, mais le supplément original fourni n'inclut pas le tableau intégral clinique↔identifiant SRA ; [audit](DEGRAEVE_WILEY_UPLOADED_SUPPLEMENT_AUDIT_20261016.md).

**Conclusion scientifique contrôlée** : Jarmusch établit l'existence d'une étude multi-modale avec vraie quantification ciblée PK, **pas** un fichier public analysable Omics→PK. Le blocage est la **traçabilité numérique des lignes individuelles**, pas une absence de métabolomique ou de PK mesurée dans le protocole.
