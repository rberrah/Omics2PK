# Rosuvastatine n40 — état des vraies matrices de métabolomique MetaboLights MTBLS3128

**Audité le 17 octobre 2026 sur les fichiers source originaux**, en CI, sans publication des identifiants biologiques ni intensités individuelles.

- Étude Silveira et al. DOI https://doi.org/10.3389/fphar.2021.752960, MetaboLights **MTBLS3128**.
- [Source originale MAF/ISA vérifiée](https://github.com/rberrah/Omics2PK/actions/runs/39024110244) par `scripts/audit_rosuvastatin40_original_maf_intensity.py`, tests anti-fausse jointure et confidentialité en CI.
- Le tableau ISA des échantillons contient **142 Sample Names** ; les deux assays s'y joignent exactement pour **69** (ESI+) et **73** (ESI−) colonnes. Ces nombres concernent des **échantillons / fichiers LC-MS**, pas des personnes biologiques indépendantes (n40 dans article).
- Les deux vraies tables d'annotation **m_*.txt/MAF** ont chacune **15 lignes de métabolites annotés**, **90 colonnes** (69 identifiants échantillons) et **94 colonnes** (73 identifiants). Les **69 et 73 colonnes de samples sont parfaitement reconnues**, MAIS contiennent **zéro cellule numérique et zéro valeur non vide** dans leurs 15 lignes : ces MAF ne sont **PAS** une matrice traitée de concentrations ou d'intensités.
- **Ne pas attribuer artificiellement zéro aux signaux manquants.** Les cellules vides d'annotation ne prouvent ni absence biologique du métabolite ni LLOQ atteint ; ce n'est pas une donnée d'expression.
- Les **fichiers originaux `mzML` prédose** déposés par les auteurs via MetaboLights/OmicsDI restent l'option scientifique pour reconstituer la matrice : extraction de pics, correction bruit et dérive, regroupement m/z–temps de rétention, imputation/limites de détection et contrôle QC, deux ionisations séparées. Documentation d'un pipeline local MZmine / XCMS / MS-DIAL ou adaptation MultiOmics, avec gel des paramètres avant validation des résultats PK.
- **Second verrou** distinct : aucune table d'AUC/Cmax individuelle, dose/concentrations rosuvastatine avec **une vraie clé volontaire** n'a été authentifiée dans les fichiers publics vérifiés. La jointure ne doit pas être créée à partir de l'ordre des fichiers ou graphiques.
- Éviter les fuites : toute sélection de métabolites, correction de batch dépendant des individus, projection latente et hyperparamètres doit être entraînée **dans les seuls individus de chaque pli**, avant heldout scoring. L'article évalue des milliers de configurations, donc ses scores annoncés ne sont pas une preuve d'amélioration validée indépendamment.

**Verdict final de ce contrôle :** `PRE_DOSE_ORIGINAL_RAW_OMICS_AVAILABLE; ISA_EXACT_SAMPLES_MATCH; MAF_HAS_NO_PROCESSED_INTENSITY; ORIGINAL_INDIVIDUAL_PK_TARGET_NOT_JOINED; E3_NUMERIC_BENCHMARK_ABSTAIN`.

Rapport de la source et du protocole : [ROSUVASTATIN40_MTBLS3128_REAL_PREDOSE_OMICS_PK_SOURCE_AUDIT_20261017.md](ROSUVASTATIN40_MTBLS3128_REAL_PREDOSE_OMICS_PK_SOURCE_AUDIT_20261017.md).
