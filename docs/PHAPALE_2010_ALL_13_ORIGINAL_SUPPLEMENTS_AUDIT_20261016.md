# Audit réel des 13 suppléments Phapale (2010), pharmacométabolomique urinaire → PK tacrolimus

**Source reçue de l'utilisateur le 16 octobre 2026** : `paphale.zip`, 1 792 400 octets, SHA-256 **`b180d49fcde75d24b6d1812be001dbfa648c9565d0f5bb7100e1bc7d5c299a77`**, contenant exactement les **13 fichiers originaux** `cptclpt2009296-sup-0001.doc` à `...-0013.doc`, originaux Word binaires OLE2 (2 688 000 octets non compressés, dossier compris).

Publication primaire : Phapale PB et al., *Clinical Pharmacology & Therapeutics* 2010;87:426–436, DOI https://doi.org/10.1038/clpt.2009.296. **29 volontaires** avec 1 256 ions LC-MS urinaires **pré-dose** et concentration tacrolimus **plasmatique jusqu'à 72h** après administration orale. C'est un protocole E3 pré-dose→PK réellement compatible avec une prédiction *prospective* ; l'étude ne constitue pas E4 omiques dynamiques, puisque les omiques ne sont pas documentées comme plusieurs visites temporelles.

## Inventaire du contenu original Word (contrôlé localement)

| Fichiers reçus | Type réellement examiné | Information |
|---|---|---|
| 0001 | Figures supplémentaires 1a/1b | AUC PLS prédite vs observée par validation interne 7 groupes ; indices VIP de 42 métabolites retenus |
| 0002–0003 | Figures supplémentaires 2–3 | Clearance ; comparaison, projections PLS, corrélation avec cortisol / 24-hydroxycholestérol, validation interne et permutation |
| 0004–0007 | Figures supplémentaires 4–7 | Chromatogrammes, distribution m/z/retention time des 1 256 ions et vérifications MS/MS d'identifications |
| 0008 | Figure supplémentaire 8 | Régression PLS à quatre métabolites et **20 permutations de la réponse AUC** |
| 0009 | Tableau 1, **image intégrée** | Résumés des composantes du premier modèle PLS AUC ; R²X/R²Y/Q² |
| 0010 | Tableau 2, **image intégrée** | Résumés des composantes du second modèle PLS AUC |
| 0011 | Tableaux 3a/3b, **images intégrées** | Critères PLS pour la clearance et le Ctrough |
| 0012 | Tableau 4, **image intégrée** | Metabolites putativement ou effectivement identifiés, rapports m/z et temps de rétention ; **pas les 29 intensités individuelles** |
| 0013 | Méthodes supplémentaires, texte | UPLC-MS/MS tacrolimus, LC-MS pré-dose, peak detection/alignment XCMS, normalisation/valeurs manquantes, PLS, sélection VIP et validation |

Les treize anciens `.doc` ont été lus avec `antiword` et convertis **uniquement dans l'environnement local privé** en documents `.docx` pour inspecter la structure interne des figures/tableaux. **Aucune table Word structurée `w:tbl` ni `word/embeddings` Excel exploitable dans les treize conversions** ; tableaux 0009–0012 insérés comme **figures raster**. Même les figures `AUC calculée versus AUC prédite` ne possèdent pas de correspondance `volontaire_id ↔ urine_id ↔ AUC/CL/C(t)` démontrée par une table source. Les fichiers .doc d'origine ne doivent pas être traités comme dépôts de données numériques patients parce qu'une figure comporte des points.

**Verdict de données :** les 13 fichiers **n'apportent pas** (i) matrice individuelle `29 × 1 256` d'intensités LC-MS, (ii) tableau dosages TAC `volontaire_id × temps × concentration` pendant 72 h, (iii) valeurs NCA individuelles AUC/Cmax/Tmax/CL avec source/censure, ni (iv) identifiant individuel permettant de joindre X aux PK mesurées. **Aucun modèle Omics→Cmax/AUC/CL ne doit donc être recalibré à partir de ces suppléments seuls**. Ce constat porte sur le contenu accessible et structuré de ces 13 documents, pas sur la disponibilité éventuelle de données détenues par les auteurs.

## Protocole original et comparaison avec les exigences de validation moderne

Les méthodes authentiques du fichier **0013** établissent :
- **XCMS 1.16.3 / R 2.9.0** ; signaux LC-MS positifs m/z 50–800 ; détection et alignement, signal/bruit 10 ; intensités `log2` par ion normalisées **par quantiles sur le jeu considéré**, nettoyage des pics instables, valeurs manquantes ;
- PLS en deux étapes : d'abord les **1 256 caractéristiques**, puis sélection selon le **VIP** (seuil décrit **VIP >1,5**) et second modèle ; normalisation des caractéristiques et transformation logarithmique des cibles PK ;
- validation décrite comme **sept parties** (6/7 pour entraîner et 1/7 pour prédire) et **20 permutations** ; les légendes de figures évoquent parfois leave-one-out et sept groupes. Les méthodes doivent être considérées comme des **validations internes**, jamais un site ou une cohorte externe ;
- Les marqueurs du dernier modèle quatre variables comprennent **cortisol, 1-méthylguanosine, acétyl-arginine, phosphoéthanolamine**. Ces noms sont des **hypothèses biologiques publiées**, pas un tableau d'abondance pour les 29 individus ; cortisol/24-hydroxycholestérol apparaissent également dans les graphes relatifs à la clearance.

**Point méthodologique pour Omics2PK :** les méthodes décrivent quantile normalization, sélection VIP et recomposition du PLS successivement, puis validation interne. Le texte consulté **ne démontre pas que le filtrage des 1 256 ions, les normalisations et la sélection VIP ont tous été intégralement recalculés dans chaque pli externe**. Ces étapes pourraient créer une fuite de données si elles utilisent les sujets de test. **Ce n'est pas une accusation que l'étude a nécessairement mal fait**, mais un risque méthodologique qu'Omics2PK doit éliminer explicitement lors d'une reproduction moderne. De plus, quatre molécules sélectionnées parmi 1 256 caractéristiques avec seulement 29 participants exigent un protocole strict et des intervalles de variation/une validation indépendante.

### Contrat de reproduction si matrices réelles reçues

1. `volunteer_id`, dose orale réelle (mg/kg et mg, poids), temps d'administration, formulation, `sampling_time_h`, concentration plasma en unités originales, censure/LOQ, calendrier complet 0–72h, réexposition et mesures exclues.
2. Tableau `volunteer_id × ion_id(m/z,RT) × intensity` **avant** imputation/quantile normalization, acquisition pré-dose certifiée ; batch, dilution urinaire/créatinine, qualité MS, adducts et ion identities ; version du pipeline XCMS.
3. Endpoint individuel source `AUC0–72`, `Cmax`, `Tmax`, `CL/F` (seulement si estimable), avec méthode de calcul, fenêtre, unités et doses. Ne pas reconstruire une clairance IV à partir d'une dose orale en supposant F=1.
4. **Validation E3 tenue hors individu entier** : clinique (poids/dose, éventuellement PGx et covariables connues) → 4 ions **préspécifiés** → sélection sur train seulement parmi les 1 256 ions, avec pipeline preprocessing entraîné sur train et projection figée sur test ; erreurs de concentrations C(t) et d'AUC NCA selon protocole, bootstrap au **niveau patient**, contrôles de permutation *emboîtés* si sélection multivariée, et évaluation d'une cohorte externe si possible.
5. Aucun lien artificiel avec les patients **Degraeve 93**, **Wang 36**, **He 109**, ni avec **mavoglurant 120** : sujets, tissus et médicaments comparables ne suffisent pas à créer un même individu.

## État Github et sécurité

Le dépôt public Omics2PK ajoute uniquement :
- `scripts/audit_local_phapale_supplements.py` : empreintes **13 fichiers et archive**, forme OLE2, texte scientifique via `antiword` **localement**, sortie de seules métadonnées agrégées et hypothèses publiées ;
- `tests/test_local_phapale_supplements.py` : tests **ZIP synthétiques uniquement** (pas de PDF, DOC, métabolites, graphiques originaux ou volontaires publiés) ;
- CI `qualified-local-phapale-supplements.yml` : vérifie absence de documents Wiley dans Git et contrats de sécurité.

Le supplément Degraeve également transmis, `cpt3077-sup-0001-supinfo(1).zip`, possède **exactement** SHA256 `03aab5deb292a03753073e493d7ae49896884be2cde22d5c2d1c9f3037488cc2` et est donc identique à l'archive Degraeve déjà qualifiée. Aucun nouvel échantillon ni fichier clinique n'y est ajouté.

**Statut final : `ORIGINAL_13_SUPPLEMENTS_AUDITED / PROSPECTIVE_STUDY_RELEVANT / INDIVIDUAL_MATCHED_PK+OMICS_DATA_ABSENT / E3_MODEL_ABSTAIN`.**

Voir également [priorité jeux publics](PUBLIC_PREDOSE_TAC_NCA_AND_MAVO_LONGITUDINAL_RESULTS_20261016.md) et [ticket #31](https://github.com/rberrah/Omics2PK/issues/31).
