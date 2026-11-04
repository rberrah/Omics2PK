# Omics2PK — 80/90% de gain prédictif, budget de signaux et modules mécanistiques
**2026-11-04 · original-source experiments executed · exploratory, research-only.**

## Objectif corrigé

Le problème scientifique n'est pas `argmax importance SHAP des 384 signaux`. C'est :
> Conserver **80–90% du gain prédictif validé** par rapport à un comparateur clinique, avec **le plus petit nombre de mesures moléculaires**, puis demander si les signaux effectivement utiles correspondent à des voies physiopharmacologiques plausibles, sans présumer qu'un pic GC-TOF est un substrat d'un transporteur.

Deux métriques **NON équivalentes** :

- `SHAP_coverage(k)=Σ(top-k mean|SHAP|)/Σ(all molecular mean|SHAP|)` : combien de contribution descriptive à la sortie **du modèle complet** est attribuée aux k premières variables (dépend corrélation, modèle, taille de couche). **Ne signifie PAS 80% de l'amélioration.**
- `G(k)=MAE_clinical−MAE_clinical+omics(k)` : **gain prédictif évalué sur données indépendantes de l'entraînement**. Sélection 80%/90% seulement via `G_inner(k)/max(0,G_inner_all_candidate)`; si aucun gain positif en validation interne, **k=0**.

Le nombre de couches omiques retenues entre RNA/protéines/métabolites dépend d'un **coût par modalité**, par défaut égal à 1 en l'absence de tarif confirmé. On privilégie le moins d'analyses effectives, pas les molécules les moins coûteuses inventées. Une mesure multi-omique ne peut apparaître que si les différentes couches sont réellement appariées sur la même unité biologique.

## Évaluation rigoureuse des deux seuils

**Interne** : sur chaque pli d'entraînement externe uniquement, entraînement XGBoost `depth2, n_estimators100` (Rotroff), imputation médiane, calcul de **Tree SHAP d'entraînement**, classement des vrais noms de signaux source, réapprentissage de k=0,1,2,3,5,8,12,20,40,80,160,384 sur 4 plis internes, estimation MAE en unités originales. Choisir le plus petit k atteignant 80 ou 90 % du **meilleur gain intérieur positif** par rapport à `k=0`. Si pas de gain, `k=0`. Un patient externe n'influence **jamais** le choix de k.

**Externe** : 32 patients Rotroff, validation Leave-One-Patient-Out (32 plis). La courbe MAE par k fixe est descriptive pour l'analyse de sensibilité; choisir le meilleur k après avoir vu cette courbe **n'est PAS une preuve validée**. Les deux stratégies adaptatives ont leur prédiction individuelle tenue hors choix, avec CI descriptif par rééchantillonnage des **32 erreurs externes déjà calculées** (pas recalcul complet / pas test confirmatoire).

### Données humaines réelles : Rotroff 2016, metformine GC-TOF prédose

Source [Rotroff DOI](https://doi.org/10.3389/fphar.2016.00135), deux classeurs éditeur déjà verrouillés par SHA, même patient n32 avec AUC 0–24 et CLr réels; 384 signaux plasma, 155 noms chimiques source, 229 inconnus numérotés. Comparateur : âge, sexe, poids, créatininémie d'origine. Pas génotype, pas de PopPK, pas d'identification chimique des inconnus. **Vraie CI source exécutée et verte** : [#41918494437](https://github.com/rberrah/Omics2PK/actions/runs/41918494437).

| XGB / n32 patients | MAE AUC 0–24 (unités numériques source) | MAE CLr (unités numériques source) |
|---|---:|---:|
| **Référence clinique k0** | **2686.783611** | **114.473311** |
| **Un seul signal SHAP classé sur entraînement**, k fixé | **2578.496343** | **109.067142** |
| k2 | 2697.736795 | 115.439249 |
| k3 | 2797.623265 | 123.736787 |
| k5 | 2820.929741 | 119.576956 |
| k8 | 2709.708054 | 118.491525 |
| k40 | 2673.904866 | 115.526259 |
| k384 complet | 2692.811405 | 115.683381 |
| **80% gain interne, sélection imbriquée** | **2851.512011** | **118.934817** |
| **90% gain interne, sélection imbriquée** | **2839.891458** | **119.468654** |

**Décision scientifique :** k1 est une **observation a posteriori** de courbe de performance, pas une taille de panel retenue indépendamment. Comparé à k0, ΔMAE k1 AUC **−108.287 [−441.511,+203.683]** et CLr **−5.406 [−19.067,+7.907]**; les deux CI descriptifs recouvrent 0. Les stratégies imbriquées sont effectivement **moins bonnes que référence clinique** : AUC 80% Δ+164.728 [−38.317,+357.573]; CLr 80% Δ+4.462 [−3.940,+13.284]. Leurs taux P20 ne justifient pas d'usage clinique.

Le sélecteur imbriqué retient `k=0` pour 15/32 plis AUC et 14/32 plis CLr; moyenne k à 80% ≈2.56 AUC, 1.66 CLr; à 90% ≈15.06 AUC et 1.72 CLr. **La variabilité de choix est forte**. En particulier il ne faut pas dire « un métabolite explique 90% de la clairance ». Le premier signal à k1 du modèle est classé sur l'entraînement, pas une molécule définitivement identifiée; les bins source dominant le SHAP d'entraînement sont `222169` (AUC), `221574` (AUC), `218765` (CLr). Ce sont des **identifiants**, non des masses m/z.

**SHAP descriptif :** le **premier signal en moyenne ne représente que 28.6% de l'importance |SHAP| moléculaire** en AUC, et 33.8% en CLr, d'après les modèles complets des 32 plis d'entraînement. Il peut néanmoins améliorer numériquement une erreur hors échantillon, ce qui est la preuve que SHAP% ≠ gain%.

**Figure produite dans l'artifact CI :** `pk_nested_information_pareto.png` (MAE vs nombre de signaux, modèles 80/90, points Pareto descriptifs), `pk_shap_endpoint_feature_heatmap.png` (Top20 signaux, SHAP moyen normalisé par endpoint; PAS heatmap patient identifiable). Artifacts [#41918494437](https://github.com/rberrah/Omics2PK/actions/runs/41918494437).

### Nombre d'omiques : 47 lignées NCI-60, doxorubicine GI50 vraie **PD in vitro, pas PK**

RNA n17,987 gènes, protéines n3,171, métabolites n280, sources originales appariées par **nom exact de lignée** mais de **laboratoires/périodes différents**. 3 répétitions × 5 plis externes de lignées, chacun 4 plis internes, 15 choix de modèles, **n47 biologiques indépendantes**, et non n141 indépendants. XGB fixe et sélection train-only des variables par variance : RNA300, SWATH protéines300, métabolites120 par pli.

**Test source réel exécuté :** [#41920041480](https://github.com/rberrah/Omics2PK/actions/runs/41920041480).

| Ajustement + XGBoost | n couches omiques | MAE GI50 source |
|---|---:|---:|
| Type tumoral seul | 0 | 0.349593 |
| Métabolomique | 1 | 0.327104 |
| Transcriptomique | 1 | 0.333638 |
| Protéomique | 1 | 0.338472 |
| ARN + protéine | 2 | 0.328180 |
| ARN + métabolites | 2 | 0.328804 |
| Protéines + métabolites | 2 | 0.330808 |
| **Trois couches fixes** | 3 | **0.323019** |
| **80% gain interne: modèle modalité(s) choisi dans l'entraînement** | **1.0 en moyenne** | **0.340917** |
| **90% gain interne: modèle modalité(s) choisi dans l'entraînement** | **1.0 en moyenne** | **0.340917** |

La règle 80/90 sélectionne **les mêmes modalités** pour les 15 plis (metabolomics=6, RNA=4, proteomics=1, RNA+protéines=2, aucune=2). AUC/résultat PK ne peut être déduit du GI50. Le gain sélectionné 80/90 vs type tumoral seul n'est pas confirmé : ΔMAE −0.008676 IC descriptif **[−0.037633,+0.021688]**. Les trois couches fixes sont numériquement meilleures mais le choix par leurs résultats externes serait posthoc, non une justification de coûts réels.

**La conclusion de méthode est négative mais utile :** attribuer 90% d'un gain *en entraînement* ne garantit pas 90% du gain *sur nouvelles unités biologiques*, surtout si les mesures et laboratoires varient.

## Architecture « NetPK » : graphes physiopharmacologiques sous preuves

Conserver quatre graphes reliés mais **non confondus** :

```text
[Signal original RNA / protéine / métabolite / protéforme / PGx]
     ├─ Mesure et identité vérifiée (tissu, matrice, temps)
     ├─ Chemical ID (ChEBI, HMDB, InChIKey) si spectralement qualifié
     └─ Gène/protéine (UniProt) si expérimentalement qualifié
             ↓ Réactions Rhea/Reactome ; interactions STRING v12.5
         [Modules enzymatiques / transporteurs / PTM / voies]
             ↓ Hypothèse fonctionnelle, tissu et compartiment explicites
     [F, CL, V, absorption, fu, partition sang/plasma]
             ↓ Modèle PopPK / NetPK source
        [AUC, CLr, C0, Tmax, Cmax observés]
```

**STRING** [API officielle](https://string-db.org/help/api/) fournit des associations protéines-protéines, physiques/indirectes/régulatrices. Ne pas qualifier une arête *STRING functional* de **transport pharmacocinétique**, d'activité catalytique ou de causalité. **Rhea** décrit réactions chimiques vérifiées, **Reactome** les voies/processus, UniProt les protéines/PTM, et littérature pharmacologique expérimentale le contexte organe-substrat. Ontologie à chaque arête : `SOURCE_ASSAY` / `SPECTRAL_CANDIDATE` / `CURATED_REACTION` / `STRING_ASSOCIATION` / `LITERATURE_DRUG_MECHANISM` / `MODEL_SHAP_ASSOCIATION` / `FUNCTION_ASSAY` / `OUTER_PK_VALIDATION`. Toutes doivent conserver DOI/URL et grade.

Le prototype `src/omics2pk/physiopharmaco_cascade.py` représente les **axes OCT1/2-MATE1/2-K/PMAT pour metformine**, **CYP3A4/5-ABCB1 pour tacrolimus**, et les couches observées séparément. Un pic `222169` sans structure reste un nœud isolé du transporteur même si sa contribution SHAP est forte. Les réseaux de protéines peuvent être informatifs mais **ne se substituent pas** au graphe chimie→enzymes et PK.

L'expression d'un gène ne présume ni protéine pleine longueur/active ni PTM, ni transport fonctionnel. Protéomique/PTM d'albumine doit être traitée séparément pour liaison du médicament et fraction libre; une mesure de HSA totale n'est pas un proxy de fu. Les modules sélectionnés **a priori** par pharmacologie seront comparés à des **graphes témoins permutés** (même taille/nombre d'arêtes) pour vérifier si les priors physiques apportent réellement du gain.

## Association réseau vérifiée sur les vrais signaux de haut rang SHAP NCI-60

[STRING v12.5 original source execution SUCCESS #41920541070](https://github.com/rberrah/Omics2PK/actions/runs/41920541070), humain (taxon 9606), `required_score=400`, `network_type=functional`, `add_nodes=0`, 10 gènes RNA/protéines rangés SHAP à partir du vrai résultat n47 + 7 gènes *a priori* du mécanisme doxorubicine.

**10 associations fonctionnelles source vérifiées** : 7 entre gènes connus de pharmacologie, **2 reliant gènes SHAP et gènes pharmacologiques**, et 1 entre deux gènes SHAP. Les deux « ponts » retenus sont `LGALS3 ↔ ABCB1` (score 0.794) et `ATP8B1 ↔ ABCG2` (score 0.490). Association RNA/protéine SHAP `IL18 ↔ LGALS3` (0.452). Autres liens attendus: `TOP2A ↔ TOP2B` (0.999), `CBR1 ↔ CBR3` (0.999), `ABCC1 ↔ ABCG2` (0.944).

Il s'agit **d'associations STRING fonctionnelles**, lesquelles peuvent être indirectes, dépendantes de co-expression, de données textuelles ou d'annotations de bases externes. Elles **ne prouvent pas** liaison protéique directe, transport doxorubicine, inhibition des transporteurs, sensibilité causale, ni gain PK. Ce réseau n'est pas dérivé de mesures contemporaines chez un même donneur. Le fait d'avoir comparé les meilleurs signaux SHAP aux gènes pharmacologiques après résultats est exploratoire et expose à un biais de sélection.

`scripts/audit_nci60_string_shap_mechanistic_network.py` source versionnée, `src/omics2pk/physiopharmaco_cascade.py` graphe typé NetPK, `scripts/plot_nci60_string_shap_cascade.py` figure statique 17 gènes par catégorie (arêtes uniquement officielles, gènes isolés conservés). **Aucun « petit métabolite inconnu » n'est attaché à STRING** avant annotation structurelle et réactionnelle. L'enrichissement de voies futur doit utiliser un arrière-plan des **gènes réellement mesurés**, et une qualification indépendante de l'absence d'enrichissement artificiel par liste pré-sélectionnée.

## Statut de validation et prochaines preuves

1. Produit aujourd'hui : réelle sélection source n32 à 80/90 % au niveau **feature**, réelle sélection source n47 à 80/90 % au niveau **omics**, figures de Pareto et heatmap, source graph evidence contract et interrogation STRING v12.5 source exécutée (10 associations fonctionnelles, zéro preuve de causalité).
2. Non réalisé/non démontré : gain clinique PK externe, molécules inconnues Rotroff identifiées, activité des protéines inférée de SHAP, coût réel de chaque plateforme, protéomique & transcriptomique et métabolomique contemporaines mesurées chez même sujet PK, génotype complet et PopPK comparateur.
3. Prochaine étape pour la recherche : **prioriser les signaux stables et peu redondants, non le seul maximum SHAP** ; établir leur identité par spectre/RI; relier les identités réellement vérifiées aux modules OCT2/MATE/OAT ou CYP3A/P-gp, puis comparer la précision avec et sans prior physio sur une cohorte humaine distincte.

## Audit additionnel : stabilité des panels sélectionnés

Une courbe de Pareto et le gain à 80–90 % ne suffisent pas à montrer que les
*identités* moléculaires sélectionnées sont reproductibles. L'audit ajouté au
code du benchmark Rotroff (sans modifier les prédictions ou la sélection)
calcule maintenant :

- la similarité **Jaccard des top-10 SHAP entraînés dans les quatre plis internes**
  de chaque pli externe, puis sa moyenne descriptive ;
- pour chaque politique 80 % et 90 %, la fréquence des **panels vides**
  (`k=0`), la similarité Jaccard **conditionnelle aux deux panels non vides**,
  le nombre de signaux partagés et la fréquence de chaque identifiant source
  sélectionné, en agrégeant seulement les modèles entraînés hors sujet test ;
- l'intersection attendue de deux panels de même taille sous des ensembles
  uniformément aléatoires, comme **repère descriptif uniquement**.

**Interprétation obligatoire :** dans la LOOCV à n=32, deux entraînements
externes partagent 30 sujets sur 31. Une forte similarité de panels peut donc
provenir de ce recouvrement mécanique ; elle n'est **ni une réplication
indépendante, ni une preuve mécanistique, ni un test statistique contre un vrai
bruit de fond omique corrélé**. Les panneaux vides ne participent jamais au
calcul Jaccard (sinon deux absences seraient artificiellement notées comme
une signature stable). Les chiffres de stabilité du benchmark original ne
doivent être cités qu'après replay réussi de la CI, sans les inventer à
partir des cas synthétiques.

Les sorties sont **diagnostiques** : elles ne servent pas à choisir `k`,
ne remplacent pas la validation externe, ne démontrent pas un gain clinical
sur la métformine et ne changent pas la conclusion négative de la sélection
adaptative dans l'échantillon actuel.

## Code et QA
- `src/omics2pk/information_budget.py` et `modalities_information_budget.py` : gain 80/90% et solutions k=0 / omic=0.
- `src/omics2pk/rotroff33_nested_information_budget.py` et `scripts/validate_rotroff33_nested_information_budget.py` : vrai PK humain source.
- `scripts/benchmark_original_nci60_nested_modality_budget.py` : vraies omiques mesurées et GI50, validation 3x5 externe + 4 interne.
- `scripts/plot_original_information_budget.py` : figures agrégées sans patient.
- `src/omics2pk/physiopharmaco_cascade.py` et `scripts/audit_nci60_string_shap_mechanistic_network.py` : réseau de preuves et consultation STRING (pas une preuve d'activité).
- Tests synthétiques `tests/test_information_budget.py`, `tests/test_modalities_information_budget.py`, `tests/test_physiopharmaco_cascade.py`, `tests/test_nci60_string_shap_mechanistic_network.py` ; GitHub Actions source clinique et source NCI, sans fichier individu committé.
