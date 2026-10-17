# Huit CYP mesurés : validation original EPA n100 avec correction des comparaisons multiples

Date : **17 octobre 2026**. **Status : exécution réelle originale vérifiée en CI, résultats E2 précliniques**.

Source principale : Liu et al., *Xenobiotica* 2021, DOI https://doi.org/10.1080/00498254.2020.1867929 ; dépôt authentique U.S. EPA DOI https://doi.org/10.23719/1519197 , classeur SHA-256 **cfc9dc334064ea61e99ca6fd39c2c8ef7c621d000a37c4da2c95f26ca85e36dd**. Réplicas indépendants de la régression Ridge, 100 véritables donneurs avec deux identifiants identiques dans les feuilles d'activité et qPCR, 94 à 100 complets selon le CYP. Le logiciel ne publie aucune ligne individuelle.

**CI original n100 vert** : https://github.com/rberrah/Omics2PK/actions/runs/39016653259. Pipeline : src/omics2pk/epa_eight_cyp_sensitivity.py , scripts/validate_epa_eight_cyp_multiplicity.py et tests/test_epa_eight_cyp_sensitivity.py . Les huit erreurs LODO reconstruisent les MAE antérieurement vérifiées à **six décimales** avant les nouvelles analyses.

## Hypothèse et méthode

Chaque valeur ci-dessous mesure la différence de MAE en **pmol/min/mg de protéines** entre modèle référence âge/sexe et modèle âge/sexe + transcription qPCR du CYP correspondant. Une différence négative est une amélioration de la prédiction du taux d'activité enzymatique in-vitro.

- Régression Ridge alpha=1 fixé ; log1p des valeurs d'ARN et de l'activité, variables standardisées uniquement chez les donneurs train, 100 validations laisse-un-donneur-de-côté (n94–100 selon disponibilité).
- Comparaisons **appairées au même donneur** pour chaque CYP. Les **huit endpoints sont liés chez les mêmes 100 donneurs**, et ne représentent donc pas huit cohortes indépendantes.
- Échantillonnage bootstrap de **5 000 tirages des mêmes 100 identifiants biologiques simultanément pour les huit CYP**, seed 20261017. Pour chaque enzyme, on ne retient que les donneurs réellement mesurés, avec leurs pertes déjà obtenues en prédiction LODO.
- Intervalle percentile 95% exploratoire, et 99,375% avec correction de **Bonferroni pour huit critères** (alpha bilatéral 0,05/8). La protection famille de tests suppose une couverture acceptable des intervalles bootstrap individuels.
- Aucune prétention d'inférence causale, d'équivalence thérapeutique ou de performance de doses chez les patients.

## Résultats originaux confirmés (donneur exclu)

| CYP mesuré | n | ΔMAE ARN − âge/sexe | IC95% exploratoire | IC99,375% Bonferroni | IC corrigé exclut 0 |
|---|---:|---:|---:|---:|---|
| CYP1A2 | 98 | **−5,589** | [−8,966 ; −2,284] | **[−10,310 ; −1,087]** | Oui |
| CYP2B6 | 100 | **−11,757** | [−17,922 ; −5,760] | **[−20,195 ; −3,224]** | Oui |
| CYP2C9 | 100 | −19,673 | [−54,219 ; +15,744] | [−66,004 ; +30,815] | Non |
| CYP2C19 | 94 | −1,530 | [−3,854 ; +0,788] | [−4,724 ; +1,706] | Non |
| CYP2D6 | 100 | −18,003 | [−32,212 ; −2,686] | **[−38,000 ; +2,821]** | Non |
| CYP2E1 | 100 | **+4,677 (dégradation)** | [−25,098 ; +37,394] | [−35,434 ; +51,672] | Non |
| CYP3A4 (activité possiblement 3A4/5) | 99 | **−337,278** | [−554,032 ; −121,585] | **[−646,051 ; −37,048]** | Oui |
| CYP4A11 (activité pouvant inclure 4A9) | 98 | −16,170 | [−41,938 ; +9,175] | [−51,636 ; +18,529] | Non |

**Résumé empirique :** 7/8 améliorations ponctuelles ; 4/8 intervalles 95% non corrigés exclus zéro ; seulement **3/8** intervalles conservateurs corrigés excluent zéro (CYP1A2, CYP2B6, CYP3A4). Le signal CYP2D6 n'est pas suffisamment robuste à la correction multiple. Pas de moyenne globale des MAE, car les plages et unités d'activité varient selon CYP/sonde ; huit CYP d'un même foie ne sont pas indépendants.

## Limites strictes

1. Il s'agit d'un **bootstrap apparié des huit séries de pertes prédictives déjà estimées**. Le modèle n'est **pas réentraîné** dans chaque réplication ; aucun intervalle présenté ne prétend valider la variabilité des paramètres NLME ou l'incertitude externe sur de nouveaux centres/donneurs.
2. Les intervalles percentiles Bonferroni n'offrent pas de garantie exacte finie de contrôle d'erreur sous distributions très asymétriques ; les extrêmes par donneur restent importants. Aucune correction causale du biais de sélection, lots, comédications ni inflammation n'est réalisable sans covariables originales correspondantes.
3. Les CYP du même donneur partagent des facteurs biologiques et des sources d'erreurs ; bootstrap synchronisé par donneur, mais ce n'est pas une validation *leave-site-out* ni une validation externe. Génotypes et exposome patient ne sont pas disponibles dans le document.
4. ARN = **qPCR ciblée**, non ARNseq, SWEET ou multi-omique. Résultat = **activité microsomale enzymatique mesurée**, et **pas PK plasma**, dose/C(t), AUC, CL/F, MIPD clinique, jumeau digital ou causalité.
5. Le cas CYP2E1 dégradé et les pertes de correction pour CYP2D6/2C9/2C19/4A11 font partie intégrante de l'évidence et ne doivent pas être masqués.

## Suites

- Renforcer l'estimation par **bootstrap refitté** si le temps de calcul le permet, ou comparaison split-grouped holdout sur un second véritable jeu de donneurs.
- Investiguer le rôle CYP2E1 et la calibration des modèles, mais ne jamais optimiser les covariables à partir du même jeu de test déclaré externe.
- Tester le mécanisme NAT2 du sulfaméthoxazole séparément sur son tableau original de n26, sans inventer de lignes et en respectant les droits de l'article (suivi issue #38).
- **E3 et E4 cliniques bloqués** jusqu'à obtention d'un appariement authentique entre omiques individuelles de pré-dose, dose/temps/C(t) et contexte clinique/PGx, qui n'est pas présent dans ce classeur EPA.

