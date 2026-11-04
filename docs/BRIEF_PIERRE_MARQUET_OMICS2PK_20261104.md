# Omics2PK — résultats humains et proposition d'étude à discuter avec Pierre Marquet
**Note de réunion scientifique — 4 novembre 2026**  
**Statut :** résultats exploratoires *in vitro* vérifiés sur des données humaines ; **aucune supériorité clinique de la multi-omique sur le génotype n'est établie**.

## New experimentally executed result — NCI-60 RNA/proteome/metabolome vs measured doxorubicin GI50 (2026-11-04)

**Original-source replayed success:** [CI #41847348137](https://github.com/rberrah/Omics2PK/actions/runs/41847348137), [scientific results report](NCI60_THREEOMICS_ORIGINAL_GI50_RESULTS_20261104.md). Same **47 source cell-line identities** (not necessarily same experimental aliquot), original **17,987 ARN / 3,171 SWATH proteins / 280 metabolites**. Fixed identical 3×5 train-only PCA/Ridge folds.

| GI50 MAE (lower is better) | Model |
|---:|---|
| **0.308916** | ARN seul (sans classe tumorale) |
| 0.364261 | Protéome seul |
| 0.382598 | Métabolome seul |
| **0.334666** | Classe tumorale + ARN |
| 0.341053 | Classe + ARN + protéines |
| 0.351537 | Classe + ARN + métabolites |
| **0.358124** | Classe + trois omiques |

La fusion à trois omiques ajoute **+0.023458 de MAE** par rapport à classe + ARN, IC descriptif 95 % **[+0.009865 ; +0.037185]** pour l'écart apparié de pertes CV (1999 resamples de lignées). **Résultat négatif scientifiquement informatif** : cette association de sources originales, obtenues dans des études et conditions de culture différentes, ne montre pas de valeur ajoutée de la multi-omique dans ce protocole. **PD cellulaire**, pas PK clinique. Ce résultat exclut toute affirmation de supériorité universelle de l'ARN ou d'inutilité générale des omiques multiples.

**Message oral pour Pierre Marquet :** nous avons maintenant non seulement des preuves ARN→fonction CYP (HLC/EPA), mais aussi un contre-exemple source-vérifié où accumuler ARN+protéines+métabolites dégrade la prédiction pharmacologique. Le verrou devient la qualité, les isoformes, les protocoles et l'identité biologique des mesures, **avant** le choix de l'algorithme.

## Message en une phrase
L'expression hépatique ciblée contient une information prédictive sur l'activité enzymatique, **mais davantage de variables omiques ne signifie pas nécessairement une meilleure prédiction** ; la preuve pertinente pour les immunosuppresseurs reste de montrer un gain sur une mesure PK indépendante chez les mêmes patients.

## Résultats déjà calculés — à présenter sans surinterprétation

| Expérience | Données humaines | Comparaison appariée, CV entre donneurs | Résultat | Niveau de preuve |
|---|---|---|---|---|
| **HLC : ARN → activité** | 387–397 foies, 10 activités CYP mesurées | Âge/sexe vs âge/sexe + sondes ARN du CYP correspondant, Ridge alpha 1, validation à 5 plis par donneur | **10/10 activités** : baisse de MAE de **3,73 à 16,29 %**. CYP3A–midazolam **851,109 → 719,247** (−15,49 %). Panel élargi CYP/POR/transporteurs **moins bon que l'ARN ciblé sur 9/10** activités. | Données d'origine acquises légitimement et appariées ; original public Yang et al. 2010 DOI **10.1101/gr.103341.109**. **Recherche interne**, pas validation externe. |
| **HLC : génotype vs ARN** | **171** foies pour CYP3A–midazolam, **172** pour CYP2C9, tous avec un SNP spécifique reconnu | Mêmes donneurs et trois répétitions de cinq plis fixes : âge/sexe, +ARN, +SNP, +ARN+SNP. | CYP3A MDZ : **706,771** (âge/sexe), **697,803** (+rs776746), **641,842** (+ARN CYP3A4/5), **641,836** (+ARN+SNP). CYP2C9 : **182,367** vs **183,260** (+rs1799853) vs **172,270** (+ARN). | Génétique **très partielle** : dosage rs776746 et rs1799853 uniquement ; ni allèles *3/*6/*7 complets ni haplotypes PGx validés. Ne démontre **pas** que l'ARN bat la pharmacogénétique complète. |
| **EPA : validation sur autre banque** | 100 foies, 8 activités CYP avec 94–100 donneurs complets | Âge/sexe vs +qPCR du CYP ; Ridge LODO. 5 000 rééchantillonnages **des pertes CV** avec correction Bonferroni sur 8 endpoints. | **7/8** baisses ponctuelles de MAE ; **3/8** IC corrigés 99,375 % excluent zéro : CYP1A2, CYP2B6 et CYP3A4. CYP2E1 se dégrade. | **Cohorte d'origine distincte de HLC**, mais qPCR ciblée, pas protéome ni véritable PK. Source ouverte EPA DOI **10.23719/1519197**. Bootstrap des pertes, pas bootstrap réajusté. |
| **Latham : protéine + fonction → CL métabolique** | **21** HLM mêmes donneurs pour tucatinib, **17** imatinib quantifiables | Moyenne d'apprentissage vs protéines CYP2C8/CYP3A vs activités sondes CYP vs combinaison, Ridge LODO. | Tucatinib (µL/min/mg) **6,246 → 2,041 (protéines) → 1,783 (protéines+activités)**. Imatinib **28,870 → 14,166 → 13,620**. IC bootstrap descriptifs du gain de combinaison par rapport à protéines seules **incluent zéro** pour les deux CL. | Mesures originales **d'activité et CL intrinsèque *in vitro***. Protéines+activité **n'est pas multi-omique** ; les deux médicaments partagent la même banque de foies. Articles DOI **10.1016/j.dmd.2025.100061** et **10.1016/j.dmd.2025.100196**. |
| **Latham : génotype → imatinib PHH** | **13** hépatocytes primaires complets sur 14 | PGx CYP2C8+CYP3A5 vs +protéines vs +sondes d'activité | MAE **2,907** PGx ; **3,004** +protéines ; **2,763** +sondes enzymatiques ; **3,263** +les deux. | **Très faible n**, données originales et tableau 4 reproduits ; rôle du CYP2C8 ; non clinique ; effet non généralisable. |
| **GSE49012 : miARN → CYP3A** | 34 foies avec miARN et activité dans GEO | miARN préspécifiés et PCA en plis d'entraînement. **AUCUNE inférence biologique valide** | Les valeurs d'activité de 22 foies malades totalisent **12,533**, incompatibles avec la moyenne publiée de 0,28 sur n28 si même mesure et vrais sous-échantillons. | **ABSTAIN_SOURCE_PARITY_FAILURE**. Contrôle qualité utile, surtout ne pas montrer comme effet omique positif/négatif. |

**Attention aux unités :** ne pas comparer les MAE absolues entre CYP ou études. La baisse en pourcentage est intra-activité, sur le même sous-ensemble de donneurs. N'employer aucun p-value non corrigé comme démonstration clinique.

## Le résultat original le plus convaincant pour ouvrir la discussion

**HLC CYP3A–midazolam, même n=171 :**
- Âge/sexe **706,771** MAE ;
- +SNP rs776746 **697,803** (gain **1,27 %**) ;
- +ARN CYP3A4 et CYP3A5 **641,842** (gain **9,19 %**) ;
- +ARN+SNP **641,836** (gain additionnel négligeable).

**Ce que cela suggère :** une mesure moléculaire d'expression peut capturer une partie du phénotype enzymatique non représentée par *ce seul* polymorphisme. **Ce que cela ne montre pas :** le gain au-delà du génotypage clinique complet, ni un gain de prédiction de l'AUC du tacrolimus.

**Deuxième question :** pourquoi une sélection large de CYP/POR/transporteurs n'améliore-t-elle que rarement l'ARN du CYP cible ? Les hypothèses à départager sont colinéarité, dimension excessive, effectif insuffisant, régulation effectivement ciblée ou artefacts de prétraitement ; ne pas conclure à l'absence d'effets biologiques des autres gènes.

## Proposition d'étude « omiques → fonction → PK » en deux étapes indépendantes

**E2 / POC fonctionnel (immédiat, sans nouveaux patients) :** HLC et EPA, éventuellement source Zhang n100 / Liu n55 si matrices authentifiées. Comparer sur individus identiques : clinique seule → PGx réellement disponible → ARN CYP cible → plus protéines indépendantes → plus metabolome endogène *préincubation* → score biologique sourcé. Endpoint **activité CYP mesurée**, idéalement **sonde sélective CYP3A5**, sinon libellé CYP3A partagé. Contrôles négatifs : permutation inter-donneurs des mesures omiques, sondes spécifiques vs partagées, témoin fonctionnel d'un médicament où CYP3A5 n'est pas déterminant.

**E3 / validation translationnelle (cohorte immunosuppresseurs) :** chez les **mêmes individus et à des temps documentés**, relier (a) dose, temps, concentrations et AUC/CL estimée par un protocole défini, (b) CYP3A5 complet (*3/*6/*7) et CYP3A4 pertinent, (c) ARN/protéines disponibles en tissu ou dans un substitut biologiquement justifié, (d) métabolomique **prédose**, hématocrite, CRP/inflammation, fonction hépatique et rénale, interactions azolées, formulation et temps post-greffe. **Tacrolimus en priorité ; MPA/MPAG comme mécanisme distinct**, sans transposer directement les biomarqueurs CYP3A.

Baselines préenregistrées :
`clinique → clinique+PGx → +ARN → +protéome → +métabolome → +score biologique sourcé → modèle ML multivue`, **mêmes sujets** et mêmes partitions patient-groupées ; apprentissage du prétraitement exclusivement chez les sujets d'entraînement, évaluation externe si une seconde cohorte est possible. Performance prioritaire : amélioration appariée de MAE de CL/AUC dans les mêmes unités, calibration et incertitude, non pas nombre de hits omiques.

**Prérequis majeur :** l'ARN/protéine circulant ne mesure pas automatiquement l'activité hépatique CYP3A ; tracer la provenance anatomique et temporelle et conserver une valeur manquante lorsque non mesuré. Les métabolites d'une incubation *sur le médicament étudié* sont des **résultats**, pas des prédicteurs prédose indépendants.

## Ce qui est demandé à Pierre Marquet

1. **Valider la cible translationnelle** : immunosuppresseurs, tacrolimus comme test de CYP3A5 ; MMF/MPA/MPAG comme second système mécanistique.
2. **Confirmer la faisabilité du chaînage BIOMARGIN ↔ ABIS / TDM** : consentements, identifiant pseudonymisé, date de prélèvement omique, dose/horodatages, génotype, concentrations, fenêtre de stabilité et co-traitements. La disponibilité réelle de ces paires n'est **pas** confirmée à ce stade.
3. **Définir le minimum expérimental** avant un plan de puissance : CYP3A5 *3/*6/*7, métabolomique prédose, un test fonctionnel sur sous-cohorte si réalisable, AUC/CL mesurée ou estimée avec un protocole indépendant ; ne pas promettre du RNA-seq/protéome hépatique à partir de simples échantillons sanguins.
4. **Demander les sources supplémentaires** Zhang n100 et Liu n55 ; [manifeste d'acquisition](PUBLIC_COHORT_SOURCE_ACQUISITION_20261104.md). En parallèle, tenter les 61 HLS9 CYP2C9/FABP1 de PXD085173 et leur activité mesurée avec appariement officiel.

**Décision attendue de la réunion :** accès aux échantillons et données **effectivement appariables**, choix d'un endpoint E3 primaire, et possibilité d'une analyse de preuve de concept *à effectif et puissance préspécifiés après qualification des données*.

## Matériel scientifique, scripts et sources

- [HLC : résultats originaux ARN, génome et activité](HLC_ORIGINAL_400_LIVER_CYP_RNA_GENOTYPE_RESULTS_20261018.md), source Yang et al. 2010 DOI 10.1101/gr.103341.109.
- [EPA 100 foies / huit CYP / correction multiple](EPA_100_LIVER_EIGHT_CYP_BONFERRONI_REAL_RESULTS_20261017.md), dépôt EPA DOI 10.23719/1519197.
- [Latham : deux médicaments et 21 HLM](HEPATIC_TWO_DRUG_CYP_PROTEIN_ACTIVITY_FUSION_REAL_20261018.md), origine source dans GitHub Actions #39147377727.
- [Latham : 13 PHH avec génotype](LATHAM_PHH_SOURCE_REPLAY_REAL_20261104.md), original-source CI #41805204594.
- [GSE49012 — source-parity *ABSTAIN*](GSE49012_ORIGINAL_MIRNA_CYP3A_SOURCE_PARITY_20261104.md), issue #138.
- [Cohortes et demandes sources en attente](PUBLIC_COHORT_SOURCE_ACQUISITION_20261104.md), issue #63.

**Formulation scientifiquement acceptable à l'oral :** « Nos premières validations sur des foies humains montrent que l'expression du CYP apporte de l'information pour prédire la fonction enzymatique, mais qu'une intégration plus large n'apporte pas mécaniquement de gain. Nous cherchons maintenant à démontrer — ou réfuter — l'intérêt incrémental des omiques sur la PK des immunosuppresseurs, au-delà de la clinique et de la pharmacogénétique. »
