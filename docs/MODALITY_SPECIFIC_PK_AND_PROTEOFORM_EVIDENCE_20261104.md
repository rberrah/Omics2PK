# Omics2PK — modalité par modalité, fonction biologique et PK : état réel au 4 novembre 2026

## Question
**Une mesure moléculaire circulante prédose améliore-t-elle la prédiction de la PK individuelle, par rapport au modèle clinique + PGx + PopPK, et quelle valeur marginale apportent séparément RNA, protéines/protéformes, métabolites et leur intégration ?**

**Règle scientifique absolue :** comparer toutes les modalités disponibles sur **les mêmes participants, le même médicament, le même endpoint PK et les mêmes plis externes**, sans attribuer à une cohorte les couches omiques qui n'y ont pas été mesurées. Des conclusions d'activité enzymatique in-vitro, de GI50 PD cellulaire et de PK plasmatique clinique ne peuvent pas être mises sur une seule échelle de MAE.

## Pourquoi une protéome/protéforme est distincte d'un petit signal métabolique
L'article El Balkhi S / Berrah R et al. *Scientific Reports* 2026, doi:10.1038/s41598-026-57614-y, montre que les isoformes d'albumine et leurs PTM, caractérisées par **LC-HR-MS top-down**, apportent une information sur la fibrose hépatique versus la mesure total d'albumine et FIB-4 (**172 malades + 82 témoins; endpoint fibrose/diagnostic, pas clairance médicamenteuse**). Les PTM peuvent modifier l'affinité albumine–ligand; pour un médicament fortement lié aux protéines plasmatiques c'est un mécanisme potentiel d'effets sur **fraction libre fu, distribution V, clairance et exposition libre**, conditionnel à cinétique du médicament et au principe de foie bien brassé. Les protéines de liaison ne doivent pas être confondues avec les transporteurs membranaires.

Rotroff *Front Pharmacol* 2016, doi:10.3389/fphar.2016.00135, **GC-TOF-MS métabolomique** prédose : 384 *features* (155 noms chimiques, 229 chiffres inconnus). Ce sont des petits signaux chromatographiques/fragmentaires inconnus et pas, en tant que tels, **229 albumines modifiées** ni des peptides de P-gp. Des molécules dérivées d'une modification protéique ou d'une protéolyse pourraient théoriquement figurer parmi de petits métabolites, mais le signal brut **ne prouve ni cette origine, ni une structure, ni une réaction ou un transporteur**. Ne jamais attribuer un signal numéroté à un métabolite/transporteur par imputation nominale.

**Metformine** : liaison albumine négligeable et élimination rénale inchangée; cibles transport OCT1/SLC22A1 (distribution hépatique), OCT2/SLC22A2 et MATE1/2-K (sécrétion rénale), PMAT/SLC29A4 et d'autres transporteurs intestinaux. **ABCB1/P-gp n'est pas la cible transport prioritaire**, contrairement au **tacrolimus**, où ABCB1 intestinal joue sur la biodisponibilité en interaction avec CYP3A4/5. Prévoir des **contrôles négatifs de voie** et séparer génotype, expression tissulaire, abondance membranaire, **état PTM** et **activité réelle**.

**Tacrolimus** : forte partition dans les érythrocytes; ne pas assimiler PTM de l'albumine à une modification directement démontrée du tacrolimus sang total. Pour un test fu/protéforme, privilégier **concentrations plasmatiques libres, fu ou liaison ex vivo** du médicament choisi; vérifier la contribution relative AAG/albumine/hématocrite.

## Ce que les données originales démontrent ou ne démontrent PAS actuellement

| Couche | Cohorte/épreuve source exécutée | Endpoint réellement observé | Résultat exploitable | Preuve que la couche améliore PK plasmatique humaine ? |
|---|---|---|---|---|
| **Génétique (PGx)** | HLC n171 CYP3A; variant source rs776746 uniquement | Activité microsomale de CYP3A (midazolam), **pas AUC** | MAE age/sex 706.771; +SNP 697.803; +RNA CYP3A4/5 641.842; +RNA+SNP 641.836 | **Non**; un SNP n'est pas génotype/diplotype exhaustif |
| **Transcriptomique** | HLC n396 foies (CYP3A midazolam) | Activité CYP3A in vitro | MAE âge+sexe 851.109 → +ARN ciblé 719.247 (–15.5 %). 10/10 activités CYP ciblées s'améliorent numériquement. | **Non**; aucun sang circulant ni PK plasmatique individuelle |
| **Protéomique** | Wegler n15 donneurs, hépatocytes humains | CLint in vitro, cinq médicaments distincts | Midazolam log-MAE référence moyenne 0.4350 → +CYP3A protéine 0.1626; bénéfice absent ou faible pour d'autres substrats et en microsomes | **Non**; preuve d'utilité **in vitro** conditionnelle au tissu/drug |
| **Protéine vs PGx** | Latham 2025, PHH n13 complets | CLint imatinib in vitro | MAE PGx 2.9066; PGx + CYP protéines 3.0037; PGx + activités probes 2.7629 | **Non**; quantités de peptides et activité CYP spécifique ≠ même chose |
| **Métabolomique** | Rotroff, 32 personnes, plasma GC-TOF prédose | **AUC metformine 0–24 et CLr réels** | Avec âge/sex/poids/créatinine : AUC MAE 2948.825; +384 2936.749; +155 noms 3146.421; +229 inconnus 2819.953; **aucun gain robuste**, pour CLr toutes les ablations étudiées sont moins bonnes | **Un véritable essai E3 existe, mais résultat négatif/non concluant**; absence de PGx et PopPK |
| **Multi-omique (RNA, protéines, métabolites)** | NCI-60 n47 lignées dans plusieurs laboratoires, source externe | **GI50 doxorubicine cellulaire = PD**, pas PK | Type tumoral + RNA 0.334666, +RNA+protéines 0.341053, +RNA+protéines+métabolites 0.358124 MAE, donc fusion pénalisante dans ce modèle fixe | **Non**, aucun endpoint pharmacocinétique humain |
| **Protéformes / PTM d'albumine** | El Balkhi–Berrah et al. n254 individus, LC-HR-MS top-down | Stade de fibrose hépatique / performance versus FIB-4 | Preuve pertinente que **qualité protéique ≠ quantité totale** | **Non**, preuve **biomarqueur pathologique**, pas mesure d'un médicament ni prédiction fu/AUC |

**Important :** aucune cohorte source actuellement qualifiée ne possède **sur les mêmes individus** transcriptome + protéome + métabolome circulants + génotype fonctionnel complet + série de concentrations plasmatiques PK humaines. Il est donc **scientifiquement impossible à ce jour** de classer 4 couches pour leur bénéfice clinique relatif sur le même endpoint PK.

## Plan de comparaison obligatoire par jeu de données source

**Endpoint préspécifié selon le médicament et le type d'échantillon :**
- Tacrolimus et médicaments CYP3A/P-gp oraux : C0/AUC/CL/F seulement si données suffisantes; un C0 unique ne sépare pas CL et F, ni le V/F sans hypothèses. Comparateur PopPK existant + dose + temps + hématocrite, DDI, âge, etc. et **génotype réellement connu**.
- Metformine : AUC/Cmax/Tmax, CLr, signature OCT1/2/MATE et fonction rénale source; exclure variables calculées **à partir du même endpoint cible** des prédicteurs. **Ne pas utiliser protéformes albumine comme mécanisme de transport metformine par défaut.**
- Médicaments fortement albumine-liés : fu, C libre, AUC libre, clairance non liée et binding assay, protéformes HSA/PTM **top-down ou bottom-up PTM localisées** si effectivement dosées.

**Modèles à comparer sur exactement mêmes sujets et mêmes plis :**
`B0` = PopPK+covariables cliniques; `B1` = B0 + PGx complète si mesurée; `T` = B1 + transcriptomique; `P` = B1 + protéines (abondance), avec **PTM/protéformes sous-couche P-PTM séparée**; `M` = B1 + métabolites; `TP`, `TM`, `PM`, `TPM` combinaisons **uniquement lorsque réellement co-mesurées**. `P-PTM` isolée et `P+P-PTM` afin de tester si la structure protéique ajoute au total protéique. Ajouter `P-functional` si mesure d'activité enzyme/transporteur, **ne pas l'étiqueter « protéomique »**.

Pour chacune : même sujet/endpoint, modalité et unité biologiques identifiées, prélèvement prédose, prévention de fuite, validation externe si source permet. Imputation/normalisation/PCA/SWEET/MultiOmics et sélection des voies **fold-local**. Générer erreurs MAE RMSE bias P20, incertitude, calibration, `delta` contre `B1`, et garde-fous mechanistic P-gp versus OCT.

**Pipeline déjà disponible :** `src/omics2pk/multimodal_pk_ablation.py::compare_modalities_pk` supporte les couches transcriptomics/proteomics/metabolomics, comparateur clinical/PGx conditionnel, Pca2, Ridge, combinaisons et plis par patient/site. **Cette capacité logicielle n'est pas une validation à trois omiques sur un jeu de données clinique**; la couche PTM peut être un sous-type du protéome avec `assay_platform=top_down_proteoform` et `modified_residue`, `proteoform_identity_level`. Les PTM de liaison doivent être reliées à une mesure expérimentale fu ou concentration libre pour une preuve PK.

## Critère de succès
Gain reproductible **out-of-sample** sur PK humaine observée **au-delà** de clinique+PGx+PopPK, avec mécanisme et tissu corrects, calibration/détection de faux positifs; éviter de sélectionner le meilleur parmi 8 bras après examen des mêmes plis. Le résultat négatif sur une cohorte n'invalide pas la couche sur un autre médicament/tissu.

Sources primaires :
- [El Balkhi–Berrah 2026](https://doi.org/10.1038/s41598-026-57614-y)
- [Rotroff 2016](https://doi.org/10.3389/fphar.2016.00135)
- [Wegler 2021](https://doi.org/10.1021/acs.molpharmaceut.1c00053)
- [Tacrolimus PK / P-gp](https://pmc.ncbi.nlm.nih.gov/articles/PMC11235362/)
- [Metformin transport OCT / MATE](https://pmc.ncbi.nlm.nih.gov/articles/PMC5709222/)
- [Albumin drug binding and modifications](https://pubmed.ncbi.nlm.nih.gov/28526195/)
- Internal source-verified notes: `docs/HLC_ORIGINAL_400_LIVER_CYP_RNA_GENOTYPE_RESULTS_20261018.md`, `docs/WEGLER_SOURCE_STRUCTURE_RESULTS_20261016.md`, `docs/ROTROFF33_ORIGINAL_NAMED_UNKNOWN_PATHWAY_PK_20261104.md`, `docs/NCI60_THREEOMICS_ORIGINAL_GI50_RESULTS_20261104.md`.
