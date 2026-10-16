# Degraeve CPT 2024 : audit du supplément original fourni — source locale autorisée

**16 octobre 2026 — source reçue et inspectée réellement, sans import de données patient dans GitHub.**

## Provenance primaire

Article : Degraeve AL et al., *Tacrolimus Pharmacokinetics is Associated with Gut Microbiota Diversity in Kidney Transplant Patients: Results from a Pilot Cross-Sectional Study*. *Clinical Pharmacology & Therapeutics* 115:104–115, https://doi.org/10.1002/cpt.3077 . Supplément éditeur `cpt3077-sup-0001-Supinfo.zip` fourni par l'utilisateur de façon autorisée. Empreinte SHA-256 **`03aab5deb292a03753073e493d7ae49896884be2cde22d5c2d1c9f3037488cc2`**, archive **3 905 302 octets**.

**Pas d'inclusion du ZIP dans le dépôt public, aucune exportation de mesures individuelles, de contenus d'objets GraphPad, d'identifiants SRA ni de données médicales.** Nouveau contrôle en lecture seule : `python -m scripts.audit_local_degraeve_supplement /chemin/autorise/cpt3077-sup-0001-supinfo.zip`. Le SHA est gelé ; les tests en CI utilisent **des ZIP synthétiques uniquement**, jamais l'archive reçue.

### Inventaire réel de l'archive

| Fichier original | Inventaire établi | Valeur pour Omics2PK |
|---|---|---|
| `2023-0506-s01.docx` | Méthodes expérimentales supplémentaires ; **aucun tableau brut clinique** | Critères de sélection, HPLC–MS/MS TAC, génotypage, protocole 16S, méthodes pSRBC et régression |
| `2023-0506-s02.pptx` | **8 figures** ; **7 objets OLE GraphPad Prism9** embarqués | Graphiques individuels/agrégés selon figure, mais **aucune table PK↔patient↔SRA exploitable et validée** |
| `2023-0506-s03.docx` | **3 tableaux** originaux : S1 15 lignes, S2 30 lignes (29 ASV annotés), S3 5 lignes | S1 associations cohortes ; S2 identité taxonomique ASV ; S3 modèles de **dose**, et non d'AUC |

## Contrôles méthodologiques confirmés dans le document

- **100 patients** ont fourni sang et selles ; **3** ont été exclus pour microbiote inexploitable et **4** pour non-respect des critères, laissant exactement **93** analysables. Critères initiaux : **18–80 ans**, BMI 17–35 kg/m², CYP3A5*3/*3 (non-exprimeurs), régime immunosuppresseur stable, sans antibiotiques durant les 3 mois précédents.
- Le sang TAC a été dosé par HPLC–MS/MS, les selles congelées rapidement, les amplicons **16S V5–V6** séquencés. C'est donc une omique **microbiote amplicon**, **ni transcriptomique humaine ni métabolomique quantitative**.
- Les covariables PK documentées sont **l'âge, l'hématocrite**, la dose, masse corporelle, temps post-greffe ; CYP3A5 constant, CYP3A4*22 et ABCB1 1199/3435 évalués sur **92/93** (une perte de matériel). LBP est disponible chez **46/93**, citrulline chez **92/93** : **ne pas imputer comme 93 valeurs observées**.
- Les taxons/ASV ont été présélectionnés par corrélations partielles de Spearman ajustées âge+hématocrite avec FDR BH, puis régression linéaire à sélection **stepwise** sur la même cohorte. Les résultats sont **exploratoires/in-sample**, pas des évaluations patients tenus à l'écart ; la sélection doit être reproduite **dans chaque pli d'entraînement** pour mesurer le bénéfice de scores microbiens hors-échantillon.
- Pour les ASV utilisés comme indicateurs présent/absent, seuil original **5 reads** ; il doit être documenté avec la profondeur séquencée et normalisation/compositionalité. ASV 1508 et 1982 sont fortement corrélés : attention à la redondance des prédicteurs.

### Ne pas confondre les deux modèles publiés

| Modèle original | Cible réellement étudiée | Covariables finales | R² de l'étude |
|---|---|---|---:|
| **Tableau supplémentaire S3** | **Dose quotidienne de tacrolimus**, pas AUC | Âge seul → hématocrite → ASV 1508 → ASV 664 | 18,5 % → 26,3 % → 34,9 % → **38,1 %** |
| **Tableau 2 de l'article principal** | **AUC/(dose/kg) estimée** (méthode LSS/ML, non NCA complète) | Hématocrite → ASV 1508 → âge → ASV 664 | 24,1 % → 40,1 % → 47,5 % → **52,1 %** |

Les pourcentages ci-dessus correspondent aux **ajustements du papier sur la cohorte analysée**. Ni `38,1 %` ni `52,1 %` ne sont des R² sur patients externes. Noter que la dose quotidienne peut être **ajustée par le clinicien selon les concentrations passées** : elle ne constitue pas un véritable phénotype PK intrinsèque indépendant de la décision thérapeutique.

### Les figures du PPTX ne constituent pas un jeu ML prêt à l'emploi

L'archive présente des barplots de familles/genres classés selon la dose, des graphiques de diversité, un nuage de points de co-abondance ASV et des comparaisons de marqueurs de barrière intestinale. Les **7 objets GraphPad Prism intégrés peuvent contenir les valeurs employées pour certaines figures**. Ces objets n'ont pas été interprétés comme des lignes patients en raison de leur format OLE propriétaire, du fait que plusieurs figures regroupent les observations et de l'absence de correspondance explicite **id patient ↔ PRJNA929904 alias ↔ dose ↔ C0/C1/C3**. Toute extraction éventuelle doit être qualifiée figure par figure, avec préservation des droits et sans pseudo-appariement par ordre des colonnes.

L'original n'apporte **aucune matrice de concentrations individuelles horodatées**, **aucune matrice ASV × patient individualisée**, **aucune correspondance explicite avec les 93 BioSamples ENA**, et aucun fichier d'estimations individuelles AUC avec indicateurs d'imputation. Il n'est donc **pas possible de réestimer un modèle Omics→PK** ou de valider NCA, PopPK/PBPK à partir de ce seul ZIP.

## Conséquences et nouvelle priorité scientifique

**E2/E3 :** les résultats du papier appuient biologiquement l'étude d'associations microbiote–TAC, mais seul un véritable jeu apparié permet de vérifier que les covariables microbiens expliquent la PK au-delà d'un modèle âge+hématocrite+clinique/PGx. La cohorte **CYP3A5 non-exprimeur seulement** ne permet pas d'estimer un effet entre génotypes CYP3A5.

**E3 :** avant de viser `CL/F`, `Ka`, `V/F`, obtenir toutes les heures de prélèvement/doses/formulations (et les observations originelles C0, C1, C3, si disponibles) ; ne pas utiliser comme vérité NCA une AUC estimée par un algorithme dont certaines valeurs ont été imputées. Evaluer `C0/(dose/kg)` sur vrais individus peut constituer un premier endpoint proxy, avec protocole de TDM bien défini.

**E4 :** le protocole est **transversal**, une seule photographie du microbiote : il ne peut démontrer des **covariables omiques variant chez un même individu dans le temps**. Pour ce volet, continuer les recherches Zhu (137 patients/693 visites) et Fannin (RNA sériel), et séparer clairement les profils moléculaires antérieurs versus post-dose.

### Demande de données minimales aux auteurs (usage privé local)

1. Tableau clinique original pseudonymisé `study_subject_id, patient_sampling_date_or_relative_time, tac_dose_mg, dosing_interval, formulation, dosing_time, C0_ng_ml, C1_ng_ml, C3_ng_ml, C0/C1/C3_is_observed, original_AUC_LSS, was_LSS_input_imputed, age, hematocrit, weight, time_post_KT, concomitants, CYP3A4_22, ABCB1_1199, ABCB1_3435`.
2. Tableau de correspondance original `study_subject_id ↔ BioSample (SAMN...), sample_alias IMM..., run SRR...` ; ne pas le déduire de l'ordre des SRA ou du rang dans la figure.
3. Matrice ASV original `study_subject_id × ASV_feature_id × count`, taille de bibliothèque, normalisation, tables d'assignations ASV (dont 1508/1982/664), métadonnées de batch ; préciser conditions d'utilisation et versions du pipeline.
4. Consentements/conditions d'accès, protocole de pseudonymisation ; héberger les données patient et la clé de correspondance **hors dépôts GitHub publics**.

Avec ces données, faire d'abord **un modèle de référence PK réel**, puis tester le gain prédictif microbiote avec sélection/normalisation apprise uniquement dans les plis et validation par patients indépendants. Aucune promesse de gain causal ou de pertinence posologique n'est faite en l'absence d'une source appariée.

**Statut : `SUPPLEMENT_QUALIFIED_AGGREGATES_ONLY / NO_PATIENT_LEVEL_PAIRED_PK`.** Le blocage concerne l'accès au **tableau clinique et sa clé de correspondance**, plus simplement l'accès aux méthodes supplémentaires.
