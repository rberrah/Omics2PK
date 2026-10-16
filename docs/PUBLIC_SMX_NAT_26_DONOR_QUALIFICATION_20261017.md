# Sulfaméthoxazole (SMX) → NA-SMX : troisième médicament, cohorte PHH indépendante

**Statut 17 octobre 2026 : 26 donneurs humains dans la Table 1 primaire ; contrat source et tests synthétiques en place. Pas encore de résultats de prédiction sur les valeurs individuelles originales SMX.**

## Publication et dépôt public

Lanphier et al. Drug Metab Dispos (2025), DOI https://doi.org/10.1016/j.dmd.2025.100195 ; https://pmc.ncbi.nlm.nih.gov/articles/PMC12861138/ ; dépôt de protéomique NAT2 PRIDE PXD066786 (https://www.ebi.ac.uk/pride/archive/projects/PXD066786).

Le jeu original comprend des hépatocytes humains primaires (PHH) de **26 donneurs individuels**. La Table 1 comporte :

- **NA-SMX**, formation en pmol/min/million cellules à 30 minutes pour n26.
- **NAT1 et NAT2 activités**, substrats-sondes PABA → NA-PABA et SMZ → NA-SMZ, analysés chez n21.
- **NAT2 génotype/phénotype**, n20 : slow 9, intermediate 6, rapid 5. Ne pas inventer six génotypes.
- Deux peptides protéomiques quantitatifs **NAT2.EFLNSHLLPK** et **NAT2.DNTDLVEFK** en pmol/mg ; leur disponibilité exacte doit être qualifiée ligne par ligne.
- Âge, sexe, IMC ; les données manquantes sont déclarées par bras, non imputées.

Les PHH SMX ne sont **pas** les mêmes donneurs ni le même système que les 21 lots HLM des articles imatinib/tucatinib. Cela donne un mécanisme NAT et un troisième médicament, mais **pas** une validation externe du modèle CYP sur les mêmes caractéristiques.

## Protocole de comparaison fixé

Validation leave-one-original-donor-out, transformation train-only, régression Ridge à alpha=1 fixée, résultats en MAE de NA-SMX dans ses unités originales. Quatre bras, comparés chacun à une moyenne train-only des **mêmes donneurs complets** :

1. Age, sexe, BMI (clinique/démographique).
2. NAT1 et NAT2 : deux activités de formation des sondes.
3. Deux concentrations de peptides NAT2.
4. Phénotype génotypique NAT2 (slow, intermediate, rapid).

Les bras ne sont pas automatiquement comparables si leurs sous-cohortes diffèrent : n20 génotypé, n21 activité, éventuellement autre n proteomics. Le module annonce les effectifs, ne comble jamais les valeurs manquantes ou non quantifiables par zéro et refuse d'entraîner sur moins de dix donneurs complets. Aucun test dupliqué, puits technique ou pool ne devient un nouveau donneur.

## Implémentation sans redistribution de valeurs originales

- src/omics2pk/smx_nat_local_validation.py : contrôle DOI original, cardinalités 26/21/20, 9/6/5 classes NAT2, unités et rejet doublons.
- scripts/validate_local_smx_nat.py : CLI locale, source CSV originale autorisée.
- tests/test_smx_nat_local_validation.py : contrats et entraînement synthétiques, aucun donneur réel en CI.
- Exécution : python -m scripts.validate_local_smx_nat /private/original_NAT_table1.csv

Colonnes CSV locales : donor_id,row_kind,source_doi,na_smx_pmol_min_million_cells,na_paba_pmol_min_million_cells,na_smz_pmol_min_million_cells,nat2_peptide_eflnshllpk_pmol_mg,nat2_peptide_dntdlvefk_pmol_mg,nat2_genotype_phenotype,age_years,sex,bmi. Les valeurs non disponibles restent vides, NA ou ND. Sortie limitée aux nombres de donneurs, MAE agrégées et empreinte de la curation ; pas de reproduction du tableau éditorial dans le repo.

## Limites

La formation NA-SMX dans des cellules PHH est une **mesure fonctionnelle in vitro**, pas clairance plasmatique du SMX, AUC, taux de DILI ou dosage pharmacothérapeutique. Dans l'article, certaines classifications NAT2 « slow » ne correspondent pas à une formation faible de NA-SMX. Les associations sont internes à une petite étude. La comparaison E2 indépendante reste exploratoire et n'établit pas Omics → PopPK patient, laquelle requiert des doses et profils concentration–temps originaux appariés à la biologie du même sujet.
