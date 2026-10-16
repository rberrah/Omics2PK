# SMX/NAT : qualification d'une source originale publique et contrôle anti-biais

Date : 17 octobre 2026. Statut : **pipeline prêt et tests synthétiques verts ; aucun résultat original SMX sur les 26 donneurs n'a été produit**.

## Données publiées existantes

Publication Lanphier et al., *Drug Metabolism and Disposition* 2025, DOI 10.1016/j.dmd.2025.100195, PMCID PMC12861138, PRIDE PXD066786. Le tableau 1 public de l'article comporte 26 donneurs indépendants de PHH avec la formation quantifiée de NA-SMX, 21 avec activités NAT1/NAT2, 20 avec NAT2 phénotype inféré du génotype (9 slow/6 intermediate/5 rapid), et les concentrations de deux peptides protéomiques NAT2 quand quantifiables. Ces ensembles **ne sont pas interchangeables**.

Les 26 donneurs sont distincts des 21 microsomes des deux études CYP/imatinib/tucatinib, et le système enzymatique NAT est distinct des CYP. Ils ne constituent **pas** une validation indépendante du même modèle CYP.

## Sources et accès techniquement testés

Un workflow GitHub public, **Original SMX article XML source audit**, a effectivement vérifié trois voies officielles depuis un environnement Python 3.11 de CI, sans afficher de données donneur :

- **Europe PMC JATS fullTextXML** sur PMC12861138 : HTTP 500 au moment du test.
- **Page PMC avec demande de report=xml** : aucune version JATS exploitable n'a été retournée (HTML non accepté comme XML).
- **Service officiel PMC OAI-PMH de récupération JATS** : HTTP 400 au moment du test.

L'état réel est **SOURCE_UNAVAILABLE_ABSTAIN**, et non un benchmark empirique ayant réussi. Les tests synthétiques passent même si les sources distantes échouent : ce sont deux qualités de preuve indépendantes. Les droits de réutilisation et de fouille de textes réservés par l'éditeur doivent être respectés ; n'utiliser les voies officielles et/ou un exemplaire téléchargé légalement pour les analyses locales. La disponibilité du tableau **en lecture sur la page publique** ne vaut pas licence de republication de toutes les lignes.

Le service PMC Author Manuscript Dataset donne accès sous conditions au XML des manuscrits ; les textes ne peuvent être réutilisés que dans le cadre des licences de leur article. La source primaire reste https://pmc.ncbi.nlm.nih.gov/articles/PMC12861138/ . Un test HTTP 400 ne démontre pas à lui seul une interdiction juridique : il établit seulement l'absence d'une extraction autorisée et exploitable dans la CI.

## Deux chemins locaux disponibles sans recontacter les auteurs

1. À partir d'une **copie originale JATS XML obtenue légalement** : exécuter
    `python -m scripts.validate_smx_nat_original_jats /private/PMC12861138_original.xml`

   Le parseur vérifie DOI, PMCID, exactement un tableau d'origine, 26 identifiants, 13 colonnes à l'ordre connu, deux peptides NAT2, génotypes n20, activités n21 et **quatre corrélations publiées** (NAT1 et NAT2 vs NA-SMX ; deux peptides NAT2 vs NA-SMZ). Il s'abstient si les contrôles échouent.

2. Si une **copie de Table 1 est consultée et reportée localement dans un CSV**, sans redistribuer la reproduction intégrale du tableau :
    `python -m scripts.validate_local_smx_nat /private/lanphier_table1.csv`

   Contrat de colonnes : `donor_id,row_kind,source_doi,na_smx_pmol_min_million_cells,na_paba_pmol_min_million_cells,na_smz_pmol_min_million_cells,nat2_peptide_eflnshllpk_pmol_mg,nat2_peptide_dntdlvefk_pmol_mg,nat2_genotype_phenotype,age_years,sex,bmi`. Pour ce second chemin CSV, **la vérification indépendante des quatre corrélations publiées reste à ajouter à l'entrée CLI**. Aucune mesure originale n'est stockée dans GitHub.

Ces commandes ne sont pas encore des résultats publiés. Ne pas indiquer que des données ont été chargées si seuls des jeux synthétiques ont été testés.

## Évaluation prévue

Les modèles NAT sont entraînés **en leave-one-original-donor-out** avec standardisation et réglages uniquement sur les donneurs d'entraînement : référence moyenne sans covariable ; âge/sexe/BMI ; NAT1/NAT2 activités ; NAT2 peptides ; classification NAT2 génétique. Chaque bras possède une baseline entraînée sur ses propres cas complets, en plus d'une **comparaison supplémentaire sur une intersection commune stricte** aux quatre modèles. Si cette intersection présente moins de 10 donneurs, l'analyse comparative s'abstient.

Les quatre NAT/CYP associations servent de contrôle de transcription, pas de critère de performance du modèle. Un résultat négatif doit être publié sans modification opportuniste des covariables ou du score. La NA-SMX formation est un résultat métabolique **in vitro**, pas une clairance clinique, une aire sous la courbe, ni une toxicité réelle.

## Prochaine porte scientifique

Achever l'acquisition locale légitime du tableau complet, vérifier l'identité des 26 donneurs et reproduire les corrélations, calculer les performances LOOCV réelles et effectuer une sensibilité aux données manquantes et aux valeurs sous limite de quantification. Ensuite, utiliser le jeu NAT comme **un autre mécanisme et une autre cohorte**, sans conclure sur la prédiction clinique Omics→PopPK. Pour la clinique, il faut toujours des données mêmes patients avec omiques pré-dose, doses, temps et concentrations.
