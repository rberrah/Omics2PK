# 17 octobre 2026 — SMX/NAT : évaluation numérique réelle du tableau 1 original (E2)

## Origine et traçabilité

**Source primaire :** Lanphier et al., *Drug Metabolism and Disposition* 2025, DOI [10.1016/j.dmd.2025.100195](https://doi.org/10.1016/j.dmd.2025.100195), PMCID [PMC12861138](https://pmc.ncbi.nlm.nih.gov/articles/PMC12861138/), Table 1, 26 hépatocytes primaires humains de **26 donneurs biologiques distincts**, formation NA-SMX à 30 min (pmol/min/million de cellules). Source de protéomique associée PRIDE PXD066786.

**Changement de statut par rapport à l'audit bloqué du 16 octobre :** les 26 lignes du tableau original ont été lues dans le contenu indexé publiquement de l'article le 17 octobre 2026, contrôlées et utilisées **localement et temporairement** pour une analyse source-originale. Les valeurs individuelles, IDs donneurs, copies des tableaux du journal ou matrices sources ne sont **pas** hébergés dans ce dépôt, conformément aux réserves de réutilisation de l'éditeur. La voie XML automatisée Europe PMC / PMC OAI demeure bloquée : le résultat scientifique a été obtenu **par transcription privée du tableau publié**, pas par succès d'accès JATS en CI. Ainsi, les résultats ci-dessous sont *recalculés et indépendamment vérifiés*, mais **ne sont pas reproductibles de bout en bout par CI publique** à partir d'un téléchargement brut à la date du rapport. Pour reproduire : utiliser une transcription locale autorisée, l'entrée `python -m scripts.validate_local_smx_nat PRIVATE.csv`, et confirmer les contrôles de corrélation d'article.

Empreinte de la **curation privée**, et **non du document éditeur** : SHA-256 `b69d215f8b05413c328eb5a2d2de7838d235dc1758df3fecbb5e89556e232823`, 26 lignes données et une ligne en-tête, UTF-8 CSV standard 12 colonnes, ordre des lignes original de Table 1, champs NA écrits `-`, floats avec précision source imprimée, colonnes dans l'ordre `donor_id,row_kind,source_doi,na_smx_pmol_min_million_cells,na_paba_pmol_min_million_cells,na_smz_pmol_min_million_cells,nat2_peptide_eflnshllpk_pmol_mg,nat2_peptide_dntdlvefk_pmol_mg,nat2_genotype_phenotype,age_years,sex,bmi`. Ni identifiants individuels ni données brutes ne sont exposés ici. La source publique contient bien **26** valeurs NA-SMX, **21** paires de substrats sondes NAT1/NAT2 et concentrations des **deux peptides NAT2** ; **20** génotypes NAT2 classés 9 Slow / 6 Intermediate / 5 Rapid ; **15** vrais donneurs disposent conjointement des quatre bras de modèle.

### Vérification indépendante des associations de l'article

Quatre corrélations de Spearman de source ont été recalculées indépendamment, en respectant les ex-aequo d'arrondi dans les valeurs publiées. Les petites divergences sont compatibles avec l'arrondi des valeurs de la Table 1 (le manuscrit utilise les données de laboratoire non arrondies).

| Données originales appariées | n | rho calculé depuis Table 1 | rho annoncé dans le texte scientifique |
|---|---:|---:|---:|
| NAT1 sonde NA-PABA vs NA-SMX | 21 | 0,572636 | 0,576 |
| NAT2 sonde NA-SMZ vs NA-SMX | 21 | 0,456790 | 0,459 |
| Peptide NAT2 **EFLNSHLLPK** vs NA-SMZ | 21 | 0,819038 | **0,823** |
| Peptide NAT2 **DNTDLVEFK** vs NA-SMZ | 21 | 0,730107 | **0,734** |

**Important correction d'audit :** la légende de Figure 8A inverse l'association peptide/rho par rapport au paragraphe Results §3.8. Les données numériques originales de Table 1 **confirment l'affectation du paragraphe §3.8** (EFLNSHLLPK ≈ 0,823 ; DNTDLVEFK ≈ 0,734). Le logiciel `SOURCE_ASSOCIATIONS` a été corrigé dans la présente livraison ; ne pas traiter la légende contradictoire comme une seconde source indépendante. Les quatre associations ont été reproduites à ±0,004 et passent le seuil logiciel de validation source ±0,035.

## Comparaison source-originale observée — leave-one-biological-donor-out

Deux implémentations indépendantes concordent jusqu'à la précision numérique utilisée : solveur fermé Ridge(alpha=1) avec centrage-réduction sur entraînement seulement, et **scikit-learn StandardScaler + Ridge(alpha=1)**. Tous les prétraitements sont effectués sur les n−1 donneurs. L'examen des résultats **n'a pas servi à optimiser alpha**, à sélectionner les variables après coup ou à imputer les cinq résultats NAT manquants.

Les erreurs sont en **pmol/min/million cellules**, cible : quantité NA-SMX formée à 30 minutes. La référence est la moyenne de la formation chez les n−1 mêmes donneurs éligibles pour chaque bras.

| Modèle | n biologiques complets | MAE moyenne train-only | MAE Ridge modèle | ΔMAE (modèle−référence) | Variation relative |
|---|---:|---:|---:|---:|---:|
| Age + sexe + IMC | 26 | 0,043046 | **0,045936** | +0,002890 | **6,7 % plus mauvais** |
| NAT1/NAT2 activité (deux sondes) | 21 | 0,042300 | **0,039033** | −0,003267 | **7,7 % meilleur** |
| NAT2 protéines (deux peptides) | 21 | 0,042300 | **0,040389** | −0,001911 | **4,5 % meilleur** |
| NAT2 phénotype inféré du génotype | 20 | 0,043284 | **0,039237** | −0,004047 | **9,4 % meilleur** |

Ces quatre lignes ne constituent **pas** une comparaison intermodèles équitable car elles utilisent n26/n21/n20. Par conséquent, le deuxième tableau est la comparaison directe conservatrice.

### Comparaison strictement appariée chez les mêmes 15 donneurs

| Modèle | n donneurs communs | MAE baseline commune | MAE modèle | ΔMAE |
|---|---:|---:|---:|---:|
| Age + sexe + IMC | 15 | 0,041914 | 0,042768 | +0,000854 |
| NAT1/NAT2 activités | 15 | 0,041914 | **0,037981** | −0,003933 |
| NAT2 protéines | 15 | 0,041914 | **0,043823** | +0,001908 |
| NAT2 phénotype inféré du génotype | 15 | 0,041914 | **0,039271** | −0,002643 |

La comparaison **à effectif et donneurs identiques** révèle que le modèle protéines NAT2 est légèrement **moins bon que la moyenne simple**, même si la comparaison chez n21 semble meilleure. Ce résultat renverse une conclusion naïve basée sur les quatre cohortes différentes. Les activités NAT sont les meilleures en MAE ponctuelle chez ces 15 donneurs, mais le classement est exploratoire.

### Incertitude, pertes et prudence

Dix mille tirages bootstrap **appariés des pertes LODO déjà calculées** (seed 20261017, valeurs de source arrondies) produisent les intervalles **exploratoires 95 % du ΔMAE**, séparément par bras :

- Donneurs complets de chaque bras : démographie n26 [−0,00433 ; +0,01006], activités n21 [−0,01198 ; +0,00533], protéines n21 [−0,00672 ; +0,00315], génotype n20 [−0,01331 ; +0,00609].
- **15 donneurs communs** : démographie [−0,00749 ; +0,00924] ; activités [−0,01641 ; +0,00956] ; protéines [−0,00096 ; +0,00501] ; génotype [−0,01032 ; +0,00590].

Tous ces intervalles recouvrent zéro : **aucune supériorité crédible démontrée** à cette échelle. Ce bootstrap n'est ni nested ni un refit complet et ne donne pas une preuve d'équivalence ou d'absence d'effet. Le statut scientifique est : **VALIDATION_IN_VITRO_NUMERICAL_EXPLORATORY_NON_SIGNIFICANT**.

## Interprétation et suites

- Ce jeu permet un **troisième médicament et deuxième famille enzymatique** NAT (distinct des 21 microsomes HLM CYP/tucatinib+imatinib), mais le système biologique, les variables et l'endpoint sont différents ; ce n'est pas une validation externe des modèles CYP ni une cohorte clinique PK.
- Le phénotype NAT2 inféré du génotype n'est pas un prédicteur fiable de chaque individu pour NA-SMX, malgré une MAE moyenne modestement favorable. Le papier rapporte six slow avec forte formation NA-SMX.
- Les peptides NAT2 ciblés sont quantitatifs, **pas** une protéomique globale et NAT1 protéines étaient BLQ. Ne pas interpréter l'activité NAT1 comme une preuve directe de mécanisme causal à elle seule.
- La preuve clinique Omics→PopPK exige encore données du même patient pré-dose omiques, dose précise, profil concentration−temps observé et suivi clinique. Aucune clairance plasmatique, AUC ni risque individuel de DILI n'est estimé.
- Prochain jalon de reproduction : accès autorisé à une table JATS/CSV intacte, exécution du CLI du dépôt et vérification du hash de curation ; ce rapport ne prétend pas que la CI publique a pu télécharger l'article automatiquement.

