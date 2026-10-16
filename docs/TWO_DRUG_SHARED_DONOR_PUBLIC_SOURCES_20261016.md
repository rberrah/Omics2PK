# 2026-10-16 — Source publique : imatinib et tucatinib, mêmes donneurs hépatiques

## Découverte confirmée dans les tables primaires

- Tucatinib : DOI https://doi.org/10.1016/j.dmd.2025.100061 ; original Table 1 https://pmc.ncbi.nlm.nih.gov/articles/PMC13095627/ ; PRIDE PXD057282.
- Imatinib : DOI https://doi.org/10.1016/j.dmd.2025.100196 ; original Table 1 https://pmc.ncbi.nlm.nih.gov/articles/PMC12990092/ ; PRIDE PXD066482.
- Les deux articles utilisent les **mêmes 21 donneurs de microsomes hépatiques humains (HLM)**. Les identifiants originaux de donneurs, numéros de lots et quantifications CYP2C8/CYP3A concordent. Ne pas compter 42 donneurs indépendants.
- Le 150-donor pooled HLM est un mélange de référence, pas un sujet supplémentaire.
- Clairance microsomale imatinib quantifiable pour seulement **17 des 21 lots** : H0509, H0502, H0028 et H0525 ne présentent pas de déplétion mesurable. Ne jamais remplacer les données manquantes/censurées par zéro.
- Formation du N-desméthyl-imatinib : 21 rapports d'aires de pics mesurés. Ce n'est ni une clairance en µL/min/mg, ni une concentration plasmatique clinique.
- Le tucatinib dispose de 21 mesures de clairance dans la publication, avec réanalyse donneur-LODO antérieure du dépôt.

## Contrat logiciel ajouté

Le nouveau module src/omics2pk/matched_tki_metabolism.py exige deux tables locales et vérifie DOI, 21 IDs biologiques identiques, mêmes protéines et activités CYP par ID (tolérance d'arrondi), exclusion du pool, quatre clairances imatinib non quantifiées, puis quatre corrélations de Spearman publiées sur les 17 clairances mesurables. Aucun rapprochement par position de ligne ou groupe démographique n'est autorisé.

Pour chacun des trois résultats **séparés** (tucatinib CL, imatinib CL et imatinib N-desméthyl produit), comparaison donor leave-one-out : moyenne train seule, ridge alpha=1 sur deux protéines CYP, ridge alpha=1 sur deux activités enzymatiques, standardisation à partir du train uniquement. Les résultats numériques réels restent NON VALIDÉS tant que les deux tables originales n'ont pas été acquises et recontrôlées localement. Le dépôt public ne contient ni copies des tableaux des articles ni données donneur-level.

Commande locale :

    python -m scripts.validate_matched_tki_local /path/actual_tucatinib.csv /path/actual_imatinib.csv

Les colonnes imatinib exigées sont donor_id, row_kind, source_doi, imatinib_cl_app_ul_min_mg, ndesmethyl_imatinib_peak_area_ratio, cyp2c8_protein_pmol_mg, cyp3a_protein_pmol_mg, cyp2c8_activity_pmol_min_mg, cyp3a_activity_pmol_min_mg. La clairance manquante est une cellule vide ou le tiret « - », jamais zéro.

## Autre jeu public prioritaire

Sulfaméthoxazole, n=26 hépatocytes humains, NAT1/NAT2, formation NA-SMX : DOI https://doi.org/10.1016/j.dmd.2025.100195 et PRIDE PXD066786. Les mesures de protéomique source et activités enzymatiques sont publiées, mais il faut qualifier la table donneur→substrat→produit avant tout modèle. Ce n'est pas la même cohorte que les 21 donneurs TKI.

## Limites d'interprétation

**E2 exclusivement :** métabolisme mesuré in vitro, et non données patient dose–temps–C(t) ou PopPK. Les études TKI partagent les mêmes donneurs, ce qui constitue une validation *entre médicaments dans un même panel* et non une validation externe sur population indépendante. Une étude du potentiel mécanistique CYP est possible ; démontrer un gain prédictif clinique Omics→PK demeure une tâche indépendante et requiert les données réellement appariées.

Aucune imputation de cible, aucun fichier patient brut public, aucun effet CYP sur la clairance plasmatique prétendument mesuré.
