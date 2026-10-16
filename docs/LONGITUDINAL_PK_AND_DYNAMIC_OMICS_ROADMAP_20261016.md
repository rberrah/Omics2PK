# Roadmap scientifique : PK longitudinale et omiques dynamiques — 16 octobre 2026

Objectif de thèse : dépasser les corrélations statiques CYP/clairance et relier **l’état biologique mesuré d’un individu aux profils PK C(t)**, puis à la PopPK avec covariables cliniques/PGx. La préparation des données se situe dans Omics2PK ; l’estimation PK et la qualification des mécanismes demeurent dans les piliers B/C conformément à l’architecture gelée.

## Progression vérifiable

| Niveau | Critères empiriques | Premier jeu public | État |
|---|---|---|---|
| L1 — PK après une dose | Vraies doses, concentrations mesurées à plusieurs temps et identités biologiques | R Theoph, 12 sujets × 11 prélèvements | Contrat PK_LONG_V1, audit sur données originales |
| L2 — PK à doses répétées et PK/PD | Historique de doses, C(t), formulation, OCC, censure et éventuellement réponse PD | nlmixr2data warfarin, mavoglurant, theo_md | Audit des formats, pas encore de NLME |
| L3 — omiques et PK longitudinales | Même patient : RNA/protéines/métabolites aux dates connues + doses et concentrations réellement mesurées | Fannin 2010 APAP (6 exposés, 3 placebo), potentiellement données autorisées tacrolimus/MPA | Données individuelles téléchargeables et correspondance non confirmées |
| L4 — démonstration PopPK | Comparaison référence clinique+PGx vs état RNA direct/SWEET et omiques parallèles sur mêmes patients et dates | Cohorte indépendante autorisée | Non qualifié |

## Contrat PK_LONG_V1

Source dans src/omics2pk/longitudinal_pk_contract.py. Un événement dose comporte l’identité, l’origine temporelle commune, la dose en **mg**, la voie, la formulation et sa provenance **observed/simulated**. Une mesure PK comporte la même identité, le temps de prélèvement, l’analyte, le type de matrice (sérum/plasma/sang total), la concentration dans l’unité déclarée et les statuts BLQ/LLOQ. Une mesure t=0 concomitante à une dose doit être étiquetée explicitement predose.

La dose antérieure « la plus proche » ne suffit **jamais** à reconstruire la PK d’un schéma multidose ; les doses précédentes doivent être connues et conservées. Si l’historique est partial/unknown, le contrôle ne le requalifie pas en complet.

Les prélèvements omiques deviennent **des états biologiques dépendants du temps**. Leur utilisation est verrouillée à la date de prédiction : uniquement une mesure antérieure, respectant une fenêtre de validité biologiquement préspécifiée. Un RNA post-dose peut être un effet du traitement et servir à prédire le futur, mais ne peut pas expliquer rétroactivement la concentration initiale.

Pour la validation, séparer les individus entiers entre entraînement et évaluation, puis distinguer les analyses intra-patient prospectives : pas de concentrations, omiques ou normalisation calculées à partir d’une visite future. Répétitions, clairance, BLQ, changement de traitement, Tmax, AUC0-t, concentrations résiduelles et occasion nécessitent un protocole PK explicite.

## Qualification des jeux

- **Theoph** : 12 sujets, 132 concentrations **réellement mesurées** après une dose orale ; source R datasets. Quantité orale originale (mg/kg) convertie explicitement avec le poids (kg) en mg. Certains échantillons avant la dose ont une concentration non nulle : valeur originale conservée, pas ramenée artificiellement à zéro.
- **theo_md** : 348 lignes, les observations après la première phase sont **simulées** ; très utile pour vérifier le code de doses répétées, mais impossible de revendiquer une validation clinique multidose réelle.
- **warfarin** : mélange PK et PD ; séparer DVID avant le calcul de la vraisemblance PK. Dans la version épinglée de l’objet RDA, 515 lignes ont été effectivement lues alors que la documentation publique annonce 519 ; ne pas inventer les quatre lignes manquantes, conserver un ticket de réconciliation.
- **mavoglurant** : 2 678 lignes documentées avec OCC/DOSE/AGE/SEX/WT/HT, à vérifier sur l’objet original et adapter sans confondre dose, MDV et DV.
- **PK-DB** : utile pour doses, profils et provenance ; API publique consultable, mais exports demandant une authentification, avec provenance individuelle ou collective selon étude.
- **Fannin et al., Hepatology 2010**, DOI 10.1002/hep.23330 : paracétamol 4 g, six sujets exposés et trois placebo ; véritables concentrations sériques mesurées dans les 12 premières heures et RNA sanguin/métabolites suivis jusqu’à plusieurs jours. Les signatures RNA notables vers 48 h surviennent **après** les concentrations PK précoces : pertinent pour dynamique PK→réponse omique/PD, pas pour prétendre qu’un RNA de 48 h prédit une AUC de 0–12 h. Matrices individuelles exactement appariables encore non vérifiées.

## Prochaines tâches scientifiques, en ordre

1. Valider d’abord les données longitudinales observées Theoph : formats, unités, qualité des profils et calcul NCA indépendant. Pas de simulation présentée comme observation.
2. Adapter ensuite des données multidose réelles (mavoglurant), puis PK/PD combinée (warfarin), avec séparation EVID/MDV/CMT/DVID/OCC et signalement des événements simulés ; qualification des données avant SAEM/nlmixr2/Monolix.
3. Rechercher un **même jeu** avec mesures moléculaires de référence et/ou répétées, et historique PK complet, sans matching probabiliste entre cohortes.
4. Tester les hypothèses à événements temporels avec covariables cliniques + PGx d’abord, puis RNA direct, SWEET qualifié, protéines/métabolites et combinaison. Si la cascade RNA→protéine→metabolite est une conséquence du traitement, tester causalité/retour et non seulement association.
5. Après validation empirique et examen IP/licence : cours FR/EN et outil de recherche à rberrah.github.io, sans données patients publiques ni prescription automatique.

Liens : https://search.r-project.org/R/refmans/datasets/html/Theoph.html ; https://nlmixr2.github.io/nlmixr2data/reference/theo_md.html ; https://nlmixr2.github.io/nlmixr2data/reference/warfarin.html ; https://nlmixr2.github.io/nlmixr2data/reference/mavoglurant.html ; https://matthiaskoenig.github.io/pkdb/api/ ; https://pmc.ncbi.nlm.nih.gov/articles/PMC2925683/ .


**Qualification originale mavoglurant après lecture RDA :** 120 sujets biologiques, **198 occasions de dose unique** (198 événements EVID=1, 2 480 concentrations EVID=0), dont **78 sujets évalués à deux occasions** ; aucun des 198 groupes sujet×occasion ne contient plus d'une dose. Il s'agit de profils PK IV répétés entre visites chez une partie des sujets, **pas** d'une série de doses successives au sein d'une même occasion. Le dataset A2121 (Wendling) concerne des infusions IV de 10 minutes ; la vignette nlmixr2 affiche la concentration plasmatique en **ng/mL**. Ne jamais coller OCC 1 et 2 sur un même axe TIME sans calendrier de visite. L'entrée multi-occasion sera la prochaine cible de l'adaptateur pour les modèles PopPK ; aucune omique appariée.
