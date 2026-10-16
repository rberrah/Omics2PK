# Protocoles E2 → NCA → E3 PopPK/PBPK → E4 omiques et PK longitudinales

## Question et hiérarchie des preuves

1. **E2 mécanistique in vitro** : protéines/RNA mesurés → activité CYP ou CLint microsomale originale, validée sur donneurs biologiques distincts. Les exemples EPA 100 foies / tucatinib 21 foies répondent à cette question, pas à la PK clinique.
2. **Passerelle NCA clinique** : partir de véritables doses et profils individuels de concentrations C(t), calculer **Cmax, Tmax, AUC0–t** (et AUCtau à l'état stationnaire seulement si l'intervalle de dose et l'historique complet sont connus). Vérifier la liaison même sujet, prélèvement omique antérieur à la dose/fenêtre prédite, unité, matrice plasma/sang total, formulation et moments de prélèvement. Comparer les données cliniques+PGx à l'ajout de RNA ou signatures omiques sur **les mêmes patients**. Utiliser les unités et fenêtres rigoureusement comparables : **AUC0-last avec dernier temps différent n'est pas une cible homogène**. Cmax et Tmax varient avec densité et timing d'échantillonnage, leur valeur discrète observée n'est pas l'optimum exact de la courbe PK réelle. Une amélioration de corrélation ne prouve pas un bénéfice prédictif.
3. **E3 PopPK** : évaluer les paramètres structurels (CL/F, V/F, Ka/MTT et éventuellement biodisponibilité F si identifiable ; pour IV CL, Vc, Q, Vp) par NLME/SAEM sur le **profil dose–C(t)** et non régression naïve des paramètres individuels calculés à partir des mêmes observations utilisées pour évaluer le modèle. Comparer le **modèle clinique+PGx de référence** à un effet biologique sur CL, absorption ou distribution **spécifié et sélectionné dans les plis d'entraînement**. Critères : erreur de prédiction concentration, AUC, calibration, résidus, incertitude, gestion des valeurs BLQ et extrapolation à des sujets externes. Les critères BIC soutiennent le choix mais ne remplacent pas le test prospectif.
4. **E3 PBPK complémentaire** : transporter des mécanismes issus de l'omique (enzymes, transporteurs, organes, dynamique physiologique) vers des paramètres de clairance intrinsèque/activité spécifiques, avec quantités protéiques, fraction libre, scaling hepatocyte/microsome→foie et phénomènes d'absorption. **RNA mesuré ≠ expression protéique ni CLint fonctionnelle par définition**. Prévoir incertitude et vérification clinique indépendante ; ne pas déduire une clairance systémique d'un log2FC.
5. **E4 longitudinal** : mesures omiques répétées réellement chez les mêmes sujets avec dates de dose et de concentration. Les signatures peuvent être des **covariables dépendantes du temps**, éventuellement influencées par le médicament. Les échantillons après dose peuvent expliquer une réponse biologique et prédire *un futur* intervalle de concentration, mais jamais l'AUC ou Cmax d'un intervalle passé. Séparer dynamique interindividuelle et intra-individuelle ; tenir compte de l'histoire complète de doses et des effets de formulation ou de greffe. **NewPK** pourra devenir un bras comparateur distinct uniquement après qualification de ces références PopPK/PBPK.

## Architecture minimale de données, sans additionner les cohortes par erreur

`MultiOmics` → `exchange/manifest.json` + `samples.csv` + `features.csv` + matrices source ; optionnel `individual_signatures.csv` issu d'une référence figée. Ces sorties ne doivent pas inclure implicitement de PK.

`Omics2PK` → lire ce **même** échange sur infrastructure locale, ajouter les profils originaux administrés/mesurés `PK_LONG_V1`, puis `NCA_TARGETS_V1` explicitement qualifiés. Champs NCA : sujet/occasion/dose/analyte/matrice, fenêtre, paramètres Cmax/Tmax/AUC0last, unité, nombre de concentrations mesurées, méthode d'intégration et empreinte source PK. Toute liaison requiert une origine temporelle *commune et documentée*, jamais construite en faisant correspondre des noms de groupes ou des cohortes étrangères.

`Pilier B / MIPD-Engine` → consommateur de covariables individuelles ainsi qualifiées ; modèle PopPK identifiable à partir du PK original, avec splitting patient/disponibilité longitudinale. `Pilier A` → prior mécanistique protéine/voie, candidat et incertitude. `Pilier C` → confirmation/abstention de la cohérence mécanistique et du sens de l'effet. **Ne pas intégrer le code de ces piliers potentiellement protégeables dans un dépôt éducatif public.**

## Validation en étapes

| Test | Minimum pour autoriser une conclusion |
|---|---|
| Omics seuls | SHA, identités, quantification, unité de temps réelle et comparaison RNA seul vs autres omiques sans pseudo-réplication |
| NCA + Omics | Association observée sur **mêmes personnes** avec vrais profils C(t), fenêtres comparables et covariable RNA mesurée avant fenêtre |
| PopPK + Omics | Modèle de référence clinique+PGx, effets PK identifiables, train/test patient-disjoint, prédiction externe/temporalement future |
| PBPK + Omics | Couplage biologique→enzyme/protéine→CLint→organe physio avec incertitude de translation et confrontation à des concentrations humaines |
| E4 | Répétition réelle de mesures moléculaires, échantillonnage prospectif sans anticipation de la PK future, validation intra/interpersonne |

**État au 16 octobre 2026 :** seuls E2 statique et des vrais jeux **PK seuls** sont disponibles et techniquement qualifiés dans le dépôt. Pas de jeu public documenté réunissant, sur les mêmes sujets, un dosage réel, C(t), RNA longitudinal et signatures multiomiques avec appariement exhaustif ; aucune amélioration PopPK/PBPK imputable à MultiOmics ne peut donc être affirmée. Le développement évite les faux « appariements » entre donneurs EPA/tucatinib et sujets Theoph/mavoglurant.

Documents connexes : [Audit MultiOmics](REAL_MULTIOMICS_OUTPUT_AUDIT_AND_PK_ADAPTER_GAPS_20261016.md), [Roadmap temporelle](LONGITUDINAL_PK_AND_DYNAMIC_OMICS_ROADMAP_20261016.md) et [jeux publics de PK](ORIGINAL_LONGITUDINAL_PK_SOURCE_AUDIT_20261016.md).
