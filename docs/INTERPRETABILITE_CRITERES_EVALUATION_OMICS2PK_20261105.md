# Omics2PK — interprétation mécanistique et choix des critères d'évaluation

**5 novembre 2026 · étude exploratoire · PR #137.**
L'objectif n'est pas seulement de prédire une concentration ou un effet, mais de
sélectionner le niveau omique utile **sans inventer un mécanisme** et en respectant
l'influence de variants authentiques.

## 1. Ce qui est biologiquement interprétable

L'ordre conceptuel est :
**génotype/variant → épissage et ARN → isoforme/protéine intacte →
activité enzymatique ou de transport → exposition PK ou effet PD**.
Les concentrations de métabolites sont des mesures contextualisées du
phénotype, **pas nécessairement des produits directs d'une seule enzyme**.
Une mutation de perte de fonction peut rendre une quantité d'ARN élevée
compatible avec une activité enzymatique très faible (épissage/NMD,
protéine tronquée, absence de site catalytique), mais cela demande une
annotation originale de la variante et des mesures fonctionnelles.

### Échelle de preuve à respecter

| Étape | Information effectivement vérifiable | Ce qu'on ne peut pas conclure |
|---|---|---|
| Génomique | Variant appelé, génotype/zygosité si source le permet, profondeur/qualité, localisation et conséquence sourcées | « Mutation » binaire CellMiner = CYP3A5 *3/*6/*7 ou activité enzymatique |
| Transcriptomique | Sonde ARN source, position/exons, expression sur échantillon identifié | ARNm = protéine active ; pseudogène/transcrit de stop = activité |
| Protéomique | Accessions, peptides détectés, isoformes et couverture si données disponibles | Protéine mesurée = pleine longueur, active et de bon compartiment |
| Métabolomique | Identité de composé validée selon méthode, matrice et temporalité | Pic inconnu = substrat / activité CYP / mécanisme de transport |
| Réseau | Associations STRING sourcées ; identités gène/chemical et types d'arêtes explicites | Association ou importance SHAP = causalité ou changement de CL |
| Pharmacocinétique | Même individu, vrais temps/doses, C0/C(t)/AUC/CL, génotype et omiques identifiés | Corrélation GI50 sur lignées = validation PK chez des patients |

Un réseau mécanistique doit distinguer explicitement une
**mesure**, une **association fonctionnelle**, une **hypothèse** et
une **relation causale établie**. Les liens de mécanismes pharmacologiques
n'entrent que s'ils sont sourcés, typés, et reliés à l'assay de
substrat/transport approprié. Aucune conversion automatique de SHAP,
de coefficients de régression ou de STRING en une relation causale.

### Deux niveaux complémentaires d'explicabilité déjà codés

1. **Explicabilité prédictive hors échantillon :**
   retirer/recombiner une modalité, réajuster uniquement sur l'entraînement,
   comparer les erreurs sur les mêmes lignées tenues à l'écart.
   Ces comparaisons indiquent si une vue apporte de l'information *au
   modèle estimé* ; ce ne sont pas des effets interventionnels. Des
   contre-exemples avec signal supplémentaire dégradant la prédiction
   sont retenus. Une attribution SHAP en présence de covariables
   corrélées ne quantifie pas directement le gain marginal en prédiction.
2. **Cohérence biologique directement observable :**
   `biological_cascade_audit.py` vérifie les sept gènes directs de
   la pharmacologie de la doxorubicine, si chaque gène possède réellement
   une mesure de mutation, d'ARN, de protéine SWATH et, lorsqu'elles
   existent, les associations d'abondances ARN/protéine en vrais
   identifiants de lignées concordantes. Il rapporte seulement des
   corrélations de Spearman **descriptives**, les nombres de mesures
   originales conjointes, les sous-groupes TP53 suffisamment représentés,
   et les **absence de preuve** (protéines absentes, mutations non
   génotypées, métabolites de structure non confirmée, flux non mesuré).
   Les omiques NCI-60 proviennent de **plusieurs études et préparations** :
   même lignée n'est pas le même prélèvement. Ce n'est pas une mesure de
   cascade intracellulaire au même instant.

**Garde-fou crucial :** le PCA3 de l'expérience multiomique mélange des
centaines de signaux. Il peut être évalué comme représentation prédictive,
mais **ses composantes ne sont pas des voies biologiques** et leurs
coefficients ne permettent pas d'attribuer directement un mécanisme à
ABCB1 ou à un métabolite. La sélection de signaux SHAP et leur
appariement à des voies doit être entraînement-seulement, confirmée par
ablation hors échantillon, et confrontée à un réseau témoin approprié.

## 2. Pourquoi MAE ne doit pas être le critère unique

Le MAE répond : « de combien me trompe-je *en moyenne dans l'unité
du critère* ? ». Il est cohérent si le coût de l'erreur absolue est
presque linéaire. **Il ne détecte pas séparément** :
- une poignée d'erreurs très graves (RMSE, quantile 90 % des erreurs
  absolues par individu) ;
- une sous- ou surestimation systématique (biais signé, pente de
  calibration observé~prédit) ;
- la qualité du classement des lignées sensibles/résistantes
  (Spearman par lignée biologique, pas sur 3 pseudo-répétitions) ;
- le bénéfice clinique d'une décision de dose ou le respect d'un
  intervalle thérapeutique ;
- l'incertitude/calibration de véritables **intervalles** de prédiction.

### Critères recommandés par contexte d'utilisation

| Donnée ou décision | Primaire proposé | Indispensables en secondaire |
|---|---|---|
| NCI-60 doxorubicine **GI50 traité par mixOmics** | MAE sur la valeur source, si la question est une valeur exacte | RMSE, Spearman (utile si classement), biais signé, P90 erreur, R² prédictif, coefficient de concordance, calibration descriptive, MAE par TP53 |
| Activité CYP **spécifique**, mesurée et quantitative | Erreur sur activité originalement vérifiée ou log-activité si préspecifié | RMSE, biais, comparateur génotype, diagnostics par variant, intervalles/calibration |
| PK humaine AUC/CL/C0/Cmax **strictement positive sur échelle vérifiée** | Erreur logarithmique absolue ou géométrique, si l'erreur multiplicative et les faibles/fortes expositions ont un poids clinique comparable | Biais en fold, RMSE(log-ratio), fraction dans un seuil fold *préspécifié*, erreurs extrêmes et sous-populations |
| PopPK/MIPD dose individuelle | Prédiction indépendante de concentration/exposition **et coût de décision clinique préspécifié** | Biais, précision, calibration, VPC/pcVPC, NPDE, couverture d'intervalles, décision de dose/target attainment |

**Exemple du danger du MAE seul :** un modèle qui se trompe un peu
sur presque toutes les lignées et un autre presque parfait sur la
majorité mais catastrophique sur une lignée peuvent avoir un MAE
proche ou inversé et un RMSE très différent. Le classement Spearman
peut aussi contredire le MAE. **Le meilleur indicateur dépend donc
de la décision que l'on cherche à améliorer**, pas du plus petit
nombre à posteriori.

### Contrats implémentés

- `src/omics2pk/predictive_diagnostics.py` :
  MAE, RMSE, biais, médiane/P90 erreurs par **lignée originale**,
  Spearman et Pearson sur **moyenne des prédictions répétées par
  lignée**, CCC, R² prédictif, pente/intercept de calibration
  *descriptive*, bootstrap des différences de MAE par lignée
  (1 999 tirages, sans réajuster les modèles).
- Pas de comptage des 3 répétitions comme 141 personnes indépendantes ;
  bootstrap apparié par identifiant biologique, non un test
  confirmatoire si modèles/variables ont été sélectionnés a posteriori.
- `multiplicative_pk_diagnostics()` n'est autorisée que
  pour `AUC`, `CL`, `C0`, `Cmax`, avec mesures et prédictions
  **originales, positives et vérifiées**. Une GI50 traitée, un
  phénotype enzymatique indirect ou un effet tumoral n'est pas une PK.
  Le seuil fold est un ratio **symétrique multiplicatif** : un seuil 1,20
  correspond à [1/1,20 ; 1,20], et **n'est pas** identique
  à P20 calculé comme `abs((prédit-observé)/observé) <= 0,20`.
- Aucun nouveau diagnostic ne modifie automatiquement
  le classement/panel issu de la **CV interne**. L'accord ou le
  désaccord des métriques est seulement documenté après prédiction.

### Références méthodologiques

- Chenel et al., *Best Practices in Physiologically Based Pharmacokinetic
  (PBPK) Modeling*, CPT:PSP 2026 : vérifier précision **et biais**,
  goodness-of-fit, intervalles, exposition (AUC/Cmax), validation externe ;
  [doi:10.1002/psp4.70295](https://doi.org/10.1002/psp4.70295).
- Berezowska et al., *Recommended approaches for integration of population
  pharmacokinetic modelling with precision dosing in clinical practice*,
  BJCP 2025 : GOF, VPC, NPDE, vérification de la précision et
  des incertitudes en situation d'usage clinique ;
  [doi:10.1111/bcp.16335](https://doi.org/10.1111/bcp.16335).
- [Article sur validation clinique externe et calibration, BMJ 2024](https://www.bmj.com/content/384/bmj-2023-074820.long).

### Validation et limites

Les sorties NCI-60 ne prouvent pas que la source GI50 a été mesurée
en unités pharmacologiques comparables à l'AUC ou sur une échelle
permettant de récupérer un ratio de concentrations. **Nous
conservons l'échelle numérique de la réponse source** et ne rapportons
pas automatiquement un MAPE, « % d'erreur de GI50 », dose, ou fold
sans qualification des unités/transformations.

Les résultats originaux ne sont des **preuves** ni de médiation
mutation→protéine→activité, ni d'amélioration de la PK clinique.
Une réplique sur **un même patient avec génotype germinal fonctionnel,
ARN/protéines/métabolites et PK observée** est nécessaire pour répondre
à la question centrale d'Omics2PK.

## Résultats réellement exécutés : métriques et interprétation biologique

Source originale indépendante et tests numériques :
[workflow NCI-60 #41974110167, SUCCESS](https://github.com/rberrah/Omics2PK/actions/runs/41974110167).
Toutes les métriques ci-dessous sont calculées sur les **47 mêmes lignées**
biologiques, avec trois prédictions hors entraînement par lignée (et
non 141 lignées indépendantes). L'endpoint est la **réponse GI50
doxorubicine traitée par la source**, sans conversion arbitraire
en unité de concentration.

### MAE, RMSE et classement ne répondent pas exactement à la même question

| Même protocole source | MAE | RMSE | Spearman par lignée | P90 erreur absolue par lignée |
|---|---:|---:|---:|---:|
| Type tumoral + ARN, **sans signal de mutation** | **0.308078** | **0.379938** | 0.479648 | **0.558996** |
| Type tumoral + mutations somatiques + ARN | 0.313205 | 0.395677 | **0.509366** | 0.670756 |
| Type tumoral + mutations + ARN + protéines | 0.328572 | 0.406734 | 0.502081 | 0.640206 |
| Type tumoral + mutations + ARN + métabolites | 0.337799 | 0.410918 | 0.463691 | 0.611291 |
| Type tumoral + mutations + trois omiques | 0.349732 | 0.426292 | 0.439292 | 0.630521 |

**Conclusion concrète :** le MAE, le RMSE et le P90 privilégient
l'ARN **sans les flags somatiques CellMiner**, alors que Spearman privilégie
numériquement l'ARN **avec ces flags**. Cela ne prouve **aucun**
avantage de classement statistiquement établi : la différence
de Spearman est petite et l'expérience n'a pas de cohorte externe.
Le choix d'un critère doit précéder la comparaison.

Au sein des **huit bras possédant tous les signaux de mutation**,
l'ARN seul est premier en MAE, RMSE et Spearman, ce qui confirme
la robustesse de sa performance **dans ce protocole seulement**.
Le biais moyen prédit−observé pour mutation+ARN est 0.001762 :
**un biais moyen faible ne suffit pas** à garantir une bonne
prédiction individuelle. Sa concordance CCC est 0.424119,
R² hors échantillon 0.231945, pente de calibration
**descriptive** 1.044527, et la P90 erreur 0.670756.

**Contrastes MAE appariés par lignée (bootstrap descriptif 1 999,
sans refit ni correction de multiplicité) :**
- Mutation+ARN − mutation seule : **−0.067973**, intervalle
  exploratoire **[−0.115179 ; −0.023317]**.
- Trois omiques+mutation − ARN+mutation : **+0.036527**,
  intervalle **[+0.003506 ; +0.070812]**.
- ARN+mutation − ARN sans mutation : **+0.005127**,
  intervalle **[−0.026375 ; +0.036135]**.
  Ce dernier est compatible avec un bénéfice comme avec une
  dégradation ; ne pas déclarer que les mutations aggravent
  **significativement** le résultat.

### Vraie cohérence inter-omique : même gène, source mesurée, pas hypothèse

Sur les **47** lignées d'origine, les deux assays sont disponibles
pour trois gènes directs pharmacologiques :

| Gène source | Corrélation ARN vs protéine SWATH (Spearman) | Ce qu'elle peut dire |
|---|---:|---|
| **CBR1** | **0.624** | Les abondances ont un lien inter-assay positif ; la corrélation ARN–GI50 est pourtant faible (−0.120), tout comme protéine–GI50 (−0.135) |
| **TOP2B** | **0.304** | Couplage modéré observé, pas preuve d'activité topoisomérase |
| **ABCB1** | **0.196** | Faible concordance entre ARN et abondance de peptides, malgré une corrélation ARN–GI50 de −0.403 et protéine–GI50 de −0.076 |

Tous ces coefficients sont **exploratoires, sans p-values confirmatoires,
non corrigés pour multiplicités ou tissu/lot**, et produits de
différents laboratoires/préparations. Une relation ARN↔protéine ne
suffit pas à prouver le transport de doxorubicine ou une modification de CL.

Sur **TOP2B**, le coefficient ARN↔protéine est 0.044 chez 24
lignées TP53 nonpositives et 0.741 chez les 23 positives.
C'est une **hypothèse exploratoire**, potentiellement influencée
par lots, types tumoraux ou hasard ; pas une interaction TP53
causale démontrée. L'ARN et la protéine CBR1 sont davantage
corrélés (0.624), mais les mutations CBR1 ne sont **pas
assayées** dans CellMiner. Le test direct
mutation CBR1 → ARN → protéine est donc **impossible** ici.

Les gènes TOP2A, ABCC1, ABCG2, CBR3 n'ont pas tous de
protéines dans le SWATH original. Pour ABCB1, seulement **4**
lignées ont un signal de mutation : effet fonctionnel de cette
mutation non estimable avec ce n. Les six métabolites « redox »
sont des **noms de source**, sans validation indépendante
de leur identité MSI-1, ni GSSG/NADPH/flux pour fermer la voie ;
ne pas attribuer une activité enzymatique à une abondance isolée.

La chaîne mécanistique gagne à être un **graphe des preuves
et des absences de preuve**, séparé des performances
prédictives. Une vraie validation sur PK humaine
nécessite le génotype fonctionnel complet, l'assay
protéique/activité, la métabolomique qualifiée et
les concentrations observées chez les **mêmes sujets**.

**Statut scientifique :** résultats descriptifs et réels,
**pas** validation d'une médiation biologique, d'une
mutation pharmacogénétique germinale CYP3A5 ni
d'une stratégie de dose pour le tacrolimus.

## Sensibilité prospective du choix des omiques : MAE ou RMSE en CV interne

**Reproduction réelle passée :** [47 lignées, originales, 17/17 CI vertes
sur le commit du code #41975269046](https://github.com/rberrah/Omics2PK/actions/runs/41975269046).
Dans chacun des plis externes, les deux stratégies évaluent exactement les
mêmes prédictions en validation interne : chacune choisit la plus petite
combinaison d'omics préservant 80 ou 90 % du gain d'erreur selon son
**critère propre** (MAE ou RMSE). La sélection est conditionnée par le
TP53 somatique source, mais **jamais optimisée sur l'outcome externe**.

| Même test original, avec mutations CellMiner | MAE finale | RMSE finale | Spearman de classement |
|---|---:|---:|---:|
| ARN fixé, sans sélection post hoc | **0.313205** | **0.395677** | **0.509366** |
| Sélection interne **MAE**, seuil 80 % | 0.329416 | 0.419691 | non requis pour la sélection |
| Sélection interne **RMSE**, seuil 80 % | 0.321363 | 0.406738 | 0.463691 |
| Sélection interne **MAE**, seuil 90 % | 0.330000 | 0.419943 | non requis pour la sélection |
| Sélection interne **RMSE**, seuil 90 % | 0.324323 | 0.408194 | 0.465888 |

Pour 80 % de conservation, **RMSE-sélection − MAE-sélection** :
- ΔMAE = **−0.008053**, intervalle bootstrap descriptif
  95 % [−0.019246 ; +0.001910] ;
- ΔRMSE = **−0.012953**, intervalle bootstrap descriptif
  95 % [−0.029722 ; +0.000104].

Les deux différences incluent zéro : **nous ne montrons pas**
de supériorité statistiquement établie du sélecteur RMSE.
Les deux méthodes restent moins performantes que le
**simple ARN fixé a priori**. Par rapport à ce bras fixe,
RMSE-selection 80 % ajoute ΔMAE **+0.008159**
(intervalle descriptif [+0.001485 ; +0.015420]) et
ΔRMSE **+0.011062** ([+0.002700 ; +0.021894]).
Les intervalles descriptifs ne sont pas corrigés pour
la multiplicité ni pour le développement du protocole au
contact d'analyses antérieures : aucune p-value confirmatoire.

L'effet apparent de l'addition de mutations somatiques sur le
**classement Spearman** après ARN (avec moins bon MAE/RMSE)
est **+0.029718**. Le bootstrap **apparié par les 47 lignées**
donne [**−0.105136 ; +0.183267**], englobant zéro :
pas de preuve que les mutations améliorent le classement.
Pour les erreurs numériques, l'ajout des flags à l'ARN
donne ΔRMSE **+0.015738**,
intervalle [−0.015017 ; +0.049662] ;
ΔMAE **+0.005127**,
intervalle [−0.026375 ; +0.036135].
Les valeurs exactes du résultat original ont été conservées,
sans rerégler les hyperparamètres après lecture de ces
performances.

**Conclusion méthodologique :** le critère doit être
défini selon la *décision* (valeur GI50 précise,
détection d'erreurs très importantes, classement de
réponse, exposition PK ou décision de dose) **avant** le test,
et non retenu pour favoriser le résultat le plus agréable.
Avec n=47 et des omiques inter-laboratoires, ni les
métriques secondaires ni les graphes biologiques ne
peuvent remplacer la validation externe.
