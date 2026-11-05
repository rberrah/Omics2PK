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
