# Omics2PK — deux expériences de falsification : réseaux et fonction

**État 2026-11-04 : expériences exécutées sur des données publiques originales, exploratoires.** Aucun résultat n'établit un gain PK clinique du réseau, un effet causal des gènes ni une fonction de transport des pics inconnus.

## 1. Contrôle négatif du réseau biologique STRING (17 gènes NCI-60)

Question : deux associations entre les 10 gènes SHAP issus du modèle de réponse à la doxorubicine (GI50) et 7 gènes pharmacologiques constituent-elles un signal spécifique, ou sont-elles attendues compte tenu du petit réseau observé ?

Méthode implémentée : `src/omics2pk/network_degree_null.py`, appelée par `scripts/audit_nci60_string_shap_mechanistic_network.py`. Énumération **exacte** des affectations SHAP/prior, en conservant (i) les 17 identités des nœuds et les **10** arêtes fonctionnelles STRING v12.5 score≥0,4, (ii) le nombre total de gènes dans chaque rôle et (iii) le nombre d'étiquettes SHAP **à chaque degré de connectivité**. Aucun ajout de voisins, aucune nouvelle interaction, aucun échantillon patient. Les identités des gènes SHAP et des gènes pharmacologiques sont post-sélectionnées à l'origine : **l'échangeabilité statistique requise pour un test confirmatoire n'est pas établie**.

**Source réellement exécutée :** [CI NCI-60 nested #41926789502](https://github.com/rberrah/Omics2PK/actions/runs/41926789502), incluant l'appel réseau original.

| Métrique | Observé | Référence conditionnelle |
|---|---:|---:|
| Ponts gène SHAP ↔ gène pharmacologique | **2** | Moyenne **3** |
| Répartitions à au moins 2 ponts | — | **10/10** |
| Somme scores STRING des ponts | **1,284** | Moyenne **2,2602** |
| Répartitions à score cumulé ≥ observé | — | **9/10** |

**Interprétation :** aucun excès de ponts n'est constaté dans ce réseau sélectionné. Le nombre limité de réaffectations possibles (**10**) et la sélection préalable des gènes interdisent de conclure sur un enrichissement statistique, encore moins sur causalité ou mécanisme pharmacologique. Une vraie comparaison *avec le transcriptome/protéome mesuré complet*, appariée sur le degré STRING global, la détection, la variance et les conditions expérimentales, reste **non réalisée**. Les associations fonctionnelles STRING ne sont ni des preuves d'interactions physiques directes, ni des activités de transport mesurées.

## 2. Fonctions présupposées versus panels témoins appariés (NCI-60, PD in vitro)

Question : des **gènes dont la pharmacologie de la doxorubicine est connue indépendamment** améliorent-ils davantage la prédiction de la réponse que des groupes témoins de gènes simplement mesurés avec des propriétés analytiques similaires ?

Méthode : `src/omics2pk/mechanism_matched_control.py` et
`scripts/benchmark_original_nci60_mechanism_vs_controls.py`. Véritables **47 lignées** appariées sur leur identité biologique entre sources mixOmics (GI50 de doxorubicine) et IntLIM (ARN), mais prélèvements/laboratoires non contemporains. Les sept gènes cibles définis sans utiliser ces résultats sont `TOP2A`, `TOP2B`, `ABCB1`, `ABCC1`, `ABCG2`, `CBR1`, `CBR3`.

Pour chaque pli de validation externe (3 répétitions × 5 plis), le témoin **type tumoral seul**, le témoin + **sept ARN pharmacologiques** et **64 procédures témoins** de sept ARN aléatoires sont ajustés sur le même entraînement. Chaque ARN témoin est choisi dans les 17 987 signaux mesurés avec détection et variance proches des ARN pharmacologiques **sur le pli entraînement seulement**. Ridge `alpha=10`, transformation/imputation entraînement seulement. Les **47 lignées** (et non 141 sujets indépendants) reçoivent des prédictions hors entraînement.

**Source réellement exécutée :** [CI originale #41927200594](https://github.com/rberrah/Omics2PK/actions/runs/41927200594), tests synthétiques et source réelle validés.

| Modèle | MAE GI50 originale (plus faible = mieux) |
|---|---:|
| Type tumoral seul | **0,347203** |
| Type tumoral + 7 gènes pharmacologiques | **0,360858** |
| 64 procédures témoins, 5e / médiane / 95e percentile de MAE | **0,333669 / 0,353923 / 0,374083** |

`Delta MAE` panel versus référence tumorale : **+0,013655** (dégradation). **48/64** procédures témoins (75 %) ont une MAE plus faible que le panel pharmacologique.

**Interprétation :** ce **panel fixe d'expression ARN** ne surpasse pas la référence du type tumoral dans cet essai; les témoins appariés montrent que sa précision n'est pas spécifique. Les 64 tirages ne sont pas 64 validations biologiques indépendantes et ne fournissent pas une p-value confirmatoire. Ce résultat ne prouve pas qu'ABCB1 ou TOP2A sont fonctionnellement inutiles ; l'ARN n'est ni un dosage de protéine active, ni une mesure de transport, ni une activité enzymatique.

## 3. Chaîne signal métabolomique → structure → fonction → PK

**Source réellement exécutée :** [CI Rotroff–BinDiscover #41927200617](https://github.com/rberrah/Omics2PK/actions/runs/41927200617). Sur **229** bins GC-TOF originaux inconnus, **4** ressemblances d'identifiants d'accession sont possibles dans les MSP publics figés, mais ne constituent **aucune identité chimique validée**. Les trois pics hautement classés par SHAP (`222169`, `221574`, `218765`) n'ont **aucune** correspondance d'accession dans ce référentiel. Les trois sont qualifiés explicitement `ABSTAIN_ORIGINAL_STUDY_EI_RI_CROSSWALK_REQUIRED`. Aucune attribution à OCT2, MATE, CYP3A, P-gp ou à une voie métabolique n'est permise à ce stade.

Pour passer à un test de mécanisme, il faut **spectre EI réellement mesuré dans l'étude + indice de rétention / dérivatisation + identifiant compatible**, puis identification chimique confirmée et **relation substrat-transporteur/enzyme sourcée**; un dosage fonctionnel propre à l'essai et un modèle PopPK indépendant seraient nécessaires pour conclure à l'amélioration clinique.

## Conclusion et critères d'avancement

- Les expériences réalisées **ne valident pas** l'idée qu'une sélection de réseau ou des gènes connus améliorent automatiquement la généralisation du modèle.
- Une comparaison de réseau de fond sur **tous les gènes effectivement mesurés** et un vrai test de covariable enzymatique/transporter vers **PK humaine** restent à réaliser lorsque les données appariées, les spectres et les fonctions seront disponibles.
- Le résultat nul est conservé. Il motive l'utilisation prudente des connaissances physiopharmacologiques et la sélection validée hors échantillon, plutôt qu'un enrichissement de voies post hoc.

**Formulation pour collègues non spécialistes :** « On cherche quelques mesures biologiques réellement utiles pour prévoir la réponse à un médicament. Dans nos premiers essais, même des gènes connus et un réseau qui semble cohérent ne prédisent pas mieux à eux seuls. Il faut vérifier les molécules mesurées et tester sur d'autres patients avant d'en faire un outil pharmacocinétique. »
