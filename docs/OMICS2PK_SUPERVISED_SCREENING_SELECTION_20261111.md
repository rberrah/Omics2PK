# Screening supervisé et sélection mécanistique : audit, conception, premier essai

**Version 1.1.0 — 2026-11-11.** Niveau de preuve : **exploratoire, exécuté sur
11 cellules**, deux cohortes, six contrôles dont un contrôle positif
obligatoire. Module : `src/omics2pk/supervised_screening_selection.py`.
Script : `scripts/run_supervised_screening_selection.py`. Artefact :
`reports/data/supervised_screening_selection.json`.

**Indépendance.** Rien ici ne lit, n'importe ni ne dépend de Thesis, du
pilier A, de NetPK, NewPK ou MetaPK. Les bases de connaissances publiques sont
utilisées directement.

---

## 1. Audit — ce qui existait déjà, et ce qui ne valait pas la peine d'être réécrit

L'instruction était de réutiliser au maximum. Sept composants existants font
déjà une partie du travail demandé, et les réutiliser évite surtout de
réinventer des garde-fous déjà payés.

| Demandé | Existait déjà | Réutilisé tel quel |
|---|---|---|
| XGBoost | `hlc_real_xgboost_validation.py` — configuration **gelée** : 150 arbres, profondeur 2, lr 0,035, min_child_weight 8, subsample 0,85, colsample 0,9, lambda 10 | **oui**, à l'identique, pour que le bras de screening reste comparable à la comparaison HLC gelée et pour qu'aucun hyperparamètre ne soit réglé sur une issue |
| SHAP | rien | **pas nécessaire** — voir §2 |
| sélection parcimonieuse | `multiomic_layer_complementarity.parsimonious_subset` — règle de parcimonie avec plancher de tolérance, et qui **déclare déjà** que la performance du sous-ensemble retenu est biaisée optimiste | la **règle** est reprise ; la courbe performance–k est nouvelle |
| sélection sur le résidu | `hybrid_mechanistic_selection.partial_rank_strength` — corrélation de rang **partielle**, résidualisée sur la base | **oui**, c'est le sélecteur `S2` |
| sélection guidée par graphe | `mechanistic_graph_selection.select_panel` + `permutation_fdr_cutoff` + `degree_strata` | disponible pour l'étape suivante ; pas dans ce premier essai, qui n'a pas besoin de STRING |
| ablation hors échantillon | `variant_mechanism_omic_selection.ablation` (ajouté le 2026-11-11) | **oui** |
| intervalles, échelle de verdicts | `fusion_headroom_diagnostics.paired_gap` avec bootstrap par grappes, et l'échelle à six barreaux | **oui** |

**Le module nouveau est donc mince** : il compose ces pièces et n'ajoute qu'une
chose — une comparaison **entre sélecteurs** à pool, plis et estimateur
identiques.

### 1.1 Environnement : XGBoost présent, `shap` absent, et c'est sans
conséquence

`xgboost 3.4.1` est installé ; `shap` ne l'est pas. **Il n'a pas à l'être** :
XGBoost calcule lui-même le **TreeSHAP exact** via `pred_contribs=True`.
Ajouter une dépendance pour obtenir exactement le nombre que la dépendance
existante produit déjà serait un mauvais échange. `scikit-learn 1.6.1` fournit
l'elastic net.

### 1.2 Données publiques — ce qui est réellement sur le disque

| Cohorte | État | Ce qui est utilisable |
|---|---|---|
| **MaveDB** | **présente**, 3 gènes | NUDT15, CYP2C9, BRAF : abondance **et** activité sur la même librairie de variants. Pool de candidats : abondance (1) + substitution (4) + mécanisme (10) + structure (6) |
| **HLC** | **présente**, accès restreint | `curatedExpressionLiver.zip`, `curatedPhenotypesynapse.zip` : 10 activités CYP mesurées *in vitro*, ~390 donneurs, ARN + âge/sexe |
| **HLC génotype** | **présent mais coûteux** | `curatedGenotype.zip` ne porte que **deux SNP vérifiés à la source** (rs776746 CYP3A5, rs1799853 CYP2C9) et ne recoupe que **171–172 donneurs** sur ~390 |
| **EPA** | **absente** | `epa_eight_cyp_sensitivity.py` et `epa_liver_cyp_validation.py` existent comme lecteurs, mais **aucune archive EPA n'est sur cette machine.** Un lecteur sans données n'est pas une cohorte |

**Trois conséquences, déclarées avant l'essai.**

1. **Le génotype HLC entre comme strate, pas comme couche.** L'inclure
   amputerait la cohorte de plus de la moitié de ses donneurs pour deux
   variables. Un gain obtenu sur 171 donneurs ne se compare pas à un gain
   obtenu sur 390, et mélanger les deux produirait une comparaison illisible.
2. **EPA est `NOT_EVALUABLE` aujourd'hui**, pour absence de données locales, et
   l'artefact le dit explicitement au lieu de l'omettre.
3. **Aucune ligne par donneur n'entre dans le dépôt.** Les conditions de la
   source HLC l'interdisent ; l'artefact ne contient que des agrégats.

---

## 2. Conception — la seule règle qui rend la réponse interprétable

> **Tout choix supervisé a lieu à l'intérieur des plis d'entraînement, y
> compris le choix du nombre de variables à retenir.**

Un panel choisi en regardant la performance hors échantillon, puis évalué sur
cette même performance, est évalué sur sa propre sélection. Le module rapporte
donc **deux nombres distincts et ne les confond jamais** :

| | Ce que c'est | Statut |
|---|---|---|
| `curve` | la perte hors échantillon en fonction de *k*, par sélecteur | **descriptif.** Son meilleur point est choisi avec le recul sur toute la grille, donc **biaisé optimiste** — jamais le résultat principal |
| `mean_loss` | la perte hors échantillon quand *k* **et** le panel sont choisis par un découpage interne de chaque pli d'entraînement | **le résultat principal.** Il estime *la méthode qui sélectionne*, pas un panel choisi après coup |

C'est le critère §8 que vous demandiez depuis plusieurs tours et que rien
n'implémentait.

### 2.1 Ce qu'une valeur SHAP n'est pas

Une valeur SHAP est une **attribution dans un modèle ajusté**, pas une preuve
qu'une variable porte de l'information. Deux sondes colinéaires se partagent un
signal et paraissent toutes deux faibles ; une variable sur laquelle le modèle
s'est appuyé peut être inutile hors échantillon. Donc **SHAP ne fait que
proposer un ordre**, et toute affirmation de valeur est tranchée par une
**ablation hors échantillon** — réajuster sans le bloc et mesurer la différence
appariée. C'est exactement ce que demandait la note de commande.

**SHAP est calculé sur des lignes que le booster n'a pas vues** : un découpage
interne du pli d'entraînement, moyenné. SHAP calculé en échantillon classe
comme importante une variable que le modèle a surajustée, c'est-à-dire
précisément la défaillance que ce harnais existe pour éviter.

### 2.2 Les sept sélecteurs

Tous reçoivent le **même pool**, les **mêmes plis** et le **même estimateur**.
Seules les colonnes retenues diffèrent, donc le contraste porte sur la
*sélection* et non sur le modèle.

| | Sélecteur | Lit l'issue ? |
|---|---|---|
| `S0` | tout, aucune sélection — ce qu'on obtient gratuitement en ne sélectionnant pas | non |
| `S1` | corrélation de rang univariée avec l'issue | oui |
| `S2` | la même, **résidualisée sur la base** | oui |
| `S3` | elastic net, pénalité choisie dans le pli | oui |
| `S4` | **TreeSHAP exact**, calculé hors échantillon dans le pli | oui |
| `S5` | ordre issu des connaissances seules | **non** |
| `S6` | **les connaissances proposent (×3), SHAP dispose** | oui |

**Le contraste principal est `S6` contre `S1`** : la méthode guidée par le
mécanisme contre la sélection statistique naïve. `S0` est le garde-fou le plus
utile : si aucun sélecteur ne le bat, la sélection ne sert à rien sur ces
données, et c'est une réponse valable.

### 2.3 La base décide de toute la réponse

« Quelles variables sont utiles » n'a pas de sens avant que « utiles **en plus
de quoi** » soit fixé. La base est donc épinglée dans **tous** les panels :

* **MaveDB : l'abondance protéique** — la mesure qu'un laboratoire ferait de
  toute façon, et déjà de loin le meilleur prédicteur unique sur les trois
  gènes. Un sélecteur doit battre cela.
* **HLC : âge et sexe** — ce dont on dispose sans aucun dosage.

### 2.4 Contrôles

| Contrôle | Ce qu'il attrape |
|---|---|
| **contrôle positif obligatoire** : cohorte où exactement 4 colonnes sur 60 portent le signal, **placées en dernier** | un harnais incapable de voir une structure clairsemée rendrait tout résultat nul illisible. Les colonnes utiles sont en fin de pool pour qu'un bug renvoyant les *k* premières échoue bruyamment |
| **copie redondante** d'une vraie colonne (ρ > 0,9) | un sélecteur qui dépense son budget sur un doublon |
| **issue permutée**, à côté de chaque résultat réel | un sélecteur ne peut pas s'abstenir par accident : donnez-lui du bruit et il renvoie les *k* colonnes les plus chanceuses |
| **connaissances permutées** | un ordre « mécanistique » qui marche aussi bien mélangé |
| **CYP2E1** (contrôle négatif) et **CYP3A4/midazolam** (spécificité de substrat) sur HLC | inchangés, repris des analyses antérieures |
| **bootstrap par grappes** | sur MaveDB l'unité est la **position**, pas la substitution. Sur HLC une ligne **est** un donneur, donc les deux coïncident, et c'est dit plutôt que supposé |

### 2.5 Les deux questions, séparées

* **Après les mesures** — quelles variables ajoutent réellement de
  l'information ? → les contrastes entre sélecteurs, et l'ablation du bloc.
* **Avant les mesures** — les connaissances peuvent-elles dire quelle couche et
  quels marqueurs vaudraient la peine d'être mesurés ? → `knowledge_agreement` :
  un ordre **aveugle à l'issue** anticipe-t-il ce que les données finissent par
  retenir, contre un null d'ordre aléatoire ? **C'est une question plus faible
  que la valeur prédictive**, et elle est rapportée comme telle : un accord avec
  ce que les données ont choisi ne prouve pas que l'un ou l'autre a raison.

### 2.6 Stabilité, rapportée à côté de la performance et jamais à sa place

Un sélecteur qui nomme un panel différent à chaque pli n'a pas trouvé de
mécanisme, quelle que soit sa perte hors échantillon. Le harnais rapporte le
Jaccard moyen entre plis, le nombre de colonnes retenues dans **tous** les
plis, et le nombre total de colonnes jamais retenues.

---

## 3. Ce qui n'est pas fait, et pourquoi

| | Pourquoi pas |
|---|---|
| **GNN** | interdit avant démonstration de nécessité, et rien ici ne l'indique |
| nouveau modèle Fusing | branche suspendue ; les résultats négatifs sont conservés |
| STRING / Reactome dans ce premier essai | `mechanistic_graph_selection` est prêt et sera branché **après** que ce harnais ait montré qu'un sélecteur bat `S0`. Brancher un graphe avant de savoir si la sélection sert à quelque chose ajouterait une source d'échec sans ajouter d'information |
| réseau complet HLC (~20 000 sondes) | le panel ciblé est primaire et déclaré ; p ≈ 20 000 contre n ≈ 390 est un régime différent, à traiter comme une sensibilité déclarée et non glissé dans le même tableau |
| PK clinique | le critère que vous avez fixé est l'**activité enzymatique d'abord**. Rien ici ne touche une clairance, une AUC ou une concentration |

---

## 4. Résultats — exécuté sur 11 cellules

Trois gènes MaveDB × deux découpages, et cinq activités CYP mesurées sur HLC
(~390 donneurs). Tous les intervalles sont groupés par l'unité de dépendance :
la **position** sur MaveDB, le **donneur** sur HLC — où une ligne *est* un
donneur, donc les deux bootstraps coïncident.

### 4.1 Sélectionner paie. Le sélecteur, non.

| Cellule | meilleur sélecteur vs **`S0` « tout prendre »** | **`S6` hybride vs `S1` univarié** |
|---|---|---|
| NUDT15, positions | +10,0 à **+11,0 %** `GAIN_PRESENT_MARGIN_UNDECIDED` | +0,27 % **réfuté** |
| NUDT15, **régions** | +12,6 à **+13,3 % `GAIN_REACHES_MARGIN`** | −0,48 % **réfuté** |
| CYP2C9, positions | +1,85 %, réfuté | +0,01 % réfuté |
| CYP2C9, régions | −0,60 %, réfuté | +0,19 % réfuté |
| BRAF, positions | −0,36 %, réfuté | −1,16 % réfuté |
| BRAF, régions | +0,55 %, réfuté | **−5,39 % `DETECTABLY_WORSE`** |
| HLC **CYP2D6** | +6,2 à **+6,8 %** `GAIN_PRESENT_MARGIN_UNDECIDED` | −0,04 % **réfuté** |
| HLC CYP3A4T | +3,3 %, non détectable | +1,87 % non détectable |
| HLC CYP2C9 | +0,13 %, réfuté | +1,31 % réfuté |
| HLC CYP3A4M | +3,2 %, non détectable | +1,19 % non détectable |
| HLC CYP2E1 *(contrôle négatif)* | **−5,6 à −8,0 % `DETECTABLY_WORSE`** | +4,48 % `GAIN_PRESENT_MARGIN_UNDECIDED` |

**Le contraste principal n'atteint la marge sur aucune des 11 cellules**, et la
seule cellule où il est positif est le **contrôle négatif** — où `S6` fait
mieux uniquement parce qu'il sélectionne moins agressivement et se nuit donc
moins. Sur une cellule il est **détectablement pire**.

**La corrélation de rang univariée fait aussi bien que XGBoost+SHAP, que
l'elastic net, que la sélection sur le résidu et que l'hybride guidé par le
mécanisme.** Sur HLC CYP2D6 tous les sélecteurs choisissent **k = 1** dans
quatre plis sur cinq : quand le signal est une seule variable, toutes les
méthodes la trouvent, et c'est pourquoi aucune ne se détache.

**C'est une réponse utile à la question posée : ne pas construire le sélecteur
compliqué.** Le gain est dans le fait de sélectionner — jusqu'à +13 % contre
l'option gratuite qui consiste à tout prendre — et non dans la finesse de la
méthode. Un XGBoost et ses valeurs SHAP exactes n'achètent rien ici que
`|ρ_Spearman|` ne donne déjà.

### 4.2 La seule chose détectablement mauvaise : les connaissances **seules**

| | `S5` connaissances seules vs `S0` | **`S6` vs `S5`** |
|---|---|---|
| NUDT15, positions | +10,00 % `GAIN_PRESENT` | +0,61 % réfuté |
| NUDT15, **régions** | **−3,63 %, réfuté** | **+15,90 % `GAIN_REACHES_MARGIN`** |
| CYP2C9, **régions** | **−2,80 %, réfuté** | +2,51 % `GAIN_PRESENT_BUT_BELOW_MARGIN` |

**Sous le découpage dur, un ordre mécanistique aveugle à l'issue est moins bon
que l'absence de sélection**, et confronter ce même ordre aux données
d'entraînement récupère **+15,90 %**, la marge atteinte. C'est l'énoncé le plus
fort que ce projet ait produit de *« la connaissance propose, les données
disposent »* — et il est exactement invisible sous le découpage par positions,
où `S5` paraît valoir +10 %. **Croire le découpage par régions**, une fois de
plus.

### 4.3 « Avant les mesures » : les connaissances anticipent sur HLC, pas sur MaveDB

Un ordre **aveugle à l'issue** anticipe-t-il ce que les données finissent par
retenir, contre un null d'ordre aléatoire ?

| Cohorte | cellules battant un ordre aléatoire | meilleur |
|---|---|---|
| **HLC**, transcrits par donneur | **7 sur 10** | CYP3A4T : recouvrement **0,550** contre un null de 0,233, **p = 0,0005** |
| **MaveDB**, annotations de variants | **0 sur 12** | NUDT15 par régions : recouvrement **0,000** contre un null de 0,088 |

**Le clivage est net et il est utile** : « les sondes de l'enzyme cible
d'abord » anticipe quelles sondes un sélecteur piloté par SHAP retient ; « les
colonnes mécanistiques et structurales d'abord » n'anticipe pas quelles
annotations de variant il retient. Pour la question posée — *peut-on dire à
l'avance quelle couche omique mesurer ?* — la réponse est **oui pour le choix
du transcrit chez des donneurs, non pour le choix de l'annotation d'un
variant**.

**Une réserve qui compte, et que le module déclarait d'avance.** Le contrôle
négatif **CYP2E1 bat aussi l'ordre aléatoire** (k = 8, p = 0,001) — alors que
sur cet endpoint sélectionner est *détectablement nuisible*. Donc
l'anticipation porte sur **ce qu'un sélecteur retient**, pas sur **l'utilité de
ce qu'il retient**. Les deux questions sont réellement distinctes, et
`what_this_module_does_not_do` le disait avant l'exécution : *« un accord avec
ce que les données ont choisi ne prouve pas que l'un ou l'autre a raison »*.

### 4.4 Parcimonie : une à deux variables suffisent presque partout

| Cellule | *k* retenant 80 % | 90 % | *k* du meilleur point *(biaisé optimiste)* |
|---|---|---|---|
| NUDT15, régions | **2** | **2** | 2 |
| NUDT15, positions | **2** | **2** | 8 |
| BRAF, les deux | **1** | **1** | 1 |
| CYP2C9, positions | 5 | 5 | 8 |
| CYP2C9, régions | 12 | 20 | 20 |

Sur HLC les plis internes choisissent **k = 1** sur CYP2D6, CYP2E1 et presque
partout ailleurs. **Les panels utiles sont minuscules** : au-delà d'une ou deux
variables il n'y a rien à gagner sur ces données, et les cellules où *k* reste
grand (CYP2C9 par régions) sont précisément celles où aucun gain n'est
détectable — un grand panel qui n'apporte rien n'est pas un panel, c'est du
bruit conservé.

### 4.5 Stabilité : basse exactement là où il n'y a rien

Jaccard moyen entre plis, de **0,29** (CYP2C9 par régions) à **0,87** (BRAF par
positions). Les cellules les moins stables sont celles sans gain détectable. Un
sélecteur qui nomme un panel différent à chaque pli n'a pas trouvé de mécanisme,
quelle que soit sa perte — et ici les deux diagnostics concordent, ce qui est
rassurant pour l'instrument et décourageant pour la méthode.

### 4.6 Contrôles

| Contrôle | Résultat |
|---|---|
| **issue permutée**, 11 cellules | **aucune n'atteint la marge** (−3,04 % à +0,76 %). Le harnais ne lit pas ses propres plis |
| **CYP2E1**, contrôle négatif | sélectionner est **`DETECTABLY_WORSE`, −8,0 %** — exactement ce qui doit arriver là où il n'y a pas de signal clairsemé : jeter des variables coûte |
| **contrôle positif** | +63,2 % contre « tout prendre », **100 % des colonnes utiles retrouvées**, perdu sous permutation (+1,1 %, réfuté) |
| **doublon redondant** | **retenu dans 5 plis sur 5.** SHAP sur un modèle d'arbres **ne déduplique pas** |

### 4.7 Verdict sur le critère principal

> *Démontrer que la méthodologie choisit mieux les informations biologiques
> utiles pour prédire l'activité enzymatique.*

**Partiellement démontré, et pas de la manière proposée.**

* **Oui** : sélectionner bat ne pas sélectionner, jusqu'à **+13 %**, marge
  atteinte sur NUDT15 par régions et gain présent sur HLC CYP2D6.
* **Non** : un screening statistique sophistiqué ne choisit pas mieux qu'une
  corrélation de rang. **0 cellule sur 11.**
* **Oui, par la négative** : les connaissances **seules** choisissent
  détectablement *mal* sous le découpage dur, et les confronter aux données
  d'entraînement récupère **+15,90 %**. La confrontation est ce qui a de la
  valeur, pas l'ordre mécanistique.

**Ce que cela implique pour la suite.** Brancher STRING, Reactome ou un GNN
ajouterait de la sophistication à l'étape dont on vient de mesurer qu'elle
n'est pas le facteur limitant. Les deux directions que ces mesures désignent
sont ailleurs : **(1)** une étape de **redondance explicite**, puisque SHAP
garde les doublons ; **(2)** le **réseau complet HLC** (~20 000 sondes contre
un panel ciblé de quelques dizaines), car c'est le seul régime testé ici où
p ≫ n et où un sélecteur pourrait enfin se distinguer d'un autre. Sur un panel
déjà choisi par un mécanisme, il n'y a presque rien à sélectionner.

---

---

## Change log

| Version | Date | Change |
|---|---|---|
| **1.1.0** | **2026-11-11** | **Exécuté. Le contraste principal est réfuté sur 11 cellules sur 11.** Sélectionner paie — jusqu'à **+13 %** contre « tout prendre », marge atteinte sur NUDT15 par régions, **+6,4 %** sur HLC CYP2D6 — mais **aucun sélecteur ne bat les autres** : XGBoost+SHAP, elastic net, sélection sur le résidu et hybride mécanistique font tous ce que fait une corrélation de rang univariée, et sur une cellule l'hybride est **détectablement pire**. Le seul effet des connaissances est **négatif** : sous le découpage dur, `S5` (connaissances seules) est à **−3,63 %** contre « tout prendre » et `S6` le bat de **+15,90 %**, marge atteinte. « Avant les mesures » : un ordre aveugle à l'issue anticipe ce que les données retiennent sur **7 cellules HLC sur 10** (p jusqu'à 0,0005) et sur **0 des 12 cellules MaveDB** — mais le contrôle négatif CYP2E1 l'anticipe aussi, donc l'anticipation porte sur *ce qui est retenu*, pas sur *son utilité*. Parcimonie : **1 à 2 variables** retiennent 80–90 % du gain presque partout. Deux défauts corrigés en route : la recherche de pénalité d'`ElasticNetCV` n'était **pas groupée** et validait sur la position qu'elle venait d'entraîner ; et `knowledge_agreement` comparait deux espaces de noms, renvoyant un recouvrement de 0,000 sur 18 cellules — un bug qui ressemblait à un résultat. |
| 1.0.0 | 2026-11-11 | Audit des sept composants réutilisables et des trois cohortes (MaveDB présente, HLC présente avec génotype à deux SNP sur 171–172 donneurs, EPA **absente** donc `NOT_EVALUABLE`). Conception : sept sélecteurs à pool/plis/estimateur identiques, tout choix supervisé dans les plis y compris *k*, TreeSHAP exact sans la dépendance `shap`, ablation hors échantillon pour trancher la valeur, courbe performance–*k* déclarée biaisée optimiste et séparée du nombre validé, six contrôles dont un contrôle positif obligatoire à signal clairsemé et doublon redondant. |
