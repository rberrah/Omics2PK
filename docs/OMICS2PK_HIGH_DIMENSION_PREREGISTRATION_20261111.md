# Pré-enregistration : trois protocoles, gelés avant exécution

**Version 1.1.0 — 2026-11-11.** Statut : **DÉCLARÉ, AUCUN NON ENCORE EXÉCUTÉ.** Le protocole 1 a été lancé puis **arrêté et sa sortie partielle supprimée sans être lue**, pour que les ajouts de la §4.1 soient antérieurs à tout résultat. Aucun résultat d'aucun protocole n'a été examiné à ce jour.

Ce document est commité **avant** l'exécution du protocole 1 et **avant tout
téléchargement** des données des protocoles 2 et 3. Rien ici ne peut changer
après lecture d'un résultat ; une version ultérieure ne peut qu'**ajouter** une
section de résultats et consigner les écarts.

---

## 0. Pourquoi geler, et pourquoi maintenant

Le commit `cc974c9` a trouvé que la sophistication du sélecteur n'achète rien
sur 11 cellules — mais à **p ≤ 33**. Ce n'est pas une réfutation en haute
dimension, et personne ne devrait la lire comme telle. Le test qui départage
réellement les méthodes n'a pas encore eu lieu.

Ce test a une propriété dangereuse : **le résultat espéré est un nul.**
« Le pipeline mécanistique + Spearman reste aussi performant que les autres »
est une affirmation d'équivalence. Un nul n'est intéressant que s'il est une
**réfutation** — borne supérieure sous la marge — et non une absence de
puissance. Choisir le filtre, la taille de panel ou le seuil après avoir vu la
courbe transformerait cette équivalence en artefact. D'où ce document.

Les deux corrections demandées sont **déjà en place** et ne sont pas
pré-enregistrées, parce qu'elles sont des corrections et non des expériences :

| Correction | État |
|---|---|
| **`k = 0`**, l'abstention, décidée dans les plis d'entraînement | faite. `FROZEN["k_grid"] = (0, 1, 2, 3, 5, 8, 12, 20, 32, 50)`. Un panel vide prédit la moyenne d'entraînement, explicitement, plutôt que de faire échouer un solveur. Le nombre de plis ayant choisi `k = 0` est rapporté par sélecteur |
| **déduplication explicite** | faite. `drop_redundant`, glouton, respectant l'ordre, calculé sur les lignes d'entraînement, seuil gelé **0,90**. La base n'est jamais écartée. Nouveau sélecteur `S7_shap_deduplicated` = `S4` + déduplication, donc le contraste apparié `S7` contre `S4` **isole la déduplication et rien d'autre** |

---

## 1. Protocole 1 — le transcriptome HLC complet. p ≫ n pour la première fois

**La priorité que vous avez fixée.** À exécuter en premier.

### 1.1 Données et dimensions, mesurées avant de geler

| | |
|---|---|
| source | HLC, accès restreint, `curatedExpressionLiver.zip` + `curatedPhenotypesynapse.zip` |
| matrice | **466 donneurs × 40 638 sondes** |
| sondes finies chez tous les donneurs | **30 218** |
| panel mécanistique | **31 sondes** (15 gènes : 10 CYP, POR, ABCB1, SLC22A2, SLC47A1, SLC47A2) |
| donneurs retenus par endpoint | activité finie et ≥ 0, âge et sexe finis — ~360 |

**Aucune ligne par donneur n'entrera dans le dépôt.** Les conditions de la
source l'interdisent ; l'artefact ne contiendra que des agrégats et des
identifiants de sonde.

### 1.2 Les deux filtres, déclarés, et tous deux aveugles à l'issue

1. **Complétude.** Une sonde non finie chez au moins un donneur retenu est
   écartée. Le nombre de sondes écartées, et combien d'entre elles
   appartenaient au panel mécanistique, sont rapportés — un pool de
   connaissances qui rétrécit en silence invaliderait le bras qui en dépend.
2. **Variance, dans le pli.** Pour les bras qui ne peuvent pas tourner sur
   30 000 colonnes, les **2 000** candidats les plus variables sont conservés,
   variance calculée sur les **lignes d'entraînement seules**.

**Ni l'un ni l'autre ne lit une activité**, donc aucun ne peut fuiter. Le filtre
de variance existe parce qu'un elastic net et un TreeSHAP exact sont
**incalculables** à 30 000 colonnes, pas parce qu'il améliore quoi que ce soit.

**Et son coût est mesuré, pas supposé nul.** `S1`, `S2` et `S8` tournent **sans
filtre sur la largeur entière** et **avec** le filtre, sur lignes et plis
identiques ; l'ablation appariée des deux donne ce que le filtre coûte. Un
filtre de variance **peut** écarter une sonde informative de faible variance :
c'est un coût de calculabilité et il sera chiffré.

### 1.3 Les bras

| | Sélecteur | Largeur |
|---|---|---|
| `S0` | tout, aucune sélection | filtré |
| `S1` | Spearman univarié | **entière** et filtrée |
| `S2` | Spearman résidualisé sur âge et sexe | **entière** et filtrée |
| `S3` | elastic net, pénalité choisie **par plis groupés** dans le pli | filtré |
| `S4` | TreeSHAP exact, calculé hors échantillon dans le pli | filtré |
| `S7` | `S4` + suppression des quasi-doublons | filtré |
| **`S8`** | **le pool mécanistique (31) délimite, Spearman ordonne dedans** | **entière** et filtré |

`S5` et `S6` ne sont **pas** reconduits ici : l'ordre mécanistique a été mesuré
détectablement mauvais sous le découpage dur, et `S8` pose la question
différente et plus intéressante — les connaissances **délimitent** au lieu
d'**ordonner**.

### 1.4 Tout ce qui est gelé

| Paramètre | Valeur |
|---|---|
| grille de budgets | **(0, 1, 2, 3, 5, 8, 12, 20, 32, 50)**, et *k* choisi par un découpage interne de chaque pli |
| filtre de variance | **2 000** |
| taille du pool pour `S8` | **31**, moins les sondes écartées pour incomplétude |
| seuil de redondance | **0,90** |
| marge, plus petit gain intéressant | **5 %** |
| graine | 20261112 |
| plis externes / internes | 5 / 3, groupés |
| perte | erreur absolue |
| unité de rééchantillonnage | le **donneur** — et une ligne *est* un donneur ici, donc bootstrap par grappes et par lignes coïncident |
| estimateur final | ridge, identique pour tous les sélecteurs, parce que le contraste porte sur la sélection |
| XGBoost | configuration gelée de `hlc_real_xgboost_validation`, non réglée |
| taille de famille | 8 |
| endpoints | **CYP3A4T, CYP2D6, CYP2C9** prédéfinis ; **CYP2E1** contrôle négatif ; **CYP3A4M** contrôle de spécificité de substrat |

### 1.5 Le contraste principal, déclaré singulier

> **`S8_pool_then_univariate` contre `S4_shap_xgboost`**, sur les trois
> endpoints prédéfinis, groupé par donneur.

C'est exactement l'affirmation à tester : *un pool mécanistique trié par
Spearman vaut-il un screening XGBoost/SHAP quand p ≫ n ?*

**Lu en deux temps** : la supériorité par l'échelle à six barreaux, et la
**non-infériorité de `S8`** par le test de la §4.1.2. Conclure que le pipeline
simple suffit exige **`NON_INFERIOR`**, pas une absence de significativité.

**Et le principal ne suffit pas sans le contrôle de la §4.1.1** : `S8` contre
`S4` ou contre `S1` ne dit pas si la connaissance a servi, seulement si un
petit pool a servi. Les **40 pools aléatoires appariés** le disent.

Secondaires déclarés, rapportés mais non décisifs : `S8` contre `S1` (les
connaissances délimitent-elles mieux que ne cherche la largeur entière) ;
`S7` contre `S4` (ce qu'apporte la déduplication) ; chaque sélecteur contre
`S0` ; `S4` contre `S1` ; le coût du filtre de variance ; et le nombre de plis
ayant abstenu.

### 1.6 La prédiction, avec des nombres

| | Prédit |
|---|---|
| **principal** : `S8` contre `S4` | **équivalence** — borne supérieure sous la marge de 5 %, dans les deux sens |
| `S4` contre `S1` | pas de gain atteignant la marge, **mais** un gain détectable est plausible ici et ne l'était pas à p ≤ 33 |
| `S8` contre `S1` | `S8` ≥ `S1` : délimiter devrait au moins ne pas nuire |
| `S7` contre `S4` | gain positif, probablement **sous** la marge |
| **CYP2E1**, contrôle négatif | **au moins 3 plis sur 5 choisissent `k = 0`** — c'est le test direct de la correction d'abstention |
| issue permutée | aucune cellule n'atteint la marge |

### 1.7 Ce qui réfute quoi

* **`S8` contre `S4` détectablement négatif** → un pool mécanistique est pire
  qu'un screening complet en haute dimension. L'automatisation du raisonnement
  NetPK par un pool ne tient pas, et il faut le dire.
* **`S8` contre `S4` avec un intervalle contenant zéro et dépassant la marge**
  → **ni équivalence ni différence** : le protocole manque de puissance et ne
  répond pas. **Ce n'est pas un succès**, et il ne sera pas rapporté comme tel.
* **`S4` contre `S1` atteignant la marge** → la sophistication paie en haute
  dimension, et la conclusion du commit `cc974c9` est bornée à p faible, comme
  vous l'aviez pressenti.
* **CYP2E1 n'abstient dans aucun pli** → la correction `k = 0` ne fonctionne
  pas là où elle était justifiée, et c'est un défaut d'implémentation, pas un
  résultat.
* **Le filtre de variance coûte un gain atteignant la marge** → les bras
  filtrés ne sont pas comparables aux bras entiers, et seuls `S1`, `S2`, `S8`
  pourront être interprétés.
* **Le pool mécanistique ne bat pas le 5ᵉ centile des 40 pools appariés** → il
  a **réduit la dimension et non appliqué des connaissances**. C'est un
  résultat négatif sur la thèse centrale du projet et il sera rapporté comme
  tel (§4.1.1).
* **`S8` contre `S4` sort `UNDECIDED` au test de non-infériorité** → le
  protocole n'autorise pas la conclusion « le pipeline simple suffit », quelle
  que soit la valeur ponctuelle (§4.1.2).

### 1.8 Ce qui sera rapporté quoi qu'il arrive

Le contraste principal, son intervalle, son verdict, les deux contrôles, le
coût du filtre et le compte d'abstentions. **Une absence de gain est un
résultat scientifique valable.** Aucun bras ne sera ajouté, aucun filtre
déplacé, aucun endpoint échangé après lecture des résultats ; tout écart
figurera en §4 avec sa raison avant qu'un contraste ne soit cité.

---

## 2. Protocole 2 — COMT : quelle **couche** mesurer, pas quelle variable

**Statut : déclaré, aucune donnée téléchargée.**

### 2.1 Pourquoi ce jeu et pas un autre

`urn:mavedb:00000661` porte, sur **la même librairie de variants**, quatre
mesures : **abondance ARN**, **abondance protéique**, **occupation polysomale
moyenne** et **occupation polysomale lourde**. C'est le seul jeu public
inventorié où plusieurs *couches omiques* sont mesurées sur les mêmes variants.

C'est votre question §« point scientifique essentiel » prise au mot : non pas
*quelles variables prédisent* mais **quelle information est nécessaire pour
décider quelle mesure réaliser**.

### 2.2 Ce qui doit être vérifié avant toute analyse

**Un jeu d'activité COMT doit exister dans le même dépôt.** S'il n'y en a pas,
le protocole est **`NOT_EVALUABLE`** et sera rapporté comme tel : sans mesure
fonctionnelle il n'y a pas d'issue, et prédire une couche depuis une autre
serait l'axe A, pas celui-ci. **Cela sera vérifié avant de geler le reste**, et
le gel de la §2.3 est conditionnel à cette vérification.

### 2.3 Le contraste principal, si le jeu est évaluable

> **Une politique d'acquisition par variant**, choisie dans les plis
> d'entraînement, contre **la meilleure couche unique fixe**.

La politique voit, pour chaque variant, **uniquement** les connaissances
— position, substitution, mécanisme présumé, structure si disponible — et
choisit **une seule** couche à acquérir. La perte est celle de la prédiction
faite avec la couche choisie. Comparateurs : toujours l'ARN, toujours la
protéine, choix aléatoire, et l'oracle par variant comme **plafond**.

**Ce n'est pas la question déjà tranchée.** `M5` et `M6` disposaient des deux
mesures et les pondéraient ; une politique d'acquisition n'en a **qu'une** et
doit choisir avant de voir quoi que ce soit. Les nuls antérieurs sur la
pondération conditionnelle **ne règlent pas** celle-ci, et le plafond de
l'oracle dira si la question a une réponse possible avant de chercher à
l'atteindre.

---

## 3. Protocole 3 — ARAF, et l'explication de l'échec BRAF

**Statut : déclaré, aucune donnée téléchargée.**

`urn:mavedb:00001314` porte, sur la même librairie, l'abondance **au départ**,
l'abondance **sous inhibition d'HSP90**, et l'**activité**. La réfutation BRAF
a été attribuée à un mécanisme de gain de fonction où l'abondance est un proxy
**inversé** ; l'abondance sous inhibition d'HSP90 teste directement cette
explication, puisque les variants activateurs déstabilisés en dépendent.

### 3.1 Le contraste principal

> **Pour ARAF, la règle structurale 3D ajoute-t-elle à l'abondance ?** Ablation
> appariée du bloc structural sur `M4`, groupée par région contiguë,
> intervalles groupés par position — identique au protocole BRAF, mot pour mot.

### 3.2 La prédiction, et pourquoi elle est risquée

| | Prédit |
|---|---|
| ρ(abondance, activité) | **négatif**, comme BRAF (−0,290) |
| activité moyenne des faux-sens | **au-dessus** de la synonyme, donc gain de fonction |
| **principal** | **réfuté**, comme BRAF |
| enrichissement près du site ATP | **< 1** |
| l'abondance sous inhibition d'HSP90 devient-elle un meilleur proxy que l'abondance de base ? | **oui** — c'est le seul élément nouveau, et la seule prédiction qui pourrait échouer indépendamment |

**Si la dernière ligne échoue alors que les autres tiennent**, l'explication
HSP90 de l'échec BRAF est fausse et l'échec reste inexpliqué. C'est le risque
réel de ce protocole, et c'est pourquoi il vaut la peine d'être couru.

### 3.3 La structure, par la même règle que BRAF

Même règle de sélection : la meilleure résolution parmi les entrées RCSB de
l'UniProt d'ARAF portant un analogue d'ATP **dans la chaîne ARAF**, en format
PDB héritable. **Le garde-fou de numérotation sera vérifié contre les
coordonnées avant d'être déclaré** — la leçon de BRAF, où le pont salin
Lys–Glu évident est rompu à 12,69 Å dans une conformation αC-out et aurait
échoué sur une numérotation correcte.

---

## 4. Écarts et amendements

Section en ajout seul.

### 4.1 Trois ajouts, antérieurs à tout résultat

Demandés après le gel de la version 1.0.0 et **avant la lecture du moindre
résultat** : le protocole 1 avait été lancé, il a été **arrêté et sa sortie
partielle effacée sans être lue**, puis relancé après ces ajouts. C'est la
seule façon dont un ajout peut rester admissible.

#### 4.1.1 Pools aléatoires appariés — le contrôle qui isole la connaissance

**L'objection est juste et le protocole 1.0.0 ne pouvait pas y répondre.**
Comparer `S8` (31 sondes mécanistiques) à `S1` (30 218 sondes) fait varier
**deux** choses à la fois : la biologie employée pour choisir, et la dimension
fouillée. Un pool mécanistique pourrait gagner en n'ayant **réduit la
dimension**, sans avoir appliqué aucune connaissance.

**Nouveau contrôle, déclaré ici.** `N_MATCHED_POOLS = 40` pools aléatoires de
la même taille, chacun apparié **sonde par sonde** sur la moyenne et la
log-variance : chaque sonde mécanistique est remplacée par un tirage uniforme
parmi ses **25** plus proches voisins dans l'espace standardisé (moyenne,
log-variance), sans remise, en excluant le pool réel et la base. Chaque pool
est ensuite trié par la même étape Spearman. La statistique est la position de
la perte hors échantillon du pool mécanistique dans la distribution des pertes
des pools appariés, avec un p empirique.

**L'appariement est aveugle à l'issue** — moyenne, variance et corrélation ne
font intervenir aucune activité — donc les tirages sont faits une fois sur les
donneurs retenus sans rien faire fuiter dans un pli.

**La redondance interne est mesurée des deux côtés mais n'est pas appariée**, et
c'est un choix, pas un oubli : un pool mécanistique est corrélé **parce que**
ses membres sont co-régulés, ce qui fait partie de ce qui en fait un mécanisme.
Forcer un pool aléatoire à la même corrélation interne y introduirait un
morceau de la structure même que la comparaison doit tester. Les deux chiffres
sont rapportés pour qu'un lecteur voie si les pools diffèrent sur ce point.

**Nouvelle condition de réfutation, et elle est sévère :** *si le pool
mécanistique ne bat pas le 5ᵉ centile des pools appariés, il a réduit la
dimension et non appliqué des connaissances.* Ce serait un résultat négatif sur
la thèse centrale du projet, et il sera rapporté comme tel.

#### 4.1.2 La non-infériorité devient un test, pas une lecture

La §1.6 annonçait une « équivalence » lue sur l'échelle de supériorité, ce qui
était insuffisant. `non_inferiority()` est ajoutée : avec
`gap = perte(titulaire) − perte(candidat)`, le verdict est **`NON_INFERIOR`**
quand la borne **inférieure** dépasse **moins** la marge, **`INFERIOR`** quand
la borne supérieure est sous moins la marge, et **`UNDECIDED`** sinon —
explicitement **pas** une équivalence. Marge inchangée à 5 %.

Le principal du protocole 1 se lit donc en deux temps : la supériorité par
l'échelle à six barreaux, et **la non-infériorité de `S8` par ce test-ci**. La
conclusion « le pipeline simple suffit » exige **`NON_INFERIOR`**, pas
l'absence de significativité.

#### 4.1.3 Protocole 2 — quatre contraintes ajoutées

| Contrainte | Pourquoi |
|---|---|
| le comparateur principal est **la meilleure politique fixe choisie sur l'entraînement**, pas la meilleure couche connue après coup | sinon le comparateur bénéficie d'une information que la politique n'a pas, et le contraste est truqué dans le sens inverse |
| **« ne rien mesurer » est une action admissible** de la politique | une politique d'acquisition qui ne peut pas s'abstenir n'est pas une politique, c'est une obligation d'achat. Même raison que `k = 0` |
| le **coût d'acquisition** entre dans la comparaison, et non la seule perte | à perte égale, acquérir une couche de moins est strictement meilleur. Les coûts seront déclarés en unités relatives explicites avant exécution, puisqu'aucun prix réel n'est disponible publiquement |
| **la couche acquise ne peut jamais être l'endpoint** | sinon la politique « acquérir l'activité » prédit l'activité par elle-même. C'est la même circularité que la règle MAPBE du cahier des charges, et elle sera imposée par une assertion, pas par une convention |

---

---

## Change log

| Version | Date | Change |
|---|---|---|
| **1.1.0** | **2026-11-11** | **Trois ajouts, antérieurs à tout résultat** (§4.1). Le protocole 1 avait été lancé ; il a été **arrêté et sa sortie partielle supprimée sans être lue**. (1) **40 pools aléatoires appariés** sonde par sonde sur la moyenne et la log-variance, 25 plus proches voisins, redondance interne mesurée mais délibérément non appariée — parce que comparer 31 sondes mécanistiques à 30 218 faisait varier la biologie **et** la dimension, et qu'un pool mécanistique pouvait gagner en ayant seulement réduit la dimension. Nouvelle condition de réfutation : ne pas battre le 5ᵉ centile des pools appariés signifie dimension réduite, pas connaissance appliquée. (2) **`non_inferiority()`**, un test séparé appliquant la marge à **moins** le seuil, parce qu'une différence non significative n'est pas une équivalence et que « le pipeline simple suffit » exige `NON_INFERIOR`. (3) **Protocole 2** : le comparateur devient la meilleure politique **fixe choisie sur l'entraînement**, « ne rien mesurer » devient une action admissible, le **coût d'acquisition** entre dans la comparaison, et **la couche acquise ne peut jamais être l'endpoint** — imposé par assertion, même circularité que la règle MAPBE. |
| 1.0.0 | 2026-11-11 | Trois protocoles déclarés. **Protocole 1**, HLC complet : 466 × 40 638 sondes, 30 218 complètes, panel mécanistique 31, deux filtres aveugles à l'issue (complétude ; 2 000 plus variables dans le pli) dont le coût est mesuré par ablation appariée contre la largeur entière, sept bras, contraste principal **`S8` contre `S4`** déclaré singulier avec lecture d'**équivalence**, six conditions de réfutation dont une qui déclare explicitement l'absence de puissance comme un échec. **Protocole 2**, COMT : conditionnel à l'existence d'un jeu d'activité, sinon `NOT_EVALUABLE` ; contraste principal une **politique d'acquisition par variant** contre la meilleure couche fixe, avec l'oracle comme plafond, distincte des nuls antérieurs sur la pondération conditionnelle. **Protocole 3**, ARAF : réplication mot pour mot du protocole BRAF, plus la seule prédiction qui peut échouer seule — l'abondance sous inhibition d'HSP90 comme meilleur proxy. Corrections `k = 0` et déduplication consignées comme déjà faites, et donc non pré-enregistrées. |
