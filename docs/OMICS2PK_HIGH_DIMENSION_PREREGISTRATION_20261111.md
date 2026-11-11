# Pré-enregistration : trois protocoles, gelés avant exécution

**Version 5.0.0 — 2026-11-11.** Statut : **PROTOCOLES 1 (§5) ET 1b (§6) EXÉCUTÉS. PROTOCOLE 2 : `NOT_EVALUABLE` (§7) — la condition déclarée en §2.2 n'est pas remplie. Protocole 3 déclaré, aucune donnée téléchargée.**

**Chronologie, parce que c'est tout ce qui donne sa valeur à ce document.**

| Commit | Ce qui était déclaré, et quand |
|---|---|
| `b1b4f23` | protocoles 1, 2 et 3 gelés — **avant** toute exécution, et sans aucune donnée COMT ou ARAF sur le disque |
| `6e3549f` | §4.1 : pools appariés, test de non-infériorité, quatre contraintes COMT. Le protocole 1 avait été lancé ; il a été **arrêté et sa sortie partielle supprimée sans être lue** |
| `fcd2a2a` | §1b.6 : diagnostic de variance — **avant** l'exécution de 1b, et sans relancer le protocole 1 |
| `83aa8c2` | §1b.7 : restriction de 1b aux trois endpoints prédéfinis — **avant** l'exécution de 1b |
| cette version | §5 et §6 : les résultats des protocoles 1 et 1b. **Rien au-dessus de la §5 n'a été modifié après lecture d'un résultat** |

Les protocoles **2 et 3 n'ont encore aucune donnée téléchargée**, et leurs
sections restent donc intégralement antérieures à tout résultat. Une version
ultérieure ne peut qu'**ajouter** des résultats et consigner les écarts.

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
* **Le pool mécanistique ne bat pas le 5ᵉ centile des 40 pools appariés** →
  **l'avantage des connaissances n'est pas établi dans cette expérience.**
  Formulation corrigée en §1b.5 : ce n'est pas « les connaissances sont
  inutiles », et avec 40 pools transductifs ce protocole ne peut de toute
  façon pas trancher confirmatoirement — le protocole 1b est là pour cela.
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

## 1bis. Protocole 1b — la réplication confirmatoire du contrôle de pool

> **AJOUT POSTÉRIEUR, 2026-11-11 — le mot « confirmatoire » de ce titre et des
> §1b.1 à 1b.2 est RETIRÉ.** Le texte déclaré n'est pas réécrit : une
> pré-enregistration qu'on corrige après coup ne vaut rien, et ce qui a été
> déclaré doit rester lisible tel quel. Mais la qualification était fausse, et
> la raison est structurelle : **1b réutilise la cohorte HLC et les mêmes
> donneurs que le protocole 1, et il a été conçu après avoir vu le signal
> CYP2D6.** Il renforce ce signal contre une distribution de pools appariés ;
> il ne le confirme pas. Une confirmation exige de **nouveaux donneurs**. Lire
> la §6.1 avant toute conclusion tirée de cette section.

**Statut : DÉCLARÉ, NON EXÉCUTÉ.** Protocole **distinct**. Il **ne modifie pas**
le protocole 1, qui tourne et doit aller à son terme sans interruption. Le
protocole 1 reste **exploratoire sur la question du pool** pour les trois
raisons ci-dessous ; le protocole 1b est le confirmatoire.

### 1b.1 Pourquoi le protocole 1 ne peut pas être confirmatoire sur ce point

**Trois limites, toutes justes, toutes acceptées.**

1. **La résolution du p empirique.** Avec 40 pools, le plus petit p atteignable
   est **1/41 = 0,0244**. Un seuil de 0,05 corrigé par Bonferroni sur les trois
   endpoints prédéfinis vaut **0,0167**. **40 pools ne peuvent donc pas, même
   en principe, franchir ce seuil** — quelle que soit la force de l'effet. Le
   protocole 1 est structurellement incapable d'une affirmation confirmatoire
   ici, et ce n'est pas une question de puissance mais d'arithmétique.
2. **L'appariement est transductif.** Les pools du protocole 1 sont appariés une
   fois sur **tous** les donneurs retenus. L'appariement n'utilise aucune
   activité, donc il ne fuite pas d'issue — mais il utilise les valeurs
   d'expression de donneurs sur lesquels le bras va être évalué. C'est un
   design **transductif**, et une validation strictement **inductive** exige
   que l'appariement n'emploie que les lignes d'entraînement.
3. **Ne pas battre le 5ᵉ centile ne démontre pas que les connaissances sont
   inutiles.** Cela démontre que *leur avantage n'est pas établi dans cette
   expérience*. La §1.7 et la §4.1.1 du protocole 1 étaient trop fortes sur ce
   point et sont corrigées ici en §1b.5 ; le **seuil reste gelé**, c'est son
   **interprétation** qui est nuancée.

### 1b.2 Ce qui change, et rien d'autre

| | Protocole 1 (en cours) | **Protocole 1b** |
|---|---|---|
| nombre de pools | 40 | **1 000** → p minimal **0,000999**, qui franchit 0,0167 |
| appariement | une fois, sur tous les donneurs — **transductif** | **dans chaque pli, sur les lignes d'entraînement seules — inductif** |
| contrôle secondaire | aucun | **pools à redondance appariée** (§1b.3) |
| tout le reste | — | **identique** : mêmes endpoints, mêmes plis, même marge, même graine, même estimateur, mêmes filtres |

Le sélecteur `S9_matched_pool_then_univariate` tire son pool **à partir des
lignes qui lui sont passées**, qui sont les lignes d'entraînement. Un test le
vérifie par espion : aucun tirage ne voit la cohorte entière.

### 1b.3 Le contrôle secondaire : le choix des gènes, ou leur covariance ?

Votre distinction, et elle sépare deux choses que le contrôle principal
confond :

* **l'intérêt du choix des gènes et des voies** — ce que mesure le contrôle
  principal, à niveau et dispersion appariés, redondance libre ;
* **l'intérêt de leur structure de corrélation** — un pool aléatoire de même
  taille **et de même cohérence interne**. S'il égale le pool mécanistique,
  l'avantage vient de la covariance et non de l'identité des gènes.

**L'échantillonnage par rejet ne peut pas construire ce pool, et c'est mesuré,
pas supposé.** Un pool mécanistique atteignant une corrélation interne moyenne
de 0,88 est co-régulé ; 24 tirages appariés sur le niveau, dans un réseau dont
la paire de sondes typique est décorrélée, ont atteint **0,05**. Appeler cela
un contrôle à redondance appariée aurait été faux.

Le pool est donc **construit** jusqu'à la cible : une sonde de départ
aléatoire, puis des membres ajoutés dont la corrélation moyenne avec le pool en
croissance est la plus proche de la cible, chaque étape choisissant parmi un
sous-ensemble aléatoire de **1 500** candidats pour que le tirage reste
aléatoire et non déterministe. Vérifié sur un réseau synthétique à 30 modules
de co-expression : cible **0,557**, obtenu **0,554**.

**Le coût est déclaré : on ne peut pas apparier taille, niveau, dispersion et
redondance simultanément.** Ce contrôle apparie la taille et la redondance ; le
niveau et la dispersion dérivent, et les quatre chiffres sont rapportés des deux
côtés.

**Et sa faisabilité dépend du réseau réel.** Si HLC ne contient aucun module de
cohérence comparable au panel mécanistique, la redondance atteinte restera loin
de la cible. **Ce serait un fait sur le réseau, pas un bug** — et il serait
rapporté comme tel, avec les deux nombres, plutôt que présenté comme un
appariement réussi.

### 1b.6 Diagnostic secondaire ajouté : le profil de variance du pool

**Amendement, déclaré avant exécution de 1b.** Ajouté comme **diagnostic
secondaire uniquement** : il ne modifie **ni** les contrastes principaux **ni**
les seuils gelés, et le protocole 1 n'est pas relancé pour lui.

#### 1b.6.1 D'où vient la question, et pourquoi ma formulation était fausse

Le protocole 1 a trouvé que le filtre de variance coûte **+13,98 %** sur
CYP2D6, soit précisément l'endpoint où le pool mécanistique bat le hasard
apparié. J'avais proposé que les deux faits puissent être liés — que le pool
mécanistique gagne en **contournant** un filtre que les autres bras subissent,
et que son avantage soit donc un artefact.

**Cette lecture est fausse pour la comparaison principale, et la raison est de
construction.** Les 40 pools appariés le sont sur la **log-variance**. Le pool
mécanistique et ses contreparties ont donc, par construction, la même
distribution de variance : il n'y a entre eux **aucune différence de variance à
expliquer**. Vérifié sur un pool synthétique délibérément à faible variance :
delta de percentile médian **−0,017**, delta de part exclue **0,000**.

**Ce que la distinction laisse debout :**

| Comparaison | Statut |
|---|---|
| `S8` contre `S4` sur CYP2D6 | **non interprétable** — `S4` subit un filtre qui détruit une partie du signal (déjà déclaré §1.7) |
| pool mécanistique contre les 40 pools appariés | **valable** — même procédure, sans filtre des deux côtés, et pools déjà appariés sur la variance |

Donc **l'enrichissement en sondes de faible variance ne peut pas expliquer le
+14,35 % comme un artefact du filtre.** Il peut en revanche révéler une
propriété plus intéressante : *les connaissances pharmacologiques préservent des
biomarqueurs informatifs qu'un filtrage non supervisé aurait éliminés.*

#### 1b.6.2 Les trois mesures, sans jamais lire l'activité

1. **Le percentile de variance** des sondes du pool mécanistique dans le
   transcriptome complet — médiane et moyenne.
2. **La proportion de ces sondes exclues** par le filtre des 2 000 plus
   variables.
3. **Les mêmes statistiques pour les 1 000 pools appariés**, avec les écarts et
   la distribution de ces écarts (intervalles à 2,5 % et 97,5 %).

#### 1b.6.3 Et la moitié prédictive, parce que des distributions ne suffisent pas

Une distribution de variance ne démontre pas qu'une sonde porte de
l'information. Le pool mécanistique est donc **scindé par le filtre** en deux
pools, qui sont ensuite évalués **comme des pools à part entière** :

* les sondes du pool que le filtre **conserve** ;
* les sondes du pool que le filtre **écarte**.

**Si la moitié écartée prédit mieux que la moitié conservée, alors le filtrage
non supervisé retire des biomarqueurs que les connaissances avaient retenus** —
et c'est une affirmation sur la prédiction, pas sur des distributions.

**La coupure est transductive et c'est déclaré** : elle emploie la variance de
tous les donneurs retenus, non celle de chaque pli. Elle ne lit aucune
activité, mais cela en fait un **diagnostic** et jamais un bras confirmatoire ;
les deux bras qui en dérivent héritent de ce statut.

#### 1b.6.4 Ce que ce diagnostic ne permettra pas de conclure

| Résultat | Lecture admissible |
|---|---|
| pool mécanistique nettement moins variable que ses appariés | **impossible par construction.** Si cela apparaît, c'est un défaut d'appariement et non un résultat — à traiter comme un bug |
| pool mécanistique beaucoup plus exclu par le filtre que la moyenne du transcriptome | explique son avantage sur les **bras filtrés**, et rien de plus |
| la moitié écartée par le filtre prédit mieux que la moitié conservée | **les connaissances préservent des biomarqueurs qu'un filtre non supervisé élimine.** Le seul énoncé intéressant du lot, et il reste exploratoire : coupure transductive, deux pools très petits, un seul réseau |
| la moitié conservée prédit aussi bien ou mieux | l'avantage du pool mécanistique ne tient pas à des sondes de faible variance, et la §1b.6.1 est confirmée dans les deux sens |

### 1b.7 Restriction de la famille, déclarée avant exécution

**Amendement, avant toute exécution de 1b.** Il correspond à la famille de
tests **déjà définie** par le seuil Bonferroni de 0,0167 : aucun changement
d'hypothèse, aucun changement de seuil.

| Endpoint | Dans 1b ? |
|---|---|
| CYP3A4T, **CYP2D6**, CYP2C9 | **oui** : 1 000 pools inductifs, diagnostic de variance, scission prédictive |
| CYP2E1 *(contrôle négatif)*, CYP3A4M *(spécificité)* | **non** : leurs résultats du protocole 1 sont conservés comme contrôles. Mille pools pour des endpoints qui ne participent pas au contraste principal multiplieraient le coût sans toucher à l'hypothèse |

**Le contrôle secondaire à redondance appariée porte sur CYP2D6 seulement**, à
**200** pools. Mesuré avant de déclarer : un pool construit coûte **2,46 s** à
cette largeur, donc 200 × 5 plis × 3 endpoints vaudrait deux heures pour un
contrôle secondaire. Le plancher de p devient 1/201 = **0,00498**, qui franchit
toujours 0,0167. CYP2D6 est la **priorité absolue** de 1b.

**La scission prédictive reste EXPLORATOIRE**, et son statut est inscrit dans
l'artefact lui-même : la coupure est transductive, et comparer les deux moitiés
ne démontre pas à elle seule que les connaissances sont responsables du signal.

**Interprétabilité biologique.** Pour chaque endpoint, l'artefact nomme les
gènes auxquels appartiennent les sondes du pool, et ceux des deux moitiés de la
scission. Un avantage de +14,35 % sans gène interprétable serait un résultat
différent d'un avantage porté par CYP2D6 lui-même.

### 1b.4 Le contraste principal du protocole 1b

> La perte hors échantillon du pool mécanistique contre la distribution des
> **1 000** pools appariés **tirés dans les plis**, sur les trois endpoints
> prédéfinis, avec un p empirique corrigé par Bonferroni à **0,0167**.

Secondaire déclaré : la même chose contre les pools à **redondance appariée**.

### 1b.5 Lecture, et ce qu'elle ne permet pas de dire

| Résultat | Ce qu'il autorise |
|---|---|
| pool mécanistique sous le 5ᵉ centile **et** p < 0,0167 | **les connaissances constituent un pool meilleur que des pools comparables.** C'est le résultat déterminant du projet |
| au-dessus du 5ᵉ centile | **l'avantage des connaissances n'est pas établi dans cette expérience.** Ce n'est **pas** « les connaissances sont inutiles » : un pool aléatoire apparié sur le niveau et la dispersion dans un réseau hépatique contient lui-même beaucoup de biologie pertinente, et le panel mécanistique pourrait être supérieur sans que ces 1 000 tirages ne le montrent |
| égalé par les pools à **redondance appariée** mais non par les pools à redondance libre | l'avantage vient de la **structure de covariance**, pas de l'identité des gènes. Conclusion différente, et plus intéressante que les deux autres |

---

## 5. Résultats du protocole 1

Exécuté le 2026-11-11 sur `reports/data/supervised_screening_fullarray.json`.
Cinq endpoints, **30 596 candidats** (2 clinique + 30 594 sondes complètes),
**10 044 sondes écartées pour incomplétude**, pool mécanistique réduit de 31 à
**25** par ce même filtre — déclaré, et rapporté parce qu'un pool de
connaissances qui rétrécit en silence invaliderait le bras qui en dépend.

### 5.1 Le contrôle principal : 1 endpoint sur 5, et au plancher arithmétique

Perte du pool mécanistique contre 40 pools appariés sur niveau et dispersion.

| Endpoint | pool méca. | moyenne des pools appariés | gain | p empirique | bat le 5ᵉ centile |
|---|---|---|---|---|---|
| CYP3A4T | 3245,4 | 3306,5 | +1,85 % | 0,390 | non |
| **CYP2D6** | **92,79** | **108,34** | **+14,35 %** | **0,0244** | **oui** |
| CYP2C9 | 249,1 | 241,4 | **−3,20 %** | 0,951 | non |
| CYP2E1 *(contrôle négatif)* | 1378,7 | 1359,6 | −1,41 % | 0,951 | non |
| CYP3A4M | 727,1 | 769,7 | +5,53 % | 0,171 | non |

**Lecture pré-enregistrée, §1b.5.** Un seul endpoint bat le hasard apparié, et
son p vaut **exactement 1/41 = 0,0244** — le plancher : **aucun** des 40 pools
aléatoires n'a fait aussi bien. Ce plancher **ne peut pas** franchir le seuil
Bonferroni de 0,0167 sur trois endpoints prédéfinis. Donc :

> **L'avantage des connaissances n'est pas établi par cette expérience.** Ce
> n'est pas « les connaissances sont inutiles » : sur CYP2D6 l'effet est grand
> (+14,35 %, 0/40) et c'est précisément ce que le protocole 1b, avec 1 000
> pools et un appariement inductif, est fait pour trancher.

**Et un résultat négatif net à ne pas enterrer :** sur **CYP2C9**, le pool
mécanistique est **pire que la moyenne** des pools appariés (−3,20 %, p = 0,951).
Pour cette activité, les 25 sondes mécanistiques ne valent pas 25 sondes
aléatoires de même niveau et de même dispersion.

### 5.2 Une condition de réfutation s'est déclenchée : le filtre de variance

| Coût du filtre, `S1` | verdict |
|---|---|
| **CYP2D6** | **+13,98 % `GAIN_REACHES_MARGIN`** — le réseau entier bat le réseau filtré |
| CYP3A4T | −0,73 %, réfuté sous la marge |
| CYP2C9 | +3,43 %, non détectable |
| CYP2E1 | −1,20 %, réfuté |
| CYP3A4M | −0,28 %, réfuté |

Sur CYP2D6, `S1` perd **107,1** sur les 2 000 sondes les plus variables contre
**92,2** sur le réseau entier, et il en va de même de `S2` (+14,22 %) et de
`S8` (+13,94 %). **Les sondes informatives pour l'activité CYP2D6 ne sont pas
parmi les 2 000 plus variables.**

La §1.7 avait déclaré la conséquence : **sur CYP2D6, les bras qui exigent le
filtre — `S3`, `S4`, `S7` — ne sont pas comparables aux bras entiers et ne
seront pas interprétés.** Sur les quatre autres endpoints le coût est sous la
marge et ils le sont.

**La leçon méthodologique mérite d'être retenue : un filtre de variance est
aveugle à l'issue mais pas aveugle à l'information.** Il ne fuite rien et il
peut tout de même détruire le signal. Je l'avais déclaré possible ; c'est
arrivé, sur l'endpoint le plus prometteur du lot.

### 5.3 La sophistication ne paie pas davantage à p = 2 000

`S4` (TreeSHAP exact) contre `S1` (Spearman), **à candidats identiques** :

| | gain | verdict |
|---|---|---|
| CYP3A4T | −0,99 % | non détectable |
| CYP2D6 | −1,87 % | réfuté sous la marge |
| CYP2C9 | +2,75 % | non détectable |
| CYP2E1 | −1,13 % | réfuté |
| CYP3A4M | +0,31 % | non détectable |

**Aucun gain atteignant la marge, sur aucun endpoint.** Votre nuance était
juste — p ≤ 33 ne prouvait rien — et à **p = 2 000** la réponse est la même.
`S4` n'a pas été évalué à p = 30 596, où il est incalculable : c'est la seule
largeur qui reste non testée pour lui, et elle le restera.

### 5.4 Non-infériorité du pipeline simple : établie sur 1 endpoint sur 3

| Endpoint | `S8` vs `S4` | `S8` vs `S1` largeur entière |
|---|---|---|
| CYP3A4T | `UNDECIDED` | `UNDECIDED` |
| **CYP2D6** | **`NON_INFERIOR`** | **`NON_INFERIOR`** |
| CYP2C9 | `UNDECIDED` (−6,70 %) | `UNDECIDED` (−6,55 %) |
| CYP2E1 | `NON_INFERIOR` | `NON_INFERIOR` |
| CYP3A4M | `UNDECIDED` | `NON_INFERIOR` |

Conclure « le pipeline simple suffit » exige `NON_INFERIOR` : **c'est acquis
sur CYP2D6 et sur le contrôle négatif, indécis ailleurs.** Trois `UNDECIDED`
sur cinq, et la §4.1.2 interdit de les lire comme une équivalence.

### 5.5 L'abstention fonctionne, et la prédiction §1.6 tient

Nombre de plis sur 5 ayant choisi **`k = 0`** sur **CYP2E1**, le contrôle
négatif où sélectionner nuisait :

| `S4` | `S7` | `S3` | `S0` | `S8` | `S1` | `S2` |
|---|---|---|---|---|---|---|
| **5** | **5** | **4** | 4 | 2 | 1 | 0 |

**Les bras à base d'arbres s'abstiennent totalement**, `S3` presque. La
prédiction « au moins 3 plis sur 5 » est vérifiée pour les arms les plus
agressifs et **échoue pour les sélecteurs univariés**, qui continuent de
prendre une variable. C'est une asymétrie que je n'avais pas prévue et elle est
rapportée comme telle : l'abstention est désormais disponible, mais tous les
sélecteurs ne s'en servent pas également.

Sous **issue permutée**, `S4` s'abstient dans 3 à 5 plis selon l'endpoint — le
comportement correct.

### 5.6 La déduplication n'apporte rien, et nuit une fois

`S7` contre `S4`, qui isole exactement cette étape : +1,31 % (non détectable),
0,00 %, −0,30 %, 0,00 %, et **−5,91 % `DETECTABLY_WORSE` sur CYP3A4M**.

**La déduplication explicite ne paie pas sur ces données.** Le doublon que le
contrôle positif synthétique faisait garder 5 fois sur 5 n'avait pas
d'équivalent coûteux ici. Je la conserve — elle est correcte et son coût est
mesuré — mais elle n'est pas une amélioration démontrée.

### 5.7 Contrôle d'issue permutée

**Aucun des 5 endpoints n'atteint la marge** : −1,53 % à +0,94 %, tous
`REFUTED_GAIN_BELOW_MARGIN`. Le harnais ne lit pas ses propres plis.

### 5.8 Ce que le protocole 1 autorise à dire

| Question | Réponse |
|---|---|
| les connaissances constituent-elles un pool meilleur que des pools comparables ? | **non établi.** 1 endpoint sur 5, au plancher arithmétique, sous le seuil corrigé. **Protocole 1b requis**, et désormais clairement justifié par l'effet CYP2D6 |
| la sophistication paie-t-elle en haute dimension ? | **non**, à p = 2 000 comme à p ≤ 33 |
| le pipeline simple suffit-il ? | **non-infériorité acquise sur CYP2D6 seulement**, indécise sur 3 endpoints sur 5 |
| l'abstention fonctionne-t-elle ? | **oui**, et les bras à arbres l'utilisent pleinement ; les univariés non |
| la déduplication aide-t-elle ? | **non**, et elle nuit une fois |

---

## 6. Résultats du protocole 1b

Exécuté le 2026-11-11, `reports/data/supervised_screening_1b.json`. Trois
endpoints prédéfinis, **1 000 pools appariés tirés dans les plis**, pool
mécanistique de **25** sondes sur **30 596** candidats.

### 6.1 Le contraste principal : CYP2D6 franchit, les deux autres non

| Endpoint | pool méca. | appariés (moyenne) | gain | aussi bons | p | Bonferroni 0,0167 |
|---|---|---|---|---|---|---|
| CYP3A4T | 3245,4 | 3291,8 | +1,41 % | 284/1000 | 0,285 | non |
| **CYP2D6** | **92,79** | **107,88** | **+13,98 %** | **0/1000** | **0,000999** | **franchi** |
| CYP2C9 | 249,1 | 240,7 | **−3,50 %** | 986/1000 | 0,986 | non |

**Aucun des 1 000 pools appariés n'égale le pool mécanistique sur CYP2D6.** Le
p est au plancher déclaré, il franchit le seuil corrigé, et l'appariement est
**inductif** — tiré dans chaque pli, sur les lignes d'entraînement seules. Le
+14,35 % transductif du protocole 1 devient **+13,98 %** : l'effet n'était pas
un artefact de l'appariement global.

> **« Confirmatoire » est retiré.** J'avais qualifié ce résultat de niveau
> confirmatoire ; c'est faux. Le protocole 1b a été **conçu après avoir observé
> le signal CYP2D6**, et il porte sur **la même cohorte HLC et les mêmes
> donneurs**. Il renforce considérablement ce signal contre une distribution de
> pools appariés — il ne le confirme pas de façon indépendante. Une
> confirmation exigerait **de nouveaux donneurs**. Le niveau de preuve correct
> est : **résultat positif robuste contre cette distribution, sur cette
> cohorte.**

**Et deux endpoints sur trois échouent, dont un dans le mauvais sens.** Sur
CYP2C9, **986 pools appariés sur 1 000 font mieux** que le pool mécanistique.
Même panel de 25 sondes, même procédure, même cohorte : il aide pour CYP2D6 et
il nuit pour CYP2C9. C'est une réussite **sur un endpoint**, pas une
méthodologie validée.

### 6.2 Le contrôle secondaire tranche, et dans le sens le plus fort

Redondance interne : cible **0,4806**, atteinte **0,4805** — appariement
essentiellement parfait, et la faisabilité craignée en §1b.6 n'était pas un
problème sur ce réseau.

| CYP2D6, 200 pools construits | |
|---|---|
| pool mécanistique | **92,79** |
| pools à redondance appariée, moyenne | 109,17 |
| aussi bons | **0/200** |
| p | **0,00498**, franchit 0,0167 |

Selon la §1b.5 : **l'avantage ne vient pas de la structure de covariance, il
vient de l'identité des gènes choisis.** C'était la lecture la plus
intéressante des trois possibles, et c'est celle que les données donnent.

### 6.3 Le diagnostic de variance : mon hypothèse est morte, la vôtre tient

**La mienne.** Le pool mécanistique n'est pas un pool peu variable : médiane au
**98,6ᵉ percentile**, appariés au 98,6ᵉ [98,6 ; 98,7], écart **+0,01**, part
exclue 32 % contre 30 %. `pool_is_lower_variance_than_its_matches` est **faux**
sur les trois endpoints. Exactement ce que la §1b.6.1 annonçait comme
impossible par construction.

**La vôtre, et elle est démontrée prédictivement.**

| CYP2D6, scission par le filtre *(EXPLORATOIRE)* | n | perte |
|---|---|---|
| sondes que le filtre **conserve** | 17 | **109,72** |
| sondes que le filtre **écarte** | 8 | **92,11** |

La moitié écartée fait **tout** le travail — 92,11 contre 92,79 pour le pool
entier — et elle contient **CYP2D6 lui-même**, avec ABCB1, CYP1A2, CYP2E1 et
POR. La moitié conservée ne contient pas CYP2D6.

**La sonde CYP2D6 est parmi les 32 % que le filtre de variance jette, et c'est
elle qui porte le signal.** Cela explique le coût de +13,98 % du filtre mesuré
en §5.2, et c'est votre formulation vérifiée : *les connaissances
pharmacologiques préservent un biomarqueur informatif qu'un filtrage non
supervisé aurait éliminé.*

**Spécifique à CYP2D6** : sur CYP3A4T et CYP2C9 la moitié écartée est **pire**
(3487 contre 3245 ; 260 contre 249). Statut **EXPLORATOIRE** maintenu —
coupure transductive, deux pools très petits, un seul réseau.

### 6.4 Interprétabilité biologique : les gènes effectivement retenus

**Compté en PLIS, pas en sondes.** Ma première version de ce tableau donnait
des fréquences de **sondes** et les présentait comme des plis — deux sondes du
même gène dans un seul pli y comptaient pour deux. Corrigé :

| Endpoint | retenu dans **tous** les plis | nombre de **plis** contenant le gène |
|---|---|---|
| CYP3A4T | **CYP3A4** | CYP3A4 5, CYP2B6 2, CYP3A5 2, CYP2A6 1, CYP2C19 1, CYP2C8 1, CYP2C9 1 |
| **CYP2D6** | **CYP2D6** | CYP2D6 **5**, CYP2E1 **3**, CYP2C19 2, **POR 2**, CYP2C8 1, CYP3A5 1 |
| CYP2C9 | **CYP3A4** | CYP3A4 **5**, CYP2C9 **4**, CYP2A6 3, CYP2C19 3, CYP2C8 3 |

**Et le détail par pli, qui change la lecture :**

| CYP2D6, pli | *k* | gènes retenus |
|---|---|---|
| 0 | 3 | CYP2D6, CYP2E1 |
| 1 | 5 | CYP2C19, CYP2D6, CYP2E1, POR |
| 2 | 8 | CYP2C19, CYP2C8, CYP2D6, CYP2E1, CYP3A5, POR |
| **3** | **1** | **CYP2D6 seul** |
| **4** | **1** | **CYP2D6 seul** |

**Deux plis sur cinq ont choisi la sonde CYP2D6 et rien d'autre.** POR
n'apparaît que dans **2 plis sur 5**, moins que CYP2E1. La cohérence biologique
de POR n'est donc pas une preuve d'information ajoutée — ce que la §6.6 mesure.

**Deux corrections à ce que j'avais écrit.**

* **CYP2C9 n'est pas ignoré par les données.** J'avais écrit que « le gène
  mécanistiquement correct n'est pas celui que les données choisissent » :
  c'est trop fort. CYP2C9 est retenu dans **4 plis sur 5**, CYP3A4 dans 5. Les
  données prennent **les deux**, avec CYP3A4 partout. L'énoncé correct est
  plus faible : un marqueur plausiblement tissulaire est retenu encore plus
  systématiquement que le gène cible, et le pool échoue tout de même.
* **CYP3A4T retient CYP3A4 dans les cinq plis et n'a aucun avantage** sur les
  pools appariés (p = 0,285). Retenir le bon gène ne suffit pas : il faut que
  ce gène porte une information que des sondes appariées ne portent pas. Cette
  lecture-là tient.

### 6.6 Ablation EXPLORATOIRE : le réseau, ou un transcrit dans un panel ?

`scripts/run_cyp2d6_network_ablation.py`,
`reports/data/cyp2d6_network_ablation.json`. **Hors de tout protocole gelé**,
spécifiée **après** avoir vu le résultat de 1b, sur la même cohorte et les
mêmes donneurs. Elle **borne** la revendication ; elle ne la confirme pas.

| Bras | n sondes | perte |
|---|---|---|
| âge + sexe seuls | 0 | 109,51 |
| **CYP2D6 seul** | **1** | **91,76** |
| CYP2D6 + POR | 2 | 91,90 |
| pool **sans** CYP2D6 | 24 | 107,06 |
| pool mécanistique complet | 25 | 92,79 |

| Contraste | gain | verdict |
|---|---|---|
| CYP2D6 seul vs âge + sexe | **+16,21 %** | **`GAIN_REACHES_MARGIN`** |
| pool complet vs âge + sexe | +15,27 % | `GAIN_REACHES_MARGIN` |
| **pool complet vs CYP2D6 seul** | **−1,12 %** | **`REFUTED_GAIN_BELOW_MARGIN`** |
| CYP2D6 + POR vs CYP2D6 seul | −0,15 % | `REFUTED_GAIN_BELOW_MARGIN` |
| **pool sans CYP2D6 vs âge + sexe** | **+2,24 %** | **non détectable** |
| **CYP2D6 seul, non-infériorité au pool complet** | +1,11 % | **`NON_INFERIOR`** |

**Le résultat CYP2D6 est, en entier, le transcrit CYP2D6.**

1. **Le pool de 25 sondes n'ajoute rien** à la sonde unique : −1,12 %, borne
   supérieure +0,35 sous le seuil, donc une **réfutation** et non un manque de
   puissance. La sonde seule est même **non inférieure** au pool complet.
2. **POR n'ajoute rien** : −0,15 %, réfuté. Sa cohérence mécanistique était
   réelle et son apport prédictif est nul. C'était exactement l'avertissement :
   *un raisonnement expert automatisé doit pouvoir identifier un mécanisme
   au-delà du gène cible sans supposer qu'il améliore la prédiction.*
3. **Le réseau privé de CYP2D6 ne porte rien de détectable** : +2,24 % sur
   24 sondes, intervalle contenant zéro. Il n'y a pas de signal de réseau
   indépendant du transcrit cible.

**Ce que la §6.1 devient donc.** « Les connaissances pharmacologiques
constituent un pool meilleur que 1 000 pools appariés » reste **vrai et
mesuré**. Mais le mécanisme de cette vérité est : **le pool contient le
transcrit du gène étudié, et les pools aléatoires ne le contiennent pas.** La
revendication honnête n'est pas « le raisonnement en réseau paie », c'est
**« mesurer le transcrit du gène qu'on étudie paie, et il faut savoir qu'il
faut le garder »**.

**Et ce dernier point n'est pas trivial** — c'est le seul endroit où les
connaissances font un travail qu'une statistique ne fait pas. Le filtre de
variance des 2 000 sondes les plus variables **écarte la sonde CYP2D6** (§6.3).
Un screening non supervisé l'aurait perdue ; le pool mécanistique la garde.
L'apport démontré des connaissances est donc modeste et précis : **savoir quel
transcrit unique conserver quand un filtre statistique l'élimine.**

### 6.5 Ce que le protocole 1b autorise à dire

| Question | Réponse |
|---|---|
| les connaissances constituent-elles un pool meilleur que des pools comparables ? | **oui sur CYP2D6** : 0/1000 inductifs, p = 0,000999, Bonferroni franchi. **Non confirmatoire** : même cohorte, protocole conçu après le signal. **Non sur CYP3A4T, et le contraire sur CYP2C9** |
| l'avantage est-il la covariance ou l'identité des gènes ? | **l'identité** : 0/200 pools à redondance appariée |
| est-ce un réseau ou un transcrit ? | **un transcrit.** Le pool de 25 sondes n'ajoute rien à la sonde CYP2D6 seule (−1,12 %, réfuté), POR n'ajoute rien (−0,15 %), et le réseau sans CYP2D6 ne porte rien de détectable (§6.6) |
| les connaissances préservent-elles des biomarqueurs qu'un filtre non supervisé élimine ? | **oui sur CYP2D6, exploratoire** — et c'est **le seul travail** que les connaissances font ici qu'une statistique ne fait pas |
| la méthodologie est-elle validée ? | **non.** Un endpoint sur trois, un échec inversé, et l'apport réduit à la conservation d'un transcrit unique |

**Obtenu avec une corrélation de rang sur 25 sondes.** Aucun GNN, aucune
fusion, aucun elastic net, aucun SHAP — conformément à l'interdiction de
complexifier avant démonstration de nécessité, et rien dans ces chiffres ne la
démontre.

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

## 7. Protocole 2 — `NOT_EVALUABLE`

**C'est le verdict déclaré en §2.2, pas une révision du protocole.** La §2.2
disait : *« Un jeu d'activité COMT doit exister dans le même dépôt. S'il n'y en
a pas, le protocole est `NOT_EVALUABLE` et sera rapporté comme tel. »*
Vérifié le 2026-11-11, avant tout téléchargement.

### 7.1 La condition n'est pas remplie

MaveDB contient **exactement quatre** jeux COMT, interrogés par l'API publique
(`api.mavedb.org`, recherche sur `COMT` et listage de l'expérience
`urn:mavedb:00000661-a`) :

| URN | n | Mesure |
|---|---|---|
| `urn:mavedb:00000661-c-1` | 3 463 | **abondance ARN** |
| `urn:mavedb:00000661-d-1` | 3 050 | **abondance protéique** |
| `urn:mavedb:00000661-a-1` | 3 317 | abondance polysomale moyenne |
| `urn:mavedb:00000661-b-1` | 3 117 | abondance polysomale lourde |

**Les quatre sont des mesures d'abondance. Il n'existe aucune mesure d'activité
COMT.** Les deux jeux polysomaux mesurent l'occupation ribosomale, c'est-à-dire
une étape de la traduction : une couche de plus en amont de la fonction, pas un
endpoint fonctionnel.

### 7.2 Ce que je ne fais pas, et pourquoi

Il serait facile de sauver l'expérience en prenant l'abondance protéique comme
« endpoint » et l'ARN comme prédicteur. **Ce serait remplacer la question par
une autre** — prédire une omique à partir d'une autre est l'axe A, et n'a rien
à dire sur le choix d'une mesure en vue de prédire une **fonction**. La
contrainte le disait avant l'exécution et elle est respectée.

De même, l'occupation polysomale n'est pas promue en proxy d'activité : rien ne
l'établit pour COMT, et un proxy non qualifié transformerait un
`NOT_EVALUABLE` en résultat apparent.

### 7.3 Inventaire : la question a-t-elle un véhicule ailleurs ?

La politique d'acquisition exige **au moins deux mesures candidates** et **un
endpoint fonctionnel indépendant**, chez les mêmes variants. Inventorié dans
MaveDB le 2026-11-11 :

| | |
|---|---|
| gènes avec un jeu d'**activité** | 26 |
| gènes avec un jeu d'**abondance** | 68 |
| gènes avec **les deux** | **4** |
| gènes avec **≥ 2 mesures candidates + activité** | **ARAF et BRAF, et eux seuls** |

`urn:mavedb:00001314` (ARAF) et `urn:mavedb:00001315` (BRAF) portent, sur la
même librairie, l'**activité** et l'abondance sous **trois conditions** —
départ, DMSO, et **inhibition d'HSP90** pour ARAF ; départ, HSP90i et
**activation de SOS1** pour BRAF.

**Mais c'est une version plus faible de la question.** Ces trois mesures sont
la **même couche omique dans trois conditions**, non trois couches. La question
devient *« quelle condition de mesure acquérir »* et non *« quelle couche
omique acquérir »*. Elle reste une décision par variant, antérieure à
l'acquisition, avec un endpoint fonctionnel indépendant et un coût réel — une
expérience sous inhibiteur d'HSP90 coûte plus qu'une mesure au départ.

### 7.4 Ce que cela laisse ouvert, et qui doit en décider

**Aucune donnée ARAF ou BRAF supplémentaire n'est téléchargée.** Le protocole 3
est gelé sur un contraste principal **différent** — la règle structurale 3D —
et y ajouter une politique d'acquisition serait modifier un protocole gelé.

Deux voies possibles, et le choix n'est pas le mien :

1. **Déclarer un protocole 4** pour la politique d'acquisition sur ARAF ou
   BRAF, écrit et commité avant tout téléchargement, en nommant explicitement
   la question affaiblie : choisir une **condition** de mesure, pas une couche.
2. **Laisser la question sans véhicule** dans les données publiques
   inventoriées, et le consigner. C'est un résultat sur la disponibilité des
   données, pas un échec de la méthode : **aucun dépôt public inventorié ne
   mesure deux couches omiques distinctes et une fonction chez les mêmes
   variants.**

Le protocole 2 reste `NOT_EVALUABLE` dans les deux cas.

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

## 8. Résultats du protocole 3 — ARAF

Exécuté le 2026-11-11. `reports/data/variant_mechanism_benchmark.json` (gène
ARAF) et `reports/data/araf_hsp90_proxy.json`. Construit C-terminal,
**n = 6 105** variants sur **310** positions, 1 092 à des positions que 9AXM
ne résout pas. Garde-fou passé : Lys336 NZ à **2,979 Å** de l'ANP.

### 8.1 Le contraste principal tient — et ma raison de l'avoir prédit est fausse

| §3.2 prédisait | mesuré | |
|---|---|---|
| **principal** : ablation structurale sur `M4`, régions | **−4,07 % [−0,0518, +0,0101]** `REFUTED_GAIN_BELOW_MARGIN` | **tenu** |
| ρ(abondance, activité) **négatif**, comme BRAF | **+0,1563** (+0,132 hors non-sens) | **ÉCHEC** |
| activité moyenne des faux-sens **au-dessus** de la synonyme | **0,788 contre 1,011 — en dessous** | **ÉCHEC** |
| enrichissement près du site ATP **< 1** | 0,638, intervalle **[0,318, 1,144]** qui contient 1 | partiel |
| **abondance sous HSP90i, meilleur proxy** | **−0,40 % [−0,0027, −0,0010]** `DETECTABLY_WORSE` | **ÉCHEC** |

Le principal est une **réfutation propre** : borne supérieure +0,0101 sous le
seuil de 0,0246. Mais **il tient pour une raison que les données réfutent.**

J'avais prédit la réfutation **parce que** ARAF serait, comme BRAF, un criblage
de gain de fonction où l'abondance est un proxy inversé. **Les deux prémisses
sont fausses** : ARAF est à **perte de fonction** — faux-sens à 0,788 contre
1,011 pour la synonyme, c'est-à-dire **moins** actifs que le type sauvage,
l'inverse exact de BRAF (1,420 contre 0,992) — et ρ est **positif**.

**Une prédiction juste tirée d'un raisonnement faux n'est pas une confirmation
du raisonnement.** C'est le résultat le plus instructif de ce protocole.

### 8.2 L'explication HSP90 de l'échec BRAF est réfutée

C'était la seule prédiction pouvant échouer seule, et elle échoue **dans le
sens opposé** :

| | ρ avec l'activité | perte hors échantillon |
|---|---|---|
| abondance, véhicule DMSO | **+0,1553** | **0,4594** |
| abondance, inhibition d'HSP90 | +0,1478 | 0,4612 |

`hsp90i_better_than_dmso` = **−0,40 %, `DETECTABLY_WORSE`**. Et le signe le
plus net : **les plis internes ont choisi `k = 0` dans 5 plis sur 5** pour le
bras HSP90i — ils ont préféré ne rien utiliser plutôt que ce prédicteur, alors
qu'ils retenaient le DMSO dans 4 plis sur 5.

Les deux abondances corrèlent entre elles à **ρ = 0,714**, donc ce n'est pas un
problème de mesure : les deux conditions mesurent largement la même chose, et
l'inhibition d'HSP90 n'ajoute pas l'information que j'avais postulée.

**Conséquence, et elle est coûteuse.** La §3.2 déclarait : *« Si la dernière
ligne échoue alors que les autres tiennent, l'explication HSP90 de l'échec BRAF
est fausse et l'échec reste inexpliqué. »* Elle échoue, et **deux** autres
lignes échouent avec elle. **L'échec BRAF est de nouveau inexpliqué**, et
l'attribution au mécanisme de gain de fonction, écrite dans
`OMICS2PK_BRAF_PREREGISTRATION_20261111.md` §8.2, ne peut plus porter la
généralisation qu'on lui faisait porter.

### 8.3 Ce que la portée de la règle devient

| Protéine | classe | mécanisme du criblage | ablation structurale, régions |
|---|---|---|---|
| NUDT15 | hydrolase Nudix | perte de fonction | **+12,09 %**, positive |
| CYP2C9 | cytochrome P450 | perte de fonction | +3,81 %, sous la marge |
| BRAF | kinase de signalisation | **gain** de fonction | **réfutée** |
| **ARAF** | **kinase de signalisation** | **perte** de fonction | **réfutée** |

**La ligne de partage n'est pas perte contre gain de fonction — c'est enzyme
métabolique contre kinase de signalisation.** ARAF le montre : perte de
fonction comme NUDT15, et la règle échoue tout de même. L'énoncé de la §8.3 du
dossier BRAF — *« pour une enzyme dont les variants perturbent la catalyse à un
site de liaison de ligand »* — survit, mais la clause qui compte n'est pas
celle que j'avais identifiée.

**Pourquoi c'est plausible plutôt que commode :** une kinase de signalisation
est régulée par autoinhibition, dimérisation et recrutement, et son activité
mesurée est une **sortie de voie** plusieurs étapes en aval du site
catalytique. Une hydrolase ou un P450 convertit un substrat au site même qu'on
mesure. La distance 3D à la machinerie catalytique dit quelque chose du second
cas et peu du premier — et cela se vérifie maintenant sur **deux** kinases,
avec des mécanismes de criblage **opposés**.

### 8.4 Le même clivage par découpage, une quatrième fois

Position-grouped, l'ablation structurale sur ARAF est **+4,24 %
[+0,0030, +0,0332]**, borne inférieure au-dessus de zéro. Région-grouped, elle
est réfutée. C'est le quatrième protéine où le bloc structural aide sous le
découpage facile et pas sous le dur. **Croire le découpage par régions** reste
la règle, et elle a maintenant quatre confirmations.

### 8.5 Contrôles

`M5_vs_M4_permuted_mechanism` est `DETECTABLY_WORSE` sous les deux découpages
(−0,44 %, −0,49 %), donc permuter le mécanisme ne crée pas d'avantage. Le
bootstrap par lignes aurait donné `GAIN_PRESENT_MARGIN_UNDECIDED` là où le
bootstrap par positions donne « non détectable » — un facteur de largeur
**×2,77** et **×3,36**, cohérent avec les autres protéines.

### 8.6 Ce que le protocole 3 autorise à dire

| Question | Réponse |
|---|---|
| la règle structurale 3D généralise-t-elle à ARAF ? | **non.** Réfutation propre sous le découpage dur |
| l'explication HSP90 de l'échec BRAF tient-elle ? | **non.** L'abondance sous HSP90i est un proxy **détectablement pire**, et les plis internes l'ont rejetée 5 fois sur 5 |
| la ligne de partage est-elle perte contre gain de fonction ? | **non.** ARAF est à perte de fonction et échoue quand même. Le partage est **métabolique contre signalisation** |
| ma prédiction principale était-elle juste ? | **oui, et pour une raison fausse.** Les deux prémisses qui la motivaient sont réfutées |

---

## 4. Écarts et amendements

Section en ajout seul.

### 4.2 Protocole 3 : la clause « format PDB héritable » est levée

**Écart, déclaré avant tout téléchargement de score ARAF.**

La §3.3 exigeait *« la meilleure résolution parmi les entrées RCSB de l'UniProt
d'ARAF portant un analogue d'ATP dans la chaîne ARAF, **en format PDB
héritable** »*. Cette dernière clause était une limite de **mon lecteur**, pas
un critère scientifique, et elle rendrait le protocole inévaluable à cause de
l'outil plutôt que des données.

**Inventaire complet : ARAF n'a que 7 entrées RCSB.** 1WXM et 2MSE sont des
structures RMN du seul domaine RBD ; 8VSL, 8VSM et 8VSN sont des complexes
14-3-3σ avec un **phosphopeptide** ARAF, sans domaine kinase ni analogue
d'ATP ; 9ZDQ (2,51 Å) porte un inhibiteur et non un analogue. **Une seule
entrée porte un analogue d'ATP : 9AXM**, 2,42 Å, ARAF/MEK1 avec **ANP et Mg**.

**La levée de la clause ne peut donc pas avoir été dirigée par le résultat :
la règle, avec ou sans elle, désigne exactement une entrée.** 9AXM n'existe
qu'en mmCIF — `.pdb` renvoie 404, `.cif` renvoie 200.

**Ce que la levée change :** rien de l'hypothèse, rien des seuils, rien du
contraste principal. Un lecteur mmCIF est ajouté, qui emploie `auth_asym_id`
et `auth_seq_id` — la numérotation d'auteur, celle à laquelle une position
UniProt se réfère ; `label_seq_id` recommence à 1 par entité et décalerait
chaque résidu en silence.

**Et il est validé contre le lecteur PDB sur la même entrée.** Sur 6PP9, les
deux lecteurs rendent 276 résidus 447–722, des comptes de ligands identiques,
les mêmes sept identités de résidus et **la même distance de 2,875 Å**. Un
lecteur qui divergerait du précédent réécrirait les résultats BRAF sans le
dire.

### 4.3 Protocole 3 : le construit ARAF, et le garde-fou vérifié

**Écart, déclaré avant tout téléchargement.** Le protocole BRAF employait le
construit **C-terminal** pour les deux dosages. **Le construit C-terminal
d'ARAF n'a pas d'abondance au départ** — seulement DMSO et HSP90i. Le
construit **N-terminal** est le seul à porter ensemble l'abondance au départ,
l'activité au départ et l'abondance sous HSP90i, c'est-à-dire les trois
mesures dont la §3.1 a besoin.

| Rôle | URN | n déclaré |
|---|---|---|
| **abondance, départ** | `urn:mavedb:00001314-e-1` | 5 882 |
| **activité, départ** | `urn:mavedb:00001314-f-1` | 6 076 |
| **abondance, inhibition d'HSP90** | `urn:mavedb:00001314-d-1` | 5 839 |

**Le garde-fou de numérotation, vérifié contre les coordonnées avant d'être
déclaré**, comme la §3.3 l'exigeait :

| Contrôle | Requis | Mesuré |
|---|---|---|
| 336 | LYS, la lysine catalytique (UniProt) | **LYS** |
| 354 | GLU, le glutamate de l'hélice αC | **GLU** |
| 382 | THR, le résidu portier | **THR** |
| 429 | ASP, accepteur de proton (UniProt) | **ASP** |
| 447, 448, 449 | ASP, PHE, GLY — le motif DFG | **ASP, PHE, GLY** |
| Lys336 NZ au plus proche atome d'ANP | ≤ 4,0 Å | **2,979 Å** |

Le DFG a été **trouvé en balayant le segment d'activation**, non supposé par
homologie. Et comme pour BRAF, **le garde-fou n'est pas le pont salin
Lys336–Glu354** : il est **rompu à 10,49 Å**, 9AXM étant lui aussi en
conformation αC-out — comme 6PP9 à 12,69 Å. Les deux structures sont donc
dans le même état inactif, ce qui rend la comparaison ARAF/BRAF **plus**
parallèle et non moins.

**MEK1 est la chaîne A et porte son propre ANP et son propre Mg**, exactement
comme dans 6PP9 : la restriction de chaîne est encore une fois décisive, pas
cosmétique.

### 4.4 Correction du construit ARAF, et ce que je savais en la faisant

**C'est une erreur de ma part, et une erreur de séquencement.** La §4.3 a
choisi le construit N-terminal parce qu'il portait les **trois mesures**
nécessaires. **Je n'ai pas vérifié les positions qu'il assaye.**

| Construit | positions distinctes | dont dans 302–574, la région que 9AXM résout |
|---|---|---|
| **N-terminal** (`e-1`, `f-1`, `d-1`) | 300 | **1 — soit 0 %** |
| **C-terminal** (`a-1`, `b-1`, `c-1`) | 313 | **273 — soit 87 %** |

**Le contraste principal du protocole 3 est donc `NOT_EVALUABLE` sur le
construit N-terminal** : les 5 188 variants appariés sont **tous** à des
positions non résolues, le bloc structural est constant, et l'ablation vaut
`0,0000` **par construction et non par mesure**. Un contraste arithmétiquement
forcé à zéro n'est pas un essai sur l'hypothèse ; c'est un instrument
inapplicable.

**Divulgation complète, parce qu'elle est nécessaire pour juger.** Avant de
décider de corriger, j'avais vu les chiffres du run N-terminal :
ρ(abondance, activité) = **+0,011**, couverture structurale **0 %**, et les
contrastes M3/M4 de ce construit. **Ce qui a motivé la correction est
exclusivement aveugle à l'issue** — quelles positions un construit assaye, et
quels résidus une structure résout. **J'aurais pu, et dû, le vérifier avant de
télécharger quoi que ce soit.** Les chiffres du run N-terminal sont conservés
au dossier plutôt que jetés, et ils ne sont pas lus comme la réponse du
protocole.

**Ce que la correction change :**

| Rôle | URN corrigée | n |
|---|---|---|
| **activité, départ** | `urn:mavedb:00001314-c-1` | 6 855 |
| **abondance, véhicule DMSO** | `urn:mavedb:00001314-a-1` | 6 857 |
| **abondance, inhibition d'HSP90** | `urn:mavedb:00001314-b-1` | 6 851 |

Le construit C-terminal n'a pas d'abondance « au départ » mais une abondance
**sous véhicule DMSO**, qui est le **contrôle véhicule** de l'expérience
HSP90i. Pour la prédiction de la §3.2 c'est un comparateur **meilleur** qu'un
« départ » sans rapport, et c'est le même construit que BRAF employait.

**Aucun seuil, aucune hypothèse et aucun contraste principal ne change.**

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
| **5.0.0** | **2026-11-11** | **Protocole 3 exécuté (§8). Le contraste principal tient, et ma raison de l'avoir prédit est réfutée.** L'ablation structurale sur `M4`, régions, donne **−4,07 % [−0,0518, +0,0101]** `REFUTED_GAIN_BELOW_MARGIN` — réfutation propre, comme prédit. Mais **trois prédictions sur cinq échouent** : ρ(abondance, activité) est **+0,156** et non négatif ; les faux-sens sont à **0,788 contre 1,011** pour la synonyme, donc ARAF est à **perte** de fonction et non à gain, l'inverse de BRAF ; et surtout **l'abondance sous inhibition d'HSP90 est un proxy `DETECTABLY_WORSE`** (−0,40 %), les plis internes l'ayant rejetée dans **5 plis sur 5**. **L'explication HSP90 de l'échec BRAF est donc réfutée et cet échec redevient inexpliqué.** Conséquence sur la portée : la ligne de partage **n'est pas** perte contre gain de fonction, puisque ARAF est à perte de fonction et échoue — c'est **enzyme métabolique contre kinase de signalisation**, vérifié sur deux kinases aux mécanismes opposés. Quatrième protéine où le bloc structural aide sous le découpage par positions (**+4,24 %**) et est réfuté sous celui par régions. |
| 4.1.0 | 2026-11-11 | **Protocole 3 préparé et ses deux écarts déclarés avant tout téléchargement de score ARAF (§4.2, §4.3).** La clause « format PDB héritable » est **levée** : c'était une limite de mon lecteur, pas un critère scientifique, et ARAF n'a que **7 entrées RCSB** dont **une seule porte un analogue d'ATP** — 9AXM, 2,42 Å, qui n'existe qu'en mmCIF. La règle désigne donc exactement une entrée avec ou sans la clause, et la levée ne peut pas avoir été dirigée par un résultat. Un lecteur mmCIF est ajouté, employant `auth_seq_id` et non `label_seq_id`, et **validé contre le lecteur PDB sur 6PP9** : 276 résidus 447–722, ligands identiques, sept identités identiques, **2,875 Å identique**. Second écart : le construit **N-terminal** est employé là où BRAF employait le C-terminal, parce que le C-terminal d'ARAF **n'a pas d'abondance au départ** ; les trois URN sont nommées. Le **garde-fou est vérifié contre les coordonnées avant déclaration** : sept identités (K336, E354, T382, D429, DFG 447-449, DFG **trouvé par balayage** et non supposé) et **Lys336 NZ à 2,979 Å de l'ANP**. Ce n'est **pas** le pont salin K336–E354, **rompu à 10,49 Å** : 9AXM est αC-out comme 6PP9 à 12,69 Å, ce qui rend la comparaison plus parallèle. |
| 4.0.0 | 2026-11-11 | **Protocole 2 : `NOT_EVALUABLE` (§7), par le verdict déclaré en §2.2 et non par révision.** MaveDB contient **exactement quatre jeux COMT et les quatre sont des mesures d'abondance** — ARN, protéine, et deux occupations polysomales ; **aucune activité**. Ni l'abondance protéique promue en endpoint, ni l'occupation polysomale promue en proxy d'activité : l'une remplacerait la question par l'axe A, l'autre transformerait un `NOT_EVALUABLE` en résultat apparent. **Inventaire** : sur 26 gènes avec activité et 68 avec abondance, **4 ont les deux**, et **ARAF et BRAF seuls** ont ≥ 2 mesures candidates plus une activité — mais ce sont **trois conditions d'une même couche**, pas trois couches, donc une version affaiblie de la question. Aucune donnée supplémentaire téléchargée, et le choix entre déclarer un protocole 4 ou consigner l'absence de véhicule n'est pas pris unilatéralement. |
| 3.0.0 | 2026-11-11 | **Protocole 1b exécuté (§6), puis borné par une ablation exploratoire (§6.6).** Premier résultat **positif robuste** du projet sur sa thèse centrale — **et non confirmatoire** : même cohorte, protocole conçu après le signal. Il porte sur un seul endpoint, et l ablation montre qu il se réduit à un transcrit. Sur **CYP2D6**, le pool mécanistique bat **0 des 1 000 pools appariés tirés dans les plis** : +13,98 %, p = **0,000999** au plancher, **Bonferroni 0,0167 franchi**, appariement **inductif**. Le +14,35 % transductif du protocole 1 n'était donc pas un artefact. **CYP3A4T échoue** (p = 0,285) et **CYP2C9 échoue à l'envers** : 986 pools appariés sur 1 000 font mieux. Le **contrôle à redondance appariée** (CYP2D6, 200 pools construits, redondance 0,4806 → 0,4805) donne **0/200**, p = 0,00498 : l'avantage vient de **l'identité des gènes**, pas de leur covariance. Le **diagnostic de variance** tue mon hypothèse comme annoncé (98,6ᵉ percentile des deux côtés, écart +0,01) et la **scission prédictive** confirme la vôtre : la moitié que le filtre écarte porte **tout** le signal (92,11 contre 109,72) et contient la sonde **CYP2D6**. **Interprétabilité** : CYP2D6 retenu dans les cinq plis, avec **POR**, donneur d'électrons obligatoire des CYP microsomaux ; mais sur CYP2C9 le gène retenu partout est **CYP3A4 et non CYP2C9**, ce qui explique l'échec. Une question déclarée non résolue : le pool de 25 sondes bat-il la sonde de l'enzyme cible seule ? À pré-enregistrer, pas à ajouter après coup. |
| 2.2.0 | 2026-11-11 | **Restriction de 1b déclarée avant exécution (§1b.7).** Les 1 000 pools inductifs portent sur les **trois endpoints predefinis** — la famille que le seuil Bonferroni de 0,0167 suppose déjà — et CYP2E1 et CYP3A4M conservent leurs résultats du protocole 1 comme contrôles négatif et de spécificité. Aucun changement d'hypothèse ni de seuil. Le controle secondaire a **redondance appariee** porte sur **CYP2D6 seulement, 200 pools** : un pool construit coûte 2,46 s à cette largeur et 200 × 5 × 3 vaudrait deux heures pour un contrôle secondaire ; le plancher 1/201 = 0,00498 franchit toujours 0,0167. La **scission prédictive** est inscrite **EXPLORATOIRE** dans l'artefact. Les **gènes** du pool et de chaque moitié sont nommés, pour que l'interprétabilité biologique du +14,35 % soit vérifiable et non supposée. Moteur dédié `pool_arm_loss` avec voisinages mis en cache par pli : 25 cibles contre 30 596 candidats recalculées à chaque tirage seraient 46 milliards de distances pour trois endpoints ; mises en cache, vingt calculs suffisent. |
| 2.1.0 | 2026-11-11 | **Diagnostic secondaire ajouté à 1b (§1b.6), avant exécution, sans toucher aux contrastes principaux ni aux seuils, et sans relancer le protocole 1.** Ma formulation initiale était fausse et la correction est enregistrée : les 40 pools appariés le sont **sur la log-variance**, donc il n'existe par construction aucune différence de variance entre eux et le pool mécanistique — vérifié à un delta de percentile médian de −0,017 sur un pool synthétique à faible variance. L'enrichissement en sondes peu variables **ne peut donc pas** expliquer le +14,35 % comme un artefact du filtre ; il peut seulement expliquer l'avantage sur les bras **filtrés**. Trois mesures déclarées, toutes aveugles à l'issue : percentile de variance du pool, part exclue par le filtre des 2 000, et les mêmes pour les 1 000 pools appariés avec écarts et distribution. Plus une **moitié prédictive** — le pool est scindé par le filtre et les deux moitiés sont évaluées comme des pools — parce qu'une distribution de variance ne démontre pas qu'une sonde porte de l'information. La coupure est **transductive**, donc diagnostic et jamais confirmatoire, et les bras qui en dérivent héritent de ce statut. |
| 2.0.0 | 2026-11-11 | **Protocole 1 exécuté (§5).** Contrôle principal : le pool mécanistique bat le hasard apparié sur **1 endpoint sur 5** (CYP2D6, +14,35 %, p = **0,0244** qui est le plancher 1/41 — 0 des 40 pools aussi bon), ce qui **ne peut pas** franchir Bonferroni 0,0167, donc **l'avantage des connaissances n'est pas établi** et le protocole 1b est requis. Sur **CYP2C9 le pool mécanistique est pire que la moyenne des pools appariés** (−3,20 %). **Une condition de réfutation pré-déclarée s'est déclenchée** : le filtre de variance coûte **+13,98 % `GAIN_REACHES_MARGIN`** sur CYP2D6, donc `S3`, `S4`, `S7` n'y sont pas interprétables — un filtre de variance est aveugle à l'issue mais **pas à l'information**. `S4` contre `S1` n'atteint la marge sur aucun endpoint à **p = 2 000**, confirmant à une largeur nouvelle ce que p ≤ 33 ne pouvait pas trancher. Non-infériorité de `S8` acquise sur **CYP2D6 seulement**, `UNDECIDED` sur 3 des 5. **L'abstention fonctionne** : sur CYP2E1 les bras à arbres choisissent `k = 0` dans **5 plis sur 5**, mais les univariés dans 0 à 1 — asymétrie non prévue, rapportée. **La déduplication n'apporte rien** et est `DETECTABLY_WORSE` une fois (−5,91 %, CYP3A4M). Issue permutée : aucune cellule n'atteint la marge. |
| 1.2.0 | 2026-11-11 | **Protocole 1b déclaré, distinct, sans modifier le protocole 1 qui tourne.** Le protocole 1 reste **exploratoire sur la question du pool** pour trois raisons acceptées : avec 40 pools le p minimal est **1/41 = 0,0244** et ne peut pas franchir un Bonferroni de 0,0167 sur trois endpoints, ce qui est une impossibilité arithmétique et non un manque de puissance ; son appariement est **transductif**, fait une fois sur tous les donneurs ; et « ne pas battre le 5ᵉ centile » ne démontre pas l'inutilité des connaissances mais seulement que leur avantage n'est pas établi ici — le seuil reste gelé, son interprétation est corrigée. Le protocole 1b porte **1 000 pools** (p minimal 0,000999), un appariement **inductif tiré dans chaque pli sur les lignes d'entraînement seules** (sélecteur `S9`, vérifié par espion), et un **contrôle secondaire à redondance appariée** qui sépare l'intérêt du choix des gènes de celui de leur covariance. Ce dernier ne peut pas être obtenu par rejet — mesuré : cible 0,88, atteint 0,05 sur 24 tirages — donc le pool est **construit** jusqu'à la cible, vérifié à 0,554 contre 0,557 sur un réseau à 30 modules, au prix déclaré d'une dérive du niveau et de la dispersion. |
| 1.1.0 | 2026-11-11 | **Trois ajouts, antérieurs à tout résultat** (§4.1). Le protocole 1 avait été lancé ; il a été **arrêté et sa sortie partielle supprimée sans être lue**. (1) **40 pools aléatoires appariés** sonde par sonde sur la moyenne et la log-variance, 25 plus proches voisins, redondance interne mesurée mais délibérément non appariée — parce que comparer 31 sondes mécanistiques à 30 218 faisait varier la biologie **et** la dimension, et qu'un pool mécanistique pouvait gagner en ayant seulement réduit la dimension. Nouvelle condition de réfutation : ne pas battre le 5ᵉ centile des pools appariés signifie dimension réduite, pas connaissance appliquée. (2) **`non_inferiority()`**, un test séparé appliquant la marge à **moins** le seuil, parce qu'une différence non significative n'est pas une équivalence et que « le pipeline simple suffit » exige `NON_INFERIOR`. (3) **Protocole 2** : le comparateur devient la meilleure politique **fixe choisie sur l'entraînement**, « ne rien mesurer » devient une action admissible, le **coût d'acquisition** entre dans la comparaison, et **la couche acquise ne peut jamais être l'endpoint** — imposé par assertion, même circularité que la règle MAPBE. |
| 1.0.0 | 2026-11-11 | Trois protocoles déclarés. **Protocole 1**, HLC complet : 466 × 40 638 sondes, 30 218 complètes, panel mécanistique 31, deux filtres aveugles à l'issue (complétude ; 2 000 plus variables dans le pli) dont le coût est mesuré par ablation appariée contre la largeur entière, sept bras, contraste principal **`S8` contre `S4`** déclaré singulier avec lecture d'**équivalence**, six conditions de réfutation dont une qui déclare explicitement l'absence de puissance comme un échec. **Protocole 2**, COMT : conditionnel à l'existence d'un jeu d'activité, sinon `NOT_EVALUABLE` ; contraste principal une **politique d'acquisition par variant** contre la meilleure couche fixe, avec l'oracle comme plafond, distincte des nuls antérieurs sur la pondération conditionnelle. **Protocole 3**, ARAF : réplication mot pour mot du protocole BRAF, plus la seule prédiction qui peut échouer seule — l'abondance sous inhibition d'HSP90 comme meilleur proxy. Corrections `k = 0` et déduplication consignées comme déjà faites, et donc non pré-enregistrées. |
