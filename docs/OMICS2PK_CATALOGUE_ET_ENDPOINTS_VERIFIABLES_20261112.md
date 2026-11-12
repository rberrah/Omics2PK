# Catalogue, recherche de données et endpoints vérifiables

**Version 1.0.0 — 2026-11-12.** Phases 1 et 2 du recentrage. **Aucune
expérience n'a été lancée.** Module : `src/omics2pk/input_catalogue.py`.
Script : `scripts/build_hypothesis_dataset_matrix.py`. Artefact :
`reports/data/hypothesis_dataset_matrix.json`.

---

## 0. Un conflit d'étiquettes à régler avant toute matrice

Le registre porte déjà des étiquettes `H1`–`H5`. **Ce ne sont pas vos
hypothèses.**

| étiquette | registre (ancien objectif) | nouvelle définition |
|---|---|---|
| H1 | valeur incrémentale des omiques (ARIANA) | omiques vs PGx + clinique + biologie, **pour la PK** |
| H2 | fiabilité de la sélection mécanistique | sélection automatique, **avant acquisition** |
| H3 | parcimonie et information minimale | omiques vs PGx pour les **fonctions biologiques** |
| **H4** | **« quand les omiques nuisent »** | **longitudinal, variabilité intra-individuelle** |
| **H5** | **robustesse biologique, analytique, temporelle** | **personnalisation posologique** |

H1 et H2 sont compatibles ; **H3, H4 et H5 n'ont aucune correspondance**.
Treize jeux portent une étiquette `H5` ancienne et six une `H4` ancienne : les
reprendre aurait fait revendiquer du longitudinal et de la posologie qui
n'existent nulle part. `assert_hypothesis_vocabulary` exige désormais qu'on
dise de quel vocabulaire on parle, et refuse sinon.

---

## 1. Le catalogue (§3 du recentrage)

`input_catalogue.py` porte les cinq volets sous forme exploitable : **16**
entrées cliniques, **12** biologiques, **10** descripteurs de contexte de
mesure, **10** classes de connaissances mécanistiques et **19** sources.

### 1.1 Trois refus codés, pas trois avertissements

**Un paramètre PK et une mesure fonctionnelle ne peuvent pas être le même
objet.** `PK_PARAMETERS` et `FUNCTIONAL_ENDPOINTS` sont disjoints par
assertion au chargement du module. `clint_in_vitro` porte des unités de
clairance et se mesure sans absorption, sans liaison et sans débit sanguin :
c'est un composant de la clairance métabolique, elle-même composant de CL.
**Conséquence immédiate et visible dans la matrice : 12 des 14 échecs en H1
sont dus à cette règle**, pas à un manque de puissance.

**L'appariement est affirmé, jamais déduit.** Le registre ne porte **aucun
champ d'appariement** — c'est exactement ce qu'il ne faut pas inférer. Chaque
jeu déclare le sien avec une **provenance** citant le document ou la notice
qui l'établit, et `unknown` rend toute analyse jointe inadmissible.

**Une méthode n'est pas une source de biologie.** `SOURCES` porte un `kind` :
`database`, `repository`, `annotation_tool`, `method`. `assert_is_a_database`
refuse DIABLO, mixOmics, l'elastic net, XGBoost/SHAP et Spearman. Appeler
DIABLO une source mettrait un choix statistique là où va un fait biologique.

### 1.2 Les connaissances mécanistiques ne sont pas réduites aux CYP

Votre remarque est encodée : dix classes, dont **transporteurs d'efflux**
(ABCB1, ABCG2, ABCC) et **d'absorption** (SLCO, SLC22, SLC47), **protéines de
liaison** (albumine, AAG, lipoprotéines), **partenaires rédox** (POR,
cytochrome b5), **récepteurs nucléaires** (PXR, CAR, AhR) et **machinerie
d'absorption** (CYP3A intestinal, transporteurs intestinaux). Chaque classe
porte le ou les paramètres PK qu'elle peut informer — `F`, `ka`, `V`, `CL`,
`CL_metabolic` — de sorte qu'une cible pharmacologique, qui n'informe aucun
paramètre PK, porte un tuple vide.

---

## 2. Recherche de données publiques (§4 du recentrage)

Six recherches, le 2026-11-12. **Les résultats négatifs sont conservés** :

| Dépôt | Requête | Retours | Sur cible |
|---|---|---|---|
| GEO | drug clearance AND expression profiling | 2 | **0** |
| GEO | liver AND cytochrome P450 activity AND Homo sapiens | 2 | **0** |
| GEO | pharmacokinetics AND transcriptome | 3 867 | surface inutilisable |
| **PRIDE** | pharmacokinetics | 25 | **1** |
| **PRIDE** | human liver microsomes / intrinsic clearance / hepatocyte CYP | 13 | **4** |
| Metabolomics Workbench | titre contenant « pharmacokinetic » | 1 | 0 (souris) |

**Le négatif GEO est informatif et vaut d'être retenu : les dépôts
transcriptomiques indexent l'essai, pas le phénotype.** Une clairance par
donneur n'apparaît presque jamais dans leurs métadonnées, donc une recherche
par mot-clé ne peut pas atteindre les jeux appariés. Les quatre retours étaient
vitamine B2/AhR, Bcrp/Pgp **chez le rat**, hépatocytes iPSC en bioréacteur et
un classifieur toxicogénomique.

**PRIDE est le filon**, pour une raison structurelle : les dépôts protéomiques
accompagnent les articles de DMPK, et ce sont ces articles qui portent la PK
dans leurs suppléments.

### 2.1 Le meilleur candidat trouvé

**PRIDE `PXD066482` — variabilité interindividuelle du métabolisme de
l'imatinib en microsomes hépatiques et hépatocytes primaires humains, phénotypes
CYP2C8 et CYP3A.** CC0, déposé 2025-07-23.

| | |
|---|---|
| médicament | **imatinib** |
| endpoint | **clairance intrinsèque** et formation de N-desméthyl-imatinib |
| couches | **génotype**, **concentration protéique**, **activité** |
| compartiments | **microsomes ET hépatocytes primaires** |
| appariement | **mêmes donneurs**, affirmé par la notice du dépôt |

**C'est la seule source inventoriée portant génotype, protéine, activité et
une clairance sur les mêmes unités.** C'est exactement le dessin de H3, et le
comparateur PGx que H3 exige est présent — ce qu'aucun autre jeu n'offre.

**Ce qui manque et doit être vérifié avant de l'appeler `READY` en pratique :**
le dépôt PRIDE porte la **protéomique** ; la clairance, les activités et les
génotypes sont dans le **supplément de l'article**, dont le caractère public
n'est **pas encore vérifié**. Le nombre de donneurs n'est pas indiqué dans les
métadonnées du dépôt. Le statut `READY` de la matrice porte donc sur le dessin,
pas sur la disponibilité effective des trois tables.

### 2.2 Les trois autres

| Accession | Ce que c'est | Pourquoi il compte | Ce qui bloque |
|---|---|---|---|
| `PXD057282` | tucatinib, **21 donneurs**, QTAP CYP2C8/3A4/3A5 + activités + déplétion | même dessin, autre médicament → **candidat de réplication** et non second jeu de découverte | **aucun génotype**, donc pas de comparateur PGx ; microsomes commerciaux, donc pas de covariables cliniques |
| `PXD060653` | protéome des enzymes de métabolisme sous **alcoolisme chronique** | une **perturbation** du protéome enzymatique, la forme dont H4 a besoin | entre groupes et non intra-individuel ; aucune PK ; appariement non affirmé |
| `PXD014131` | **homogénat, microsomes, culot, hépatocytes** du même organe | la seule source inventoriée touchant l'axe **« quel compartiment mesurer »** | aucun endpoint fonctionnel, aucune PK, appariement non affirmé |

---

## 3. La matrice (§4 et Phase 2)

15 jeux × 5 hypothèses = 75 cellules, chacune calculée et accompagnée de la
raison de son statut.

| Hypothèse | `READY` | `PARTIAL` | `NOT_EVALUABLE` |
|---|---|---|---|
| **H1** — omiques au-delà de PGx + clinique, pour la PK | **1** — ARIANA métoprolol 47 | 0 | 14 |
| **H2** — sélection automatique | **5** — les 4 jeux MaveDB + `PXD066482` | 1 — HLC | 9 |
| **H3** — fonctions biologiques vs PGx | **5** — les 4 MaveDB + `PXD066482` | 1 — HLC | 9 |
| **H4** — longitudinal | **0** | 0 | **15** |
| **H5** — posologie | **0** | 0 | **15** |

### 3.1 Deux cellules que j'avais d'abord annoncées `READY` à tort

**ARIANA en H2.** Un seul drapeau « comparateur disponible » par jeu ne peut
pas répondre à cinq questions différentes : ARIANA porte la base clinique et
PGx dont H1 a besoin, et aucune règle posologique à battre pour H5. Le drapeau
est devenu un dictionnaire par hypothèse. Et une contrainte documentée est
désormais codée : **ARIANA a exactement deux variables omiques candidates, et
choisir deux sur deux n'est pas une sélection.** H2 exige au moins cinq
candidats.

**ARIANA en H5.** Corrigée par la même modification : aucune règle posologique,
aucune atteinte de cible, 47 concentrations.

### 3.2 H4 et H5 sont vides, et la raison n'est pas la même

**H4** : aucun jeu inventorié n'a d'**occasions répétées**. Ce n'est pas un
problème de puissance ni d'accès, c'est l'absence du dessin.

**H5** : aucun jeu ne porte le comparateur minimal — **la règle posologique
déjà en usage** — ni une atteinte de cible. Une amélioration posologique ne
peut pas être revendiquée faute de l'objet à améliorer.

### 3.3 Ce que l'ancien registre aurait fait croire

Treize jeux étiquetés `H5` ancien (robustesse) et six `H4` ancien (« quand les
omiques nuisent »). En reprenant les étiquettes, la matrice aurait annoncé
**six jeux longitudinaux et treize jeux posologiques**. Il y en a **zéro de
chaque**.

---

## 4. Endpoints vérifiables, hypothèse par hypothèse

Ce que vous avez demandé avant toute nouvelle expérience. Chaque ligne est un
endpoint précis, un comparateur minimal et une condition de réfutation. **Rien
n'est lancé.**

### H1 — un seul jeu, et il est déjà gelé

| | |
|---|---|
| jeu | ARIANA, 47 patients, métoprolol |
| endpoint | **concentration** par patient (ni AUC ni clairance : la §5 du cahier interdit de l'appeler ainsi) |
| comparateur minimal | clinique + biologie + PGx chez les mêmes patients |
| couche omique | métabolomique ciblée, **2 variables** |
| réfutation | borne supérieure du gain sous la marge de 5 % |
| statut | **le contraste principal H1 Ridge M3 vs M1 est déjà gelé et ne doit pas être modifié.** Il n'y a pas de nouvelle expérience H1 à lancer sur les données inventoriées |

**Conséquence à énoncer clairement : H1, votre hypothèse prioritaire, n'a
qu'un seul véhicule public et il est déjà consommé.** Toute avancée en H1
demande de nouvelles données, pas une nouvelle méthode.

### H2 — cinq jeux, et la distinction de tâche est décisive

| | |
|---|---|
| jeux | NUDT15, CYP2C9, BRAF, ARAF (MaveDB) ; `PXD066482` sous réserve §2.1 |
| endpoint | score fonctionnel de variant, ou CL<sub>int</sub> pour l'imatinib |
| comparateur minimal | **Spearman** sur le même pool, et `k = 0` admissible |
| réfutation | aucun sélecteur ne bat Spearman d'un gain atteignant la marge |
| déjà mesuré | **0 cellule sur 11** à p ≤ 33 et à p = 2 000. La sophistication ne paie pas |
| **ce qui reste à faire** | **la tâche prospective** : décider de mesurer **sans accéder à la valeur** de la mesure choisie |

**La seconde tâche n'a pas de véhicule.** Elle exige **au moins deux couches
candidates et un endpoint fonctionnel indépendant** chez les mêmes unités.
Inventorié : 26 gènes MaveDB avec activité, 68 avec abondance, **4 avec les
deux**, et aucun avec deux **couches** distinctes. ARAF et BRAF offrent deux
**conditions** d'une même couche, ce qui est une autre question.

### H3 — la voie la plus praticable, et une correction de portée à intégrer

| | |
|---|---|
| jeux | les 4 MaveDB `READY` ; HLC `PARTIAL` ; `PXD066482` sous réserve |
| endpoints | activité enzymatique, CL<sub>int</sub>, score fonctionnel — **jamais appelés CL, F, V ou AUC** |
| comparateur minimal | **PGx seule** chez les mêmes unités |
| réfutation | borne supérieure du gain sur la PGx sous la marge |
| **non fait, et c'est le manque le plus net** | **aucun endpoint de transporteur, de liaison ni d'absorption n'a été testé.** Tout ce qui est fait porte sur des enzymes |

**Votre remarque « H3 n'est pas limité aux enzymes » n'est pas encore honorée
par les données.** Les dix classes mécanistiques sont au catalogue ; les
mesures ne couvrent que les phases I et II. `PXD066482` n'y change rien :
CYP2C8 et CYP3A sont des enzymes.

### H4 — rien, et ce qui serait nécessaire

| | |
|---|---|
| endpoint | **variation intra-individuelle** de CL, F ou exposition entre visites |
| comparateur minimal | variation intra-individuelle de l'état clinique |
| exigence de dessin | **occasions répétées chez les mêmes individus**, avec omiques à chaque occasion |
| inventorié | **aucun jeu** |
| à chercher | cohortes TDM longitudinales avec prélèvement omique répété ; `PXD060653` montre la **forme** (perturbation du protéome enzymatique) sans le dessin |

### H5 — rien, et ne pas le revendiquer

| | |
|---|---|
| endpoint | erreur de prédiction d'exposition et atteinte de cible |
| comparateur minimal | **la règle posologique déjà en usage** |
| inventorié | **aucun jeu** |
| position | à garder comme finalité, sans revendication anticipée |

---

## 5. Ce que je recommande, et ce que je ne fais pas

**Je ne lance rien.** Trois constats me paraissent devoir précéder toute
décision :

1. **H1 n'a qu'un véhicule public et il est consommé.** La priorité la plus
   haute du programme est bornée par les données, non par la méthode.
2. **H3 est la seule voie où du nouveau travail est possible aujourd'hui**, et
   `PXD066482` est le seul jeu inventorié portant le comparateur PGx qu'elle
   exige. Vérifier que le supplément de l'article est public est **la seule
   action de faible coût et de fort rendement** actuellement disponible.
3. **Le trou le plus intéressant est un trou de contenu, pas de méthode** :
   rien dans l'inventaire ne mesure un **transporteur**, une **liaison
   plasmatique** ou une **absorption**, alors que ce sont les mécanismes qui
   informent `F`, `ka` et `V` — les paramètres que la clairance métabolique ne
   touche pas.

**Prochaine action proposée, à votre arbitrage** : vérifier la disponibilité
publique du supplément de `PXD066482` et, si elle est acquise, pré-enregistrer
un protocole H3 avec la PGx comme comparateur minimal. C'est la première fois
que ce comparateur serait disponible.

---

## Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-12 | Phases 1 et 2. Catalogue structuré (16 entrées cliniques, 12 biologiques, 10 descripteurs de contexte, 10 classes mécanistiques dont transporteurs, liaison et absorption, 19 sources typées). Trois refus codés : PK et fonctionnel disjoints, appariement affirmé avec provenance et jamais déduit, méthode distincte de source. Conflit d'étiquettes `H1`–`H5` ancien/nouveau documenté et gardé par assertion. Six recherches publiques dont les négatifs GEO. Quatre candidats PRIDE, dont `PXD066482` qui porte génotype + protéine + activité + clairance sur les mêmes donneurs. Matrice 15 × 5 calculée : **H1 1 `READY`, H2 et H3 5 `READY`, H4 et H5 0**. Deux cellules annoncées `READY` à tort puis corrigées par un comparateur par hypothèse et un plancher de pool pour H2. Aucune expérience lancée. |
