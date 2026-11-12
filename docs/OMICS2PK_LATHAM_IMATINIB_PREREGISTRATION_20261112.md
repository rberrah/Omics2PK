# Pré-enregistration : la protéomique ajoute-t-elle à la PGx sur une clairance mesurée ?

**Version 2.0.0 — 2026-11-12.** Statut : **EXÉCUTÉ UNE FOIS (§8). Rien n'est établi : la seule cellule atteignant la marge échoue à sa propre condition de robustesse pré-déclarée.** La §7.1 était **fausse** et n'est pas effacée : la §7.4 donne la cause réelle, un bug de mon extraction. Cohorte primaire **n = 13**. Aucun modèle n'a encore été ajusté.

Ce document est commité **avant l'extraction de la moindre valeur par donneur**.
Seules la **structure** du tableau source, les effectifs et le régime de licence
ont été vérifiés ; aucune clairance, aucune concentration et aucune activité
n'a été regardée. Lire les valeurs avant d'écrire le protocole le rendrait sans
valeur.

---

## 1. Pourquoi ce jeu, et ce qu'il est

Latham *et al.*, *Drug Metabolism and Disposition* 2025,
DOI `10.1016/j.dmd.2025.100196`, `PMC12990092`, protéomique brute déposée sous
**`PXD066482`** (CC0).

**C'est la première source inventoriée portant les quatre couches chez les
mêmes unités.** Structure du **tableau 4**, vérifiée le 2026-11-12 :

| Rôle | Colonnes |
|---|---|
| clé d'appariement | **Donor ID** |
| covariables cliniques | âge, sexe, origine rapportée, **IMC** |
| **PGx** | génotypes **CYP2C8**, **CYP3A4**, **CYP3A5**, et statut métaboliseur CYP3A5 |
| **protéomique** | concentrations **CYP2C8**, **CYP3A4**, **CYP3A5**, **CYP3A** |
| **activités** | CYP2C8 (N-déséthyl-amodiaquine), CYP3A (1-OH-midazolam), CYP3A (6β-OH-testostérone) |
| **endpoints** | **clairance de l'imatinib** (µL/10⁶ cellules/min), formation de N-desméthyl-imatinib |

**14 lignes de données.** Le tableau 1 porte 22 lignes pour les microsomes, mais
**sans génotype individuel** : les microsomes ne peuvent donc pas porter un
contraste dont la PGx est le comparateur.

### 1.1 Le régime de licence, et ce qu'il autorise

`isOpenAccess: N`, `license: None`, `restricted-by: pmc`. Le bloc de licence du
dépôt PMC indique : *« This file is available for text mining. It may also be
used consistent with the principles of fair use under the copyright law. »*

**Donc : extraction pour analyse autorisée, redistribution non.** Les valeurs
par donneur sont écrites **hors du dépôt**, comme pour HLC ; seuls des agrégats,
des pertes et des intervalles entrent dans les artefacts. **Aucune ligne par
donneur ne figurera dans un fichier commité.**

---

## 2. Ce que ce jeu n'est pas, et la phrase à ne jamais écrire

**La clairance mesurée ici est une clairance intrinsèque en hépatocytes
primaires.** Elle se mesure sans absorption, sans liaison plasmatique, sans
débit sanguin hépatique et sans excrétion rénale.

> **Elle n'est pas une clairance systémique, et aucun résultat de ce protocole
> ne sera extrapolé à CL, à F, à V ou à une AUC chez des patients.**

Le catalogue le code plutôt que de l'espérer : `clint_in_vitro` est dans
`FUNCTIONAL_ENDPOINTS`, disjoint de `PK_PARAMETERS` par une assertion au
chargement du module. Un contraste H1 déclaré sur cet endpoint est donc
**refusé automatiquement**, et c'est voulu.

**Conséquence de nomenclature, pour éviter l'ambiguïté :** le contraste ci-
dessous est appelé **H1-in-vitro**. Ce n'est **pas** l'hypothèse H1 du
programme, qui exige un paramètre PK chez des individus. H1 garde son unique
véhicule, ARIANA, et son contraste gelé.

---

## 3. Les contrastes, déclarés

### 3.1 Primaire — H1-in-vitro

> **clinique + PGx + protéomique** contre **clinique + PGx**, endpoint
> **clairance de l'imatinib**, chez les donneurs complets.

* **clinique** : âge, sexe, IMC.
* **PGx** : génotypes CYP2C8, CYP3A4, CYP3A5 et statut métaboliseur CYP3A5,
  encodés en dosage d'allèle ou en statut, **gelé avant exécution** : le
  statut métaboliseur CYP3A5 tel que rapporté par les auteurs, et pour
  CYP2C8/CYP3A4 un indicateur « non-référence » par gène.
* **protéomique** : concentrations CYP2C8, CYP3A4, CYP3A5. **CYP3A total est
  exclu** du bloc protéomique : c'est une somme de CYP3A4 et CYP3A5 et son
  inclusion ferait entrer deux fois la même mesure.
* Lecture : échelle à six barreaux, marge **5 %**, appliquée à la borne
  **inférieure**.

### 3.2 Secondaire — H3

> **protéomique + PGx** contre **PGx seule**, endpoints **activité CYP2C8** et
> **activité CYP3A**.

Deux endpoints d'activité, donc **taille de famille 2** sur ce bras. L'activité
CYP3A est mesurée par deux sondes (midazolam et testostérone) : **le midazolam
est l'endpoint déclaré**, la testostérone est une sensibilité, et les deux ne
comptent pas comme deux validations indépendantes.

### 3.3 Déclaré non évaluable — H2

**H2 est `NOT_EVALUABLE` sur ce jeu.** Trois protéines CYP sont quantifiées, et
le catalogue exige un pool d'au moins **cinq** candidats indépendants : choisir
trois sur trois n'est pas une sélection. Ce seuil a été fixé avant ce jeu, à
partir du cas ARIANA, et il n'est pas assoupli pour l'occasion.

### 3.4 Secondaire — la formation de métabolite

La formation de N-desméthyl-imatinib est un second endpoint fonctionnel,
déclaré **secondaire** et non une réplication : c'est la même incubation chez
les mêmes donneurs.

---

## 4. Tout ce qui est gelé

| Paramètre | Valeur |
|---|---|
| cohorte | donneurs d'hépatocytes primaires du tableau 4, **complets sur les trois concentrations protéiques** |
| n | **13** : les donneurs complets sur le bloc protéomique et l'endpoint. GKP sort **par les données** (§7.4), jamais par une liste |
| analyse d'influence | **obligatoire** : 13 réexécutions, chaque donneur retiré à son tour |
| appris dans le pli | centrage et mise à l'échelle, **et rien d'autre** |
| fixé a priori | ordre CYP3A5 PM < IM < NM ; diplotype de référence `*1/*1` |
| contrôle d'abstention | **synthétique et obligatoire**, où `k = 0` est la bonne réponse (§5.1 bis) |
| validation | **leave-one-donor-out**, chaque donneur étant son propre pli |
| unité de rééchantillonnage | le **donneur** — et une ligne *est* un donneur, donc bootstrap par grappes et par lignes coïncident, ce qui est dit plutôt que supposé |
| estimateur | **ridge**, pénalité 1,0, identique à tous les bras |
| perte | erreur absolue |
| marge | **5 %** du bras de référence |
| taille de famille | **4** : primaire, deux activités, formation de métabolite |
| graine | 20261112 |
| transformation | log de la clairance et des concentrations, déclarée avant exécution ; les valeurs nulles ou négatives seraient consignées en §7 |
| bras d'abstention | `k = 0` admissible partout, c'est-à-dire clinique + PGx seules |

---

## 5. La prédiction, et ce qui réfute

| | Prédit |
|---|---|
| **primaire** | **`NO_DETECTABLE_GAIN_UNDERPOWERED`** |
| H3, activité CYP2C8 | gain positif de la protéomique sur la PGx, marge probablement indécise |
| H3, activité CYP3A | idem |
| nombre de plis choisissant `k = 0` | au moins 1 sur 13 |

**La prédiction principale est un non-résultat, et c'est délibéré.** Avec
**13 donneurs**, un intervalle apparié sur une différence de pertes sera large ;
l'issue la plus probable est qu'il contienne zéro **et** dépasse la marge, ce
qui ne répond à rien. **Cette issue ne sera pas rapportée comme un soutien, ni
comme une réfutation.**

### 5.1 Ce qui réfute quoi

* **Borne supérieure du gain primaire sous la marge** → réfutation propre : la
  protéomique n'ajoute pas à la PGx d'un montant qui compterait, **sur cet
  endpoint et à cet effectif**.
* **Gain primaire atteignant la marge** → résultat positif, et il resterait
  **exploratoire** : 13 donneurs, un médicament, un laboratoire, aucune
  réplication indépendante.
* **Intervalle contenant zéro et dépassant la marge** → **le protocole ne
  répond pas.** Déclaré ici pour qu'il ne soit pas présenté autrement.
* **La PGx seule battant clinique + PGx + protéomique de façon détectable** →
  la protéomique nuit, et il faudra le dire.
* ~~**Aucun pli ne choisissant `k = 0`** → défaut d'implémentation.~~
  **Corrigé en §7.5 : cette règle était fausse.**

### 5.1 bis Deux règles corrigées avant exécution

**« Aucun pli ne choisit `k = 0` » n'est pas un défaut.** Ma §5.1 l'affirmait et
c'est faux : cela peut simplement signifier que le modèle **retient** les
variables supplémentaires, ce qui est une réponse et non une panne.
L'abstention ne se vérifie pas en observant qu'elle se produit sur des données
réelles, mais par un **contrôle synthétique où elle est la bonne réponse** —
une cohorte dont les variables supplémentaires sont du bruit pur. Ce contrôle
est obligatoire avant que le compte d'abstentions de ce protocole ne soit
interprété dans un sens ou dans l'autre.

**Tout ce qui est appris l'est sur les plis d'entraînement seuls.** À 13
donneurs, les catégories PGx rares sont instables : **un seul porteur
`CYP3A4 *1/*15`** et **un seul `NM` CYP3A5**. Dans le pli qui retire ce
donneur, la catégorie est **absente de l'entraînement**. Le codage est donc
conçu pour que cela produise une colonne **constante** et non une erreur, et :

| Fixé a priori, donc sans fuite | Appris dans le pli |
|---|---|
| l'ordre phénotypique CYP3A5 **PM < IM < NM** | le centrage |
| le diplotype de référence par gène (`*1/*1`) | la mise à l'échelle |

**Rien d'autre n'est appris**, et un test vérifie que retirer un donneur de
l'entraînement change bien le centrage — sinon la standardisation ne serait pas
dans le pli du tout.

### 5.2 Ce que ce protocole ne pourra pas établir, quoi qu'il montre

**Il ne peut pas démontrer un gain généralisable.** Treize donneurs, un
substrat, deux enzymes, un seul laboratoire, et un endpoint in vitro. C'est une
**validation biologique supplémentaire**, pas la démonstration clinique que
cherche Omics2PK — et la différence est la distance entre une clairance
hépatocytaire et une clairance systémique.

---

## 6. Ce qui sera rapporté quoi qu'il arrive

Le contraste primaire, son intervalle, son verdict, le compte d'abstentions,
les deux contrastes H3 et le nombre de donneurs effectivement complets. Aucun
bras ne sera ajouté et aucune marge déplacée après lecture des valeurs ; tout
écart figurera en §7 **avant** qu'un contraste ne soit cité.

---

## 8. Résultats

Exécuté une fois, le 2026-11-12, après validation du contrôle d'abstention.
`reports/data/latham_imatinib_protocol.json`. **Cohorte n = 13**, GKP exclu
par les données. Aucun comparateur, aucune marge et aucun modèle n'a été
changé.

### 8.1 Les cinq contrastes

| Endpoint | gain | verdict | influence |
|---|---|---|---|
| **primaire — clairance de l'imatinib** | **−14,25 %** | `NO_DETECTABLE_GAIN_UNDERPOWERED` | **5 / 13 changent** |
| H3 — activité CYP2C8 | +26,08 % | non détectable | **0 / 13** |
| **H3 — activité CYP3A (midazolam)** | **+54,14 %** | **`GAIN_REACHES_MARGIN`** | **1 / 13 change** |
| secondaire — formation de métabolite | +4,70 % | non détectable | 1 / 13 |
| sensibilité — CYP3A (testostérone) | +56,72 % | `GAIN_REACHES_MARGIN` | **5 / 13 changent** |

### 8.2 Ce que ma propre règle d'influence impose

La §7.2 déclarait : *« si au moins une exclusion change le verdict, le
contraste est rapporté comme non concluant — pas comme positif avec une
réserve »*, et la §7.4.2 a maintenu cette lecture. Elle s'applique donc :

| Endpoint | Lecture contraignante |
|---|---|
| primaire | **non concluant**, et déjà non détectable |
| H3 CYP2C8 | le **non-détectable est robuste** : aucune exclusion ne le change |
| **H3 CYP3A midazolam** | **non concluant**, malgré +54,14 % et une marge atteinte, parce que retirer **BGY** l'annule |
| métabolite | **non concluant** |
| sensibilité testostérone | **non concluant** |

> **Rien n'est établi par ce protocole.** La seule cellule atteignant la marge
> échoue à sa propre condition de robustesse.

**Et je ne relâche pas la règle maintenant.** Elle a été écrite pour neutraliser
une contradiction de source qui n'existait pas, et elle est probablement **trop
stricte** pour la fragilité ordinaire d'un effectif de 13 : à cette taille, un
donneur qui fait basculer un verdict est attendu, pas alarmant. Mais
l'assouplir **après avoir vu qu'elle mord** est exactement le geste que tout
cet appareil existe pour empêcher. Elle est donc contraignante ici, et la §9
l'enregistre comme une leçon de conception pour la prochaine pré-enregistration.

### 8.3 Les prédictions, et celle qui a été dépassée

| §5 prédisait | mesuré | |
|---|---|---|
| **primaire : `NO_DETECTABLE_GAIN_UNDERPOWERED`** | exactement cela | **tenu** |
| H3 CYP2C8 : gain positif, marge indécise | +26,08 % positif, **non détectable** | partiel |
| H3 CYP3A : gain positif, marge indécise | **+54,14 %, marge atteinte** | **dépassé**, puis annulé par l'influence |
| au moins 1 pli à `k = 0` | **4 plis sur 13** au primaire, **0** sur les bras H3 | voir §8.4 |

Le gain primaire est **négatif** : la protéomique dégrade la prédiction de la
clairance de 14 % par rapport à clinique + PGx seules. Les **cinq** exclusions
qui changent le verdict le font toutes dans le même sens —
`REFUTED_GAIN_BELOW_MARGIN` — donc la tendance du primaire est cohérente :
**sur cet endpoint et à cet effectif, ajouter la protéomique n'aide pas et
nuit peut-être.**

### 8.4 La règle corrigée a servi immédiatement

Les bras H3 ne choisissent `k = 0` dans **aucun** pli. Mon §5.1 d'origine en
aurait fait un défaut d'implémentation ; la correction de la §5.1 bis dit que
cela signifie que le modèle **retient** les protéines. C'est cohérent avec les
gains les plus élevés du run étant précisément sur ces bras, et le contrôle
synthétique avait déjà démontré que l'abstention fonctionne quand elle est la
bonne réponse. **Sans cette correction, j'aurais signalé un bug là où il y a un
résultat.**

Et le primaire, lui, s'abstient dans **4 plis sur 13** — l'abstention se
produit donc bien, sur l'endpoint où la protéomique n'aide pas.

### 8.5 Ce que ce run autorise à dire

| Question | Réponse |
|---|---|
| la protéomique ajoute-t-elle à clinique + PGx pour la clairance hépatocytaire ? | **non concluant.** Point négatif à −14,25 %, non détectable, et 5 exclusions sur 13 le basculent |
| la protéomique ajoute-t-elle à la PGx pour l'activité CYP3A ? | **non concluant**, malgré +54,14 % et la marge atteinte : un seul donneur l'annule |
| pour l'activité CYP2C8 ? | **non**, et c'est le seul verdict **robuste** du run |
| H2 ? | **`NOT_EVALUABLE`**, trois protéines contre un plancher de cinq |
| une conclusion clinique ? | **aucune.** Clairance hépatocytaire, 13 donneurs, un substrat, un laboratoire |

---

## 9. Leçons de conception pour la prochaine pré-enregistration

1. **Une condition de robustesse binaire est trop stricte à petit n.** « Toute
   exclusion qui change le verdict rend non concluant » est raisonnable pour
   neutraliser un donneur douteux, et punitive comme critère général : à 13
   donneurs, 1 bascule sur 13 est le bruit attendu. La prochaine fois, déclarer
   un **seuil** — par exemple « plus de 20 % des exclusions » — et le déclarer
   **avant** de savoir s'il mord.
2. **Vérifier les planchers de taille du harnais avant de geler un effectif.**
   `pool_arm_loss` refusait n < 20 et exigeait 10 lignes par pli interne ; à
   n = 13 le protocole était inexécutable et je l'ai découvert après l'avoir
   gelé. Les deux planchers sont désormais des paramètres déclarés qui
   apparaissent dans l'artefact.
3. **Décider la manquance par classe de caractères, jamais par une liste.**
   Voir §7.4.

---

## 7. Écarts

Section en ajout seul. **Consigné avant tout ajustement de modèle.**

### 7.1 ~~La source se contredit sur l'effectif protéomique~~ — **FAUX, voir §7.4**

> **Cette section est conservée telle quelle et elle est erronée.** Sa
> conclusion — « la contradiction n'est pas résoluble depuis les sources » —
> était un artefact de **mon** extraction, non une propriété de l'article.
> Elle reste lisible parce qu'effacer une interprétation fausse rendrait la
> trace inutilisable : ce qui suit est ce que je croyais le 2026-11-12 avant
> la vérification externe. La cause réelle est en §7.4.

### 7.1 bis Ce que je croyais, et qui était faux

Le §4 annonçait « n attendu **13 sur 14** d'après le rapport externe ; **à
recompter à l'extraction** et à consigner si l'écart apparaît ». Il apparaît,
et il est dans la source elle-même.

| Source | Ce qu'elle dit |
|---|---|
| **texte de l'article** | « P450 activity measures were from n=14 donors, and **P450 protein quantitation was from n=13 donors** » |
| **tableau 4** | **14 lignes**, et les trois concentrations CYP2C8, CYP3A4, CYP3A5 sont **renseignées sur les 14**, sans aucune cellule vide ni marqueur de manquant |

Vérifications faites pour trancher, en vain :

* **aucune ligne n'est un pool.** Les 14 identifiants sont des codes de donneur
  à trois lettres ; le pool à 5 donneurs mentionné dans les méthodes ne figure
  pas dans le tableau.
* **aucune colonne du tableau 4 n'a de valeur manquante**, sur les 18.
* Le texte ne nomme pas le donneur concerné.

**La contradiction n'est donc pas résoluble depuis les sources disponibles.**

### 7.2 Comment je la traite : une analyse d'influence, pas un choix arbitraire

Retirer « le » donneur exigerait de savoir lequel, et le choisir au hasard
introduirait une décision arbitraire dans le contraste principal.

**Déclaré ici, avant tout ajustement :**

1. **Le contraste principal tourne sur les 14 donneurs** du tableau 4, qui est
   la seule source portant des valeurs par donneur.
2. **Une analyse d'influence accompagne obligatoirement le résultat** : le
   contraste est recalculé **14 fois, chaque donneur retiré à son tour**, et le
   verdict de chacune des 14 exclusions est rapporté.
3. **Lecture déclarée d'avance :** si **aucune** des 14 exclusions ne change le
   verdict, la contradiction de la source est **sans effet** sur la conclusion
   et c'est dit ainsi. Si **au moins une** le change, alors le résultat dépend
   d'un donneur dont la quantification protéique n'est peut-être pas soutenue
   par le texte de l'article, et **le contraste est rapporté comme non
   concluant** — pas comme positif avec une réserve.

Cette analyse est presque gratuite ici : la validation est déjà
leave-one-donor-out, donc 14 réexécutions coûtent 14 fois un ajustement sur
13 lignes.

### 7.4 La cause réelle : mon extraction lisait un tiret comme une valeur

**Il n'y a aucune contradiction dans la source.** Le tableau 4 **porte bien**
des valeurs manquantes, marquées par un tiret, et mon contrôle de complétude
les a comptées comme présentes.

| | |
|---|---|
| donneur concerné | **GKP**, et lui seul |
| cellules absentes | **CYP3A4**, **CYP3A5** et **CYP3A total** |
| ce que GKP garde | CYP2C8 à **1,3 pmol/mg**, les trois activités, la clairance, le métabolite, l'âge, l'IMC |
| caractère employé | **`U+2212 MINUS SIGN`**, et non le trait d'union ASCII |

**Pourquoi mon contrôle a échoué.** Il testait l'appartenance à une liste de
jetons écrite à la main, qui contenait le trait d'union, le tiret demi-cadratin
et le cadratin — **mais pas `U+2212`**, qui est classé `Sm` (symbole
mathématique) et non `Pd` (tiret). Les trois cellules ont donc été déclarées
présentes, et j'en ai tiré que l'article se contredisait.

**Le décompte des auteurs était exact depuis le début :** 14 donneurs pour
l'activité, **13 pour le panel protéomique**.

#### 7.4.1 Le correctif, et pourquoi il n'est pas une liste plus longue

Allonger la liste de jetons aurait reproduit la même faute à la prochaine
table. La manquance est désormais décidée par **catégorie Unicode** —
`unicodedata.category == "Pd"` couvre tous les tirets, et `U+2212` est nommé
explicitement **à côté**, parce que c'est lui qui a causé le bug. Trois règles
tiennent le reste :

1. **Un tiret est manquant ; un zéro est une mesure.** `parse_number` rend
   `None` pour une manquance et `0.0` pour un zéro, et un test vérifie qu'une
   manquance ne devient jamais `0.0`.
2. **Un jeton non reconnu lève une exception**, au lieu d'être traité comme
   manquant. Élargir silencieusement l'ensemble des manquants est exactement
   ce qui a produit la première erreur.
3. **Les colonnes numériques sont vérifiées après conversion** : chaque valeur
   est `None` ou un flottant, jamais une chaîne qui aurait survécu.

33 tests, dont une non-régression nommée sur GKP et la vérification qu'**aucun
identifiant de donneur n'est codé en dur comme exclusion** : GKP sort par les
données, pas par une liste.

#### 7.4.2 La cohorte, recomptée par bloc

| Bloc | n complet | manquant |
|---|---|---|
| clinique (âge, IMC) | **14** | — |
| trois activités | **14** | — |
| **trois concentrations CYP** | **13** | **GKP** |
| CYP3A total | 13 | GKP |
| clairance imatinib | **14** | — |
| formation de métabolite | **14** | — |

**Cohorte primaire : n = 13**, définie comme les donneurs complets sur le bloc
protéomique **et** l'endpoint, sans exclusion arbitraire. **GKP reste
disponible** pour tout bras n'exigeant que CYP2C8 ou les activités.

**L'analyse d'influence de la §7.2 est conservée** et porte désormais sur les
**13** donneurs complets : 13 réexécutions, chaque donneur retiré à son tour.
Sa justification change — ce n'est plus une contradiction de source à neutraliser
mais la fragilité ordinaire d'un effectif de 13 — et sa lecture reste la même.

### 7.3 Extraction : ce qui est sorti et où

Tableau 4 extrait depuis le XML PMC par le canal eutils vers
`Latham2025_imatinib_PHH/latham2025_table4_phh_per_donor.csv`, **hors du
dépôt**. 18 colonnes, 14 lignes, largeur uniforme. Conformément à la §1.1,
**aucune ligne par donneur n'entrera dans un artefact commité.**

---

## Change log

| Version | Date | Change |
|---|---|---|
| **2.0.0** | **2026-11-12** | **Exécuté une fois. Rien n'est établi.** Primaire : la protéomique ajoutée à clinique + PGx donne **−14,25 % [−0,2584, +0,0487]** `NO_DETECTABLE_GAIN_UNDERPOWERED` sur la clairance hépatocytaire — **exactement la prédiction de la §5** — avec **5 exclusions sur 13** qui le basculent, toutes vers `REFUTED_GAIN_BELOW_MARGIN`. H3 activité CYP3A (midazolam) atteint la marge à **+54,14 % [+0,1526, +0,8017]**, mais **retirer le donneur BGY l'annule**, donc ma règle d'influence pré-déclarée le rend **non concluant** — et je ne la relâche pas après avoir vu qu'elle mord. H3 activité CYP2C8 : +26,08 %, non détectable, et **0 exclusion sur 13** ne le change : **le seul verdict robuste du run**. La correction de la §5.1 bis a servi immédiatement : les bras H3 ne s'abstiennent dans aucun pli, ce que l'ancienne règle aurait déclaré bug alors que c'est un résultat. §9 ajoutée : trois leçons de conception, dont qu'une condition de robustesse binaire est punitive à petit n et doit être un seuil déclaré d'avance. |
| 1.2.0 | 2026-11-12 | **Ma §7.1 était fausse, et la cause est un bug de mon extraction (§7.4).** Le tableau 4 porte bien trois valeurs manquantes, pour le seul donneur **GKP** (CYP3A4, CYP3A5, CYP3A total), marquées par **`U+2212 MINUS SIGN`** — absent de ma liste de jetons écrite à la main, parce qu'il est classé `Sm` et non `Pd`. **Le décompte des auteurs était exact** : 14 pour l'activité, 13 pour la protéomique. Le correctif décide la manquance par **catégorie Unicode**, rend `None` pour un tiret et `0.0` pour un zéro, **lève** sur un jeton inconnu, et vérifie les colonnes après conversion ; 33 tests dont une non-régression sur GKP et la garantie qu'aucun identifiant n'est codé en dur. **Cohorte primaire n = 13**, GKP exclu par les données et conservé pour les bras n'exigeant que CYP2C8 ou les activités ; analyse d'influence sur les 13. Deux règles corrigées (§5.1 bis) : **l'absence de pli à `k = 0` n'est pas un défaut** et se vérifie par un contrôle synthétique, et **tout ce qui est appris l'est dans le pli** — seuls le centrage et l'échelle, l'ordre CYP3A5 et le diplotype de référence étant fixés a priori. §7.1 **conservée et barrée**. |
| 1.1.0 | 2026-11-12 | **Écart consigné avant tout ajustement (§7).** La source se contredit : son texte dit que la quantification protéique porte sur **n = 13** donneurs, et son tableau 4 porte **14 lignes dont les trois concentrations sont renseignées sans aucun manquant**. Vérifié qu'aucune ligne n'est le pool à 5 donneurs, qu'aucune des 18 colonnes n'a de valeur manquante, et que le texte ne nomme pas le donneur concerné : **la contradiction n'est pas résoluble depuis les sources**. Traitée non par un retrait arbitraire mais par une **analyse d'influence obligatoire sur les 14 exclusions**, avec une lecture déclarée d'avance : si au moins une exclusion change le verdict, le contraste est **non concluant** et non positif-avec-réserve. Tableau 4 extrait **hors du dépôt**. |
| 1.0.0 | 2026-11-12 | Déclaré, **aucune valeur lue**. Structure du tableau 4 vérifiée : 14 donneurs d'hépatocytes primaires avec Donor ID, âge, sexe, IMC, génotypes CYP2C8/3A4/3A5, statut métaboliseur CYP3A5, concentrations CYP2C8/3A4/3A5/3A, trois activités et la clairance de l'imatinib. Licence PMC : extraction pour analyse autorisée, **redistribution non**, donc valeurs par donneur hors du dépôt. Primaire **H1-in-vitro** nommé ainsi pour qu'il ne soit pas confondu avec H1 du programme, qui exige un paramètre PK chez des individus. CYP3A total exclu du bloc protéomique comme somme de CYP3A4 et CYP3A5. **H2 déclaré `NOT_EVALUABLE`** : trois candidats contre un plancher de cinq fixé avant ce jeu. Validation leave-one-donor-out, marge 5 %, famille 4. **La prédiction principale est un non-résultat**, et la §5.1 déclare d'avance qu'un intervalle contenant zéro et dépassant la marge ne répond à rien. |
