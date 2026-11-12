# Pré-enregistration : la protéomique ajoute-t-elle à la PGx sur une clairance mesurée ?

**Version 1.0.0 — 2026-11-12.** Statut : **DÉCLARÉ, AUCUNE VALEUR LUE.**

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
| n attendu | **13 sur 14** d'après le rapport externe ; **à recompter à l'extraction** et à consigner si l'écart apparaît |
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
* **Aucun pli ne choisissant `k = 0`** → l'abstention ne fonctionne pas à ce
  très petit effectif, et c'est un défaut d'implémentation à corriger, non un
  résultat.

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

## 7. Écarts

*Aucun. Section en ajout seul.*

---

## Change log

| Version | Date | Change |
|---|---|---|
| 1.0.0 | 2026-11-12 | Déclaré, **aucune valeur lue**. Structure du tableau 4 vérifiée : 14 donneurs d'hépatocytes primaires avec Donor ID, âge, sexe, IMC, génotypes CYP2C8/3A4/3A5, statut métaboliseur CYP3A5, concentrations CYP2C8/3A4/3A5/3A, trois activités et la clairance de l'imatinib. Licence PMC : extraction pour analyse autorisée, **redistribution non**, donc valeurs par donneur hors du dépôt. Primaire **H1-in-vitro** nommé ainsi pour qu'il ne soit pas confondu avec H1 du programme, qui exige un paramètre PK chez des individus. CYP3A total exclu du bloc protéomique comme somme de CYP3A4 et CYP3A5. **H2 déclaré `NOT_EVALUABLE`** : trois candidats contre un plancher de cinq fixé avant ce jeu. Validation leave-one-donor-out, marge 5 %, famille 4. **La prédiction principale est un non-résultat**, et la §5.1 déclare d'avance qu'un intervalle contenant zéro et dépassant la marge ne répond à rien. |
