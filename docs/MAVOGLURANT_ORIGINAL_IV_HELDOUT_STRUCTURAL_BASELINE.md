# Mavoglurant original humain : comparateur IV 1C/2C, validation inter-sujets

**Étape technique E3 préparatoire ; aucune preuve omique→PK ni modèle PopPK/NLME.**

## Donnée d'entrée authentifiée

La source est `nlmixr2data::mavoglurant`, commit `f2cfb01ade88d4d30e66f56b19efce2a2c6e26bd`, source `mavoglurant.rda` de SHA256 `ba1d75e20bcc09f20a55a05a9e4424244efa15b2dfcef6079cd27385d5876174`. Chargement en mémoire éphémère CI uniquement ; aucune concentration individuelle, aucun ID, aucune dose par personne ni donnée clinique brute déposée dans le repo ou dans les logs.

Le contrat source antérieur impose 120 patients, 198 occasions et 2 480 lignes EVID=0 candidates, élimine MDV=1/conc. non quantifiée, et abstient 2 occasions à doublons temporels sans créer de pseudo-moyennes. Il reste **196 profils admissibles/2 403 mesures**, dont 76 sujets avec deux occasions. La nouvelle évaluation emploie exactement ces 196 profils, en revalidant les comptages, la dose AMT et la durée d'administration AMT/RATE de chaque occasion. Aucun échantillon t=0 quantifié n'est attendu.

## Hypothèse et modèle

- **1C IV perfusé linéaire**, paramètres positifs fixes CL (L/h) et V (L). La concentration est la convolution analytique de la perfusion de durée effectivement indiquée dans le dossier EVID=1 avec la décroissance exponentielle.
- **2C IV perfusé linéaire**, paramètres positifs fixes CL, Vc, Q et Vp, solution biexponentielle de la masse centrale pendant/après la perfusion.
- Conversion explicite de mg/L vers ng/mL ; temps depuis la dose en heures. La source ne doit **pas** être traitée en absorption orale, bolus instantané ou perfusion de durée arbitraire.
- Ajustement numérique de paramètres **fixes partagés** en minimisant les erreurs de concentration en échelle log sur l'entraînement, avec poids normalisés par sujet et multi-initialisations. Le code **n'est pas une estimation populationnelle NLME à effets aléatoires** : pas d'ETA, de variance interindividuelle ou interoccasion, d'EBE, ni de SAEM/FOCEi.

## Anti-fuite et métrique

La répartition déterministe SHA256 en cinq folds affecte **un sujet biologique et ses deux occasions au même fold**. Chaque fold fait l'ajustement sur les sujets d'entraînement seulement ; les sujets tests ne contribuent pas au choix de paramètres et ne font l'objet d'aucun réajustement post-hoc. La métrique principale est la **racine de l'erreur quadratique moyenne en log concentration, calculée à poids égal par personne**, agrégée sur les plis. Le biais moyen log signé et la différence 2C−1C sont descriptifs, pas une preuve statistique de supériorité.

Les unités, échantillons et le split sont testés sur données synthétiques sans inclure de données originales. La CI optionnelle `pk-structural` exécute ensuite la validation sur l'original public épinglé. Toute variation SHA, effectifs et événements doit conduire à l'abstention.

## Résultat obtenu sur **la source humaine réelle**

Le workflow [GitHub Actions #38968850119](https://github.com/rberrah/Omics2PK/actions/runs/38968850119) a réussi le 16 octobre 2026, avec **6 tests analytiques**, **1 test de fit synthétique**, puis le calcul sur le fichier source original dont le SHA-256 a été confirmé.

| Score hors-échantillon, par sujet humain (5 folds) | 1C | 2C |
|---|---:|---:|
| Racine de l'erreur quadratique moyenne des log concentrations observées | **0,786014** | **0,383383** |
| Biais signé moyen par sujet en log concentration | −0,000563 | −0,000180 |
| Sujets biologiques jamais vus pendant l'ajustement (cumul des plis) | 120 | 120 |

Différence de RMSLE 2C−1C = **−0,402631**. Le 2C présente une meilleure capacité prédictive que le 1C dans cette comparaison descriptive interne patient-disjointe. L'ajustement n'a reçu aucun ID de patient tenu à l'écart, aucune concentration de test et aucune omique. La proximité des biais globaux avec zéro ne démontre **pas** une absence de biais dans les sous-groupes ni une calibration des intervalles individuels ; aucune erreur AUC individualisée et aucun intervalle d'incertitude n'ont été évalués dans cette itération.

**Limitation critique :** la métrique log pondère les concentrations tardives de faible amplitude, les modèles sont sans IIV/IOV ni erreur résiduelle explicitement estimée, et il s'agit de validation **interne entre individus**, pas d'une cohorte externe. Ne pas employer la diminution d'erreur 1C→2C comme preuve de contribution biologique des signatures omiques.

## Interprétation et suite

Un résultat CI vert prouve la reproductibilité technique d'un comparateur structural *pooled* sur de vrais C(t) humains ; ce n'est **pas** la preuve d'une prédiction PK par des omiques, d'un PopPK individuel, ni d'une validation externe indépendante. Les comparateurs E3 requis ensuite sont (1) NLME / SAEM clinique à effets aléatoires et interoccasion, (2) covariables cliniques et PGx disponibles, (3) effets omiques uniquement chez les **mêmes** patients avec horodatage antérieur à l'horizon prédit, et (4) test externe/temps futur avec incertitude de prédiction et BLQ.

Commande :

```sh
python -m pip install -e ".[pk-structural]"
python -m unittest discover -s tests -p 'test_iv_pk_structural.py'
python -m unittest discover -s tests -p 'test_mavoglurant_structural_fit.py'
python -m scripts.benchmark_mavoglurant_iv_structural
```

Résultats agrégés et SHA original vérifiés dans le workflow cité ci-dessus. Les sorties individuelles ne sont pas publiées.
