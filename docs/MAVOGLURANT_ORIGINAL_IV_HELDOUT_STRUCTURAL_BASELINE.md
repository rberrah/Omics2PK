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

## Interprétation et suite

Un résultat CI vert prouve la reproductibilité technique d'un comparateur structural *pooled* sur de vrais C(t) humains ; ce n'est **pas** la preuve d'une prédiction PK par des omiques, d'un PopPK individuel, ni d'une validation externe indépendante. Les comparateurs E3 requis ensuite sont (1) NLME / SAEM clinique à effets aléatoires et interoccasion, (2) covariables cliniques et PGx disponibles, (3) effets omiques uniquement chez les **mêmes** patients avec horodatage antérieur à l'horizon prédit, et (4) test externe/temps futur avec incertitude de prédiction et BLQ.

Commande :

```sh
python -m pip install -e ".[pk-structural]"
python -m unittest discover -s tests -p 'test_iv_pk_structural.py'
python -m unittest discover -s tests -p 'test_mavoglurant_structural_fit.py'
python -m scripts.benchmark_mavoglurant_iv_structural
```

Résultats quantitatifs à reporter **après** preuve de succès du workflow ; ne pas inventer de RMSLE avant l'exécution.
