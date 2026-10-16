# Étape E3 — modèle PopPK NLME pilote du mavoglurant

**Prototype méthodologique, non clinique, sans effet omique.** Cette étape vient après la comparaison 1C/2C IV perfusés à effets fixes du dépôt.

## Modèle

Deux compartiments IV perfusés, doses mg, temps h, volumes L et concentrations ng/mL. Pour le sujet biologique i : CL_i = CL_pop × exp(ETA_i), ETA_i suit une loi normale centrée de variance omega_CL². Le log de la concentration mesurée possède une erreur normale indépendante de variance sigma_log² conditionnellement à ETA_i. **Les deux occasions du même sujet partagent la même ETA**, sans transformer 196 visites en autant de personnes.

Les 6 paramètres positifs CL, Vc, Q, Vp, omega_CL et sigma_log sont ajustés par vraisemblance **marginalisée au niveau du sujet** via quadrature de Gauss-Hermite (9 nœuds par défaut), avec multidépart et limites physiologiquement larges. C'est un **prototype NLME** et non une implémentation SAEM, FOCEi ou un modèle IIV/IOV complet. Seule la clairance possède une IIV ; aucune IOV n'est estimée.

## Source, comparateur et séparation

Source officielle nlmixr2data::mavoglurant figée, RDA SHA256 : ba1d75e20bcc09f20a55a05a9e4424244efa15b2dfcef6079cd27385d5876174. Extraction inchangée : **120 personnes, 196 occasions retenues, 2 403 concentrations positives MDV=0** ; 53 MDV1 exclus et 2 occasions à temps dupliqués écartées. Le script refuse un changement du SHA ou des comptages. Les véritables fichiers individuels ne sont pas publiés dans GitHub.

Cinq plis déterministes **disjoints par personne**, incluant toutes ses occasions. Deux modèles estimés sur le train exclusivement : (a) modèle fixe 2C antérieur, (b) NLME 2C avec effet aléatoire individuel CL. Pour chaque personne tenue à l'écart, la prédiction NLME est **populationnelle a priori**, moyenne géométrique de C(t) intégrée sur la distribution théorique des ETA : **pas d'EBE, pas de mise à jour sur ses concentrations observées**. La métrique principale est la RMSLE à poids égal par personne. La log-vraisemblance prédictive **jointe par personne** est disponible pour le seul NLME, sans comparer cette densité à une densité résiduelle inexistante dans le modèle fixe.

## Résultat humain original, score tenu à l'écart des 120 personnes

Le [workflow source SHA-vérifié #38971073439](https://github.com/rberrah/Omics2PK/actions/runs/38971073439) a terminé avec succès : quatre tests synthétiques et un benchmark réel 120 humains / 196 occasions / 2 403 concentrations.

| Métrique CV cinq plis, poids égal par personne | 2C fixe | 2C NLME IIV(CL) |
|---|---:|---:|
| RMSLE log-concentrations populationnelles | **0,383383** | **0,384246** |
| Biais signé moyen log | −0,000180 | −0,010975 |

**Différence NLME − fixe = +0,000863 RMSLE**. Il n'y a donc **pas de gain de prédiction populationnelle ponctuelle démontré** avec la seule IIV(CL) marginalisée. La log-vraisemblance prédictive jointe du NLME est 5,17405 par humain avec signe négatif, **non directement comparable** au comparateur fixe sans modèle d'erreur probabiliste défini. Ni la comparaison 1C→2C, ni l'ajout de IIV ne valident une signature omique. Le contrôle de stabilité 7/9/15 nœuds est ajouté au workflow avant qualification finale.

## Ce qui reste à valider

Il ne faut **pas présumer** une RMSLE meilleure que le 2C fixe : un modèle de variabilité peut aider la vraisemblance/calibration sans améliorer la prédiction populationnelle ponctuelle. Les contrôles nécessaires incluent la stabilité 7/9/15 nœuds, estimations des erreurs/ETA, IOV, shrinkage, VPC, NPDE, BLQ, robustesse aux résidus et un **modèle SAEM/FOCEi indépendant** de référence. Le modèle actuel ne fournit pas d'incertitude paramétrique et ne convient pas directement à l'administration orale du tacrolimus.

**Aucune omique n'est présente dans cette source** : ni score SWEET, ni covariable multiomique, ni validation E3/E4 omique→PopPK n'est revendiquée. Les données individuelles dose/temps/C(t) avec biomarqueur pré-dose et identifiant commun demeurent le verrou pour Phapale, He ou Degraeve.

## Exécution

    python -m pip install -e ".[pk-structural]"
    python -m unittest discover -s tests -p test_mavoglurant_iv_nlme.py -v
    python -m scripts.benchmark_mavoglurant_iv_nlme

**Scores réels : à reporter uniquement après vérification du workflow CI sur source SHA authentifiée ; aucun chiffre humain NLME simulé.**
