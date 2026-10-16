# Sources réelles E3/E4 : prédose métabolomique→PK tacrolimus et NCA longitudinale observée — 16 octobre 2026

## 1. Résultat E3 PK : mavoglurant n=120 humains, 198 occasions

Source originale **`nlmixr2data::mavoglurant`** (Wendling et al. 2016 / A2121 IV), commit figé `f2cfb01ade88d4d30e66f56b19efce2a2c6e26bd`, source `data/mavoglurant.rda` SHA256 `ba1d75e20bcc09f20a55a05a9e4424244efa15b2dfcef6079cd27385d5876174`. La source initiale comporte **2 678 lignes** : 198 doses IV et 2 480 lignes EVID0. **EVID0 ne signifie pas forcément concentration réellement quantifiée**, car certaines lignes ont `MDV=1`. Référence publique : https://nlmixr2.github.io/nlmixr2data/reference/mavoglurant.html ; la vignette officielle filtre les observations par `EVID==0,MDV==0,DV>0` : https://nlmixr2.github.io/nlmixr2lib/articles/Elmokadem_2023_mavoglurant_pbpk.html .

### Résultats réellement calculés en CI sur observations sources

| Quantité | Résultat calculé |
|---|---:|
| Sujets humains biologiques | **120** |
| Occasions originales `sujet×OCC` | **198** |
| Occasions avec une dose IV réelle / chacune | 198 |
| Sujets à deux occasions initialement | **78** |
| Lignes EVID0 candidates à une mesure | **2 480** |
| Lignes MDV1 explicitement exclues | **53** |
| Occasions avec temps dupliqués non résolus, exclues de NCA | **2** |
| Occasions véritablement admissibles à la NCA descriptive | **196** |
| Concentrations réellement utilisées, échelle ng/mL | **2 403** |
| Sujets avec **deux** occasions admissibles | **76** |
| Médiane prélèvements quantifiables / occasion admissible | **12** |
| Cmax **observée** médiane / occasion admissible | **556,5 ng/mL** |
| Tmax **observé** médian depuis administration | **0,20 h** |
| AUC du **premier au dernier temps effectivement observé** (médiane entre profils) | **1 096,619875 ng·h/mL** |
| Premier/dernier prélèvement médian | **0,20 h / 24,20 h** |

Calcul : `python -m scripts.validate_mavoglurant_120_real_occasion_nca`, source originale RDA chargée **en mémoire éphémère de CI**, sans export d'identifiants, de paramètres NCA par sujet ou de concentrations. Evidence : [GitHub Actions, source exacte](https://github.com/rberrah/Omics2PK/actions/runs/38965413511). Les tests synthétiques de contrats refusent les effets de dose inventés, les temps dupliqués, les mesures `MDV=1`, les durées de perfusion aberrantes.

**Limites scientifiques non négociables :** il ne s'agit pas de `AUC0_last` strict si la mesure t=0 n'existe pas (aucun profil ne possède de t0 quantifié) ; notre cible est **AUC(t_premier,t_dernier)**. La fenêtre réelle est variable et ne permet pas une comparaison brute des AUC entre occasions sans alignement. Cmax et Tmax sont **des maxima observés** sous un calendrier de prélèvement et non les paramètres d'une fonction continue estimée ; le Tmax réel peut précéder la première observation. **AUC∞ et CL IV non inférées**. Les 198 occasions ne font pas 198 patients indépendants ; 76 personnes ont deux occasions admissibles. Aucune omique dans ce jeu. Le prochain essai E3 sera un véritable modèle structurale PopPK/PBPK *à partir des données C(t) IV originales*, avec split par sujet (et ses deux occasions), variance inter-sujet/occasion et confrontation aux résidus. Ne pas alimenter un modèle avec l'AUC calculée à partir de ses données de validation.

## 2. Découverte prioritaire E3 humain : Phapale et al., tacrolimus, 2010

Publication *Clinical Pharmacology & Therapeutics* **DOI 10.1038/clpt.2009.296** : https://doi.org/10.1038/clpt.2009.296 .

- **29 volontaires sains originaux**, LC-MS sur **1 256 ions métabolites urinaires pré-dose**.
- Après administration **orale** de tacrolimus, **concentrations plasmatiques jusqu'à 72 h** sur les **mêmes volontaires** ; paramètres PK évalués et prédits par PLS à partir des métabolites pré-dose. C'est un protocole « état moléculaire antérieur → vraie PK future » où l'ordre temporel est biologiquement valide.
- L'éditeur signale **13 pièces supplémentaires d'origine `.doc`**, numéros `cptclpt2009296-sup-0001.doc` à `...-0013.doc` : https://ascpt.onlinelibrary.wiley.com/doi/abs/10.1038/clpt.2009.296 . **Le contenu des 13 documents, les colonnes individuelles 29×1 256, les concentrations 0–72 h et la correspondance d'identité n'ont PAS été qualifiés**.
- L'étude **n=29 et p=1256** demande sélection de variables entièrement nested et des intervalles d'incertitude honnêtes ; si 1 256 métabolites sont projetés in-sample en PLS puis évalués sur les mêmes sujets, on n'obtient pas une preuve de transportabilité. Benchmark envisagé : dose/PK clinique seul → urine métabolites pré-dose → scores biologiques qualifiés avec 29 folds LODO + validation indépendante possible ; comparer Cmax/AUCinf/Tmax selon les données réellement disponibles, pas extrapoler une AUC sans série d'origine.
- **Action prochaine** : inspecter les 13 suppléments Wiley ou demander aux auteurs la matrice urinaire quantifiée et la vraie dose/temps/PK du même volontaire. Le PDF seul ne suffira pas forcément ; ne pas transformer figures / PLS loadings en valeurs individuelles.

## 3. Autre candidat E4 réel : He et al., 2022 n=109 renal tacrolimus

*Individualized Dosage of Tacrolimus for Renal Transplantation Patients Based on Pharmacometabonomics*, DOI https://doi.org/10.3390/molecules27113517 , PMC9182099.

- **109 greffés rénaux**, métabolomique plasmatique prélevée au **J7, M1, M3**, avec prélèvements de **C0 tacrolimus avant la dose matinale** selon méthodes ; 2 695 caractéristiques LC-MS initiales, **1 514** après prétraitement, KNN et filtres QC.
- La publication vise d'abord la **réponse pharmacodynamique** des lymphocytes T, avec 73 réponses « normales », 16 réponses insuffisantes et 20 réponses excessives ; cela ne fournit **aucun modèle PopPK C(t)** ni les doses horaires. Cependant la **répétition des métabolites, l'orientation pré-dose et l'existence d'une C0** en font une **piste E4 forte si les valeurs individuelles et dates sont accessibles**.
- Supplement original **5,2 Mo** avec données MS brutes déclaré sur https://www.mdpi.com/article/10.3390/molecules27113517/s1 ; la tentative `PMC .../bin/molecules-27-03517-s001.zip` retourne **HTML**, pas l'archive. Ne pas conclure absence de données ; il faut récupérer la vraie pièce via éditeur ou l'utilisateur et auditer **identifiants sujet×D7/M1/M3 et C0 réelle**. Rien ne prouve encore que toutes les 327 visites (109×3) soient effectivement complètes.
- Le document contient des résultats concernant T-lymphocytes et dose moyenne, mais **C0/D individuelle n'est pas vérifiée dans le supplément**. Aucune affirmation modèle Omics→C0 jusqu'à source originale appariée.

## 4. Autre contrôle public : Burghelea et al., 2022 n=42 transplantés

DOI https://doi.org/10.3390/biomedicines10051157, `PMC9138871`. Metabolomique lipidique n42 **sélectionnés sur groupes TAC C0 extrêmes : 19 haut >8 ng/mL, 23 bas <5 ng/mL**, parmi 135 patients suivis. 336 composés LC-MS répertoriés ; supplément d'environ 249 Ko mentionné mais tentative de lien public PMC `.../bin/biomedicines-10-01157-s001.zip` renvoie **HTML** et non source ZIP. Ce plan est **cas-témoins dépendant du résultat C0**, pertinent pour classification exploratoire, MAIS **biais de sélection** si utilisé comme jeu de régression continue C0/AUC sur population générale, et pas de séries C(t) ni de véritable Cmax/Tmax. Échelle unité de TAC sang vs serum LC-MS à préserver distinctement.

## 5. Priorités pour Omics2PK

1. **E3 vrai PK sain** : Phapale n29 métabolome pré-dose → concentrations TAC 72 h, **originales numeric, IDs, doses, times indispensables**. Article+13 suppléments et demande d'accès aux sources.
2. **E4 109 greffés longitudinal** : He n109 D7/M1/M3 omiques + dosages C0 réels dans protocole, **source matrie quantifiée et triple visite identitaire à obtenir** ; validation patient-disjointe et chronologie.
3. **E4 n137/693** : Zhu 2022 plasma métabolomique → C0/dose répété ; fichier original individuel manquant.
4. **PopPK socle réel et numérique** : Mavoglurant n120×198 réel, 196 occasions NCA admissibles, 76 sujets à deux temps, **pas omique**. Construire et vérifier modèles 1C/2C et PBPK sur C(t) avec splits biologiques ; aucun prétendu gain omique.
5. **PBPK métabolites** : APAP parent/AG/AS 3 études / 32 observations source, seulement comparaison de courbes agrégées, pas données patients.

**Statut scientifique : E2 biologique→fonction in vitro validé ; NCA descriptive vraie PK humaine désormais calculée sur 196 occasions ; E3/E4 Omics→PK humain toujours **non démontré** sans données individuelles appariées.**
