# Omics2PK — PREUVE POSITIVE réelle : CYP2D6 + solanidine plasmatique, n=47

**Date : 2026-11-06. Statut : reproduction indépendante du modèle sur les données individuelles originales d'une étude publiée, mais aucune validation sur une population externe.**

## Source, provenance biologique et question

- Sarömba JA et al., *Frontiers in Pharmacology* 2026, DOI https://doi.org/10.3389/fphar.2026.1812007.
- Dépôt original des auteurs https://github.com/katjajust/Solandine-Metoprolol-Correlations à la révision **b6814967c0be7eae6f1843d0ec84e1987509bacc**, source unique `Solanidine_metoprolol_ratios_PGx.xlsx`; **SHA256 exact** `7397b13e43cd546b93f83a20b2dc84b6d49e97633173d56711c39973beedfe37`.
- Le classeur original contient **47 sujets appariés × 16 colonnes** sans absence sur les champs cliniques/genétiques/métaboliques requis. Aucun identifiant ni valeur individuelle n'est émis ou copié dans notre dépôt.
- Biologie : **plasma EDTA** prélevé chez les patients à jeun, **avant la dose suivante**, après au moins une semaine de métoprolol oral BID stable (95 ou 190 mg/j). Ce n'est **ni du tissu hépatique**, ni une analyse avant introduction du métoprolol.
- Génotype : deux allèles CYP2D6, copie CNV et activité prédite (**score AS**), estimées à partir d'une puce PGx + Stargazer. Le modèle utilise AS et CNV, plutôt qu'un codage indépendant de tous les diplotypes en n47.
- Biomarqueurs moléculaires **présélectionnés par le mécanisme CYP2D6** : **ln(4-OH-solanidine / solanidine)** (mesure LC-MS/MS ciblée) et **ln(SSDA / solanidine)** (métabolite de solanidine nommé, **semi-quantitatif sans étalon analytique authentique**). La solanidine est un alcaloïde d'origine végétale/alimentaire, **pas un stéroïde endogène**.
- Résultat biologique/cinétique indépendant : **concentration mesurée de métoprolol plasmatique ng/mL**, échantillon avant la prochaine prise. Ce n'est **pas une AUC**, ni une CL/F du patient, et l'horaire de la dose précédente peut être imparfaitement contrôlé.
- Mesure source `lnMR(OH-metoprolol/metoprolol)` expressément **exclue de toutes les covariables prédictives** parce que le dénominateur contient le médicament cible, créant une fuite par construction. La colonne `c(metoprolol)weighted`, dérivée de la concentration cible, également **exclue**. L'étude ne reconstruit pas la PK par une concentration drug-derived qui révèle l'issue.

## Notre replay sur les 47 vrais patients

Modèles préspécifiés à partir de l'étude : covariables cliniques **âge, poids, sexe, dose journalière de métoprolol, inhibition CYP2D6** ; ajout du **score CYP2D6 et de la CNV** ; ajout indépendant de chacun des deux rapports identifiés ou des deux ensemble.

- Même cohorte **47/47** sur tous les bras.
- Validation **LOPO** (un patient entièrement exclu) × 4 plis internes pour la **pénalisation Ridge** `alpha = {0.1, 1, 10, 100}`, sans sélection de caractéristiques par résultats de test.
- Entraînement sur le **log de la concentration observée**, optimisation sur **MAE en log-concentration interne** ; évaluation réelle de la MAE **ng/mL** et de la MAE log sur les sujets exclus, échelle de concentration directe.
- Standardisation et encodage des covariables catégorielles (sexe, inhibition, CNV) **fittés dans les plis d'entraînement uniquement**. Comparateurs source génétiques et métaboliques ne changent pas de population.
- Bootstrap apparié **1 999** tirages d'erreurs externes patients ; intervalles descriptifs **non corrigés pour la multiplicité** et **sans réestimation complète du modèle dans chaque bootstrap**. Les comparaisons ayant été motivées par le papier, ce n'est **pas une étude confirmatoire indépendante**.

### Résultat original vérifié

| Modèle | MAE concentration métoprolol (ng/mL) | MAE log | P20 | Gain de MAE versus clinique+PGx | IC bootstrap descriptif du ΔMAE (ng/mL) |
|---|---:|---:|---:|---:|---|
| Clinique seule | 51,644507 | 0,766727 | 10,6 % | — | — |
| Clinique + **CYP2D6 AS + CNV** | **47,086152** | 0,731443 | 14,9 % | référence | — |
| + **OH-solanidine / solanidine** | **41,022264** | 0,627696 | 21,3 % | **12,88 %** | [−14,553754 ; +1,956251] |
| + **SSDA / solanidine** | **36,289047** | 0,495339 | 27,7 % | **22,93 %** | [−21,699908 ; −1,187597] |
| + **deux rapports** | **34,873052** | 0,476897 | 29,8 % | **25,94 %** | [−25,033511 ; −2,131299] |
| Clinique + OH-solanidine sans PGx | 40,042033 | 0,623203 | 19,1 % | 14,96 % vs clinique+PGx | [−19,837386 ; +3,591875] |

**Exécution authentifiée** : [GitHub Actions #42179099625 — SUCCESS](https://github.com/rberrah/Omics2PK/actions/runs/42179099625). Script `scripts/benchmark_ariana_original_cyp2d6_solanidine_vs_pgx_20261106.py`. Original author-supplied source sheet + SHA audit [#42178678974 — SUCCESS](https://github.com/rberrah/Omics2PK/actions/runs/42178678974).

**Interprétation :** sur les 47 vrais patients, l'utilisation de l'un ou deux rapports enzymatiquement motivés **améliore l'erreur prédictive moyenne hors échantillon** comparativement à l'estimateur clinique + score génétique disponible ; il existe donc bien une **preuve de concept positive** pour Omics2PK, **sans pharmacocinétique complète, sans CYP3A ni tacrolimus**.

Le gain du duo a un intervalle d'erreur appariée négatif, mais **aucune significativité confirmatoire** ne peut être réclamée en raison des petits effectifs, comparaisons multiples et choix après connaissance d'une étude publiée. La MAE absolue reste élevée malgré l'amélioration ; ce n'est pas une recommandation posologique.

### Les variables ont-elles été choisies ou écartées et pourquoi ?

**Admissibles avant modèle :** le génotype CYP2D6 (AS+CNV), puis deux **substrat/produit alimentaires liés à CYP2D6**, disponibles dans le plasma au bon instant. Les biomarqueurs ont un mécanisme biologique nommé, indépendamment de l'amélioration MAE.

**Exclues explicitement :** métoprolol métabolite parent ratio (contient la variable réponse par construction), concentration pondérée du métoprolol (dérive de la cible), allèles bruts sur une forte dimension de dummies rares (capacité n47 ; les allèles restent disponibles pour futur modèle), biomarqueurs non identifiés ou postdose non contemporains. Pas de filtrage SHAP sur les données de test.

**Résultat du modèle :** le rapport OH-solanidine améliore numériquement l'erreur, mais sa CI descriptive croise zéro. Le SSDA améliore davantage dans ce replay, mais son identité structurale est moins certaine que son nom ne le suggère (absence de standard authentique). La combinaison est meilleure en MAE, sans garantie externe.

### Limites et expérience suivante

1. La **solanidine exogène** dépend de l'alimentation ; qualité et accessibilité d'un dosage sur les patients non consommateurs à tester.
2. Vérifier la **quantification du SSDA par étalon**, la précision du ratio, le lot chromatographique et les effets de l'inhibition (phénoconversion CYP2D6).
3. Évaluer un vrai **modèle PGx CYP2D6 diplotypes + activité résiduelle** à effectif suffisant et rechercher de la stabilité entre allèles.
4. Comparer avec un modèle PK **AUC ou CL/F indépendante**, avec posologie et temps contrôlés, et **validation externe** (autre hôpital, autre étude) ; cette cohorte est déjà utilisée pour développer l'hypothèse.
5. **Ne pas généraliser automatiquement de CYP2D6/métoprolol au CYP3A/tacrolimus.** C'est une justification de la méthode physiopharmacologique, pas l'endpoint final de la thèse.

**Message finançable** : des biomarqueurs plasmatiques fonctionnels, préspécifiés suivant une biotransformation CYP, peuvent capter une information qui échappe au score génétique et améliorer une prédiction réellement hors-sujet d'exposition mesurée. Le succès dépend de leur identité, du site/temps de prélèvement et de l'hétérogénéité médicamenteuse. L'indication cible et l'AUC restent à valider.
