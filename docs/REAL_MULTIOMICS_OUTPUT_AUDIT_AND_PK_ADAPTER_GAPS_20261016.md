# Audit réel des sorties MultiOmics → Omics2PK (16 octobre 2026)

**Sources relues dans `rberrah/MultiOmics@33f6853c4d852ac67f5c7dafa51e7213d2dc3bd6`** : `src/lib/multiomics/downstream-exchange.js`, `src/lib/multiomics/research-export-package.js`, `docs/multiomics/DOWNSTREAM_EXCHANGE_V1.*`, `static/multiomics/demo_metadata.csv`, tests natifs de l'export ZIP et backend R longitudinal. Les modifications de ce PR concernent uniquement **Omics2PK** ; aucune correction n'est écrite dans MultiOmics sans coordination avec son agent.

## 1. Sortie effective actuelle

Le dossier de recherche `multiomics_dossier_recherche_v1.zip` contient, selon le type d'analyse :
- `exchange/manifest.json`, `exchange/samples.csv`, `exchange/features.csv`, `exchange/matrices/{transcriptomics,proteomics,metabolomics}.csv` (uniquement modalités présentes) : **valeurs et identifiants par assay**, pas seulement les gènes significatifs. Format `DOWNSTREAM_EXCHANGE_V1`; noms des variables conservés, échelles et namespaces déclarés, SHA-256 ; sans méthodes SWEET ni PK exécutées.
- `tableaux/{omic}.csv` : tableaux exhaustifs des résultats de chaque *analyse exploratoire* (effets, q, etc.) ou **source R native** uniquement si l'exécution R et la complétude ont été qualifiées. **Effet de groupe/q-value ≠ covariable biologique par patient**.
- `reference_r/*.csv` selon les méthodes complètes réellement disponibles ; `provenance/analysis_complete.json`, `provenance/manifest.json`, `methods/methods_audit.json` : traçabilité de l'analyse originale. Certaines tables R dépassant les limites de retour restent archivées localement, **hors ZIP**. La reproduction de résultats R dépend de l'environnement R/Bioconductor original ; JS exploratoire ≠ résultat de référence R.
- `source/projet_reproductible.json`, `reproduction/` conditionnel, `figures_browser/` et `figures_r/` : données sources et rapports/replays, mais **aucun score individuel RNA/SWEET qualifié et figé explicitement pour validation PK prospective dans le schéma V1**.

L'export peut être `downstream_exchange.status: not_available` si le CSV source n'est pas au format original attendu. Il faut regarder le **manifest externe** et non supposer que tout ZIP contient `exchange/`.

## 2. Test d'intégration du VRAI producteur JS

Le workflow Omics2PK a repris le producteur source MultiOmics actuel à SHA figé, généré de vrais fichiers V1 à partir des **démos synthétiques fournies par MultiOmics**, puis passé `inspect_exchange` et le nouvel `inspect_pk_readiness`.

| Observation directement vérifiée | Valeur |
|---|---:|
| Sujets biologiques distincts | **8** |
| Prélèvements distincts | **16** |
| Analyses / colonnes d'omics | **49** |
| RNA | **17** analyses |
| Protéomique | **16** analyses |
| Métabolomique | **16** analyses |
| Analyses RNA techniques supplémentaires (non nouveaux prélèvements) | **1** |
| Temps numériques utilisables dans le V1 exporté | **0 / 49** |
| Originaux doses / C(t) / données PK dans cet échange | **0** |
| Scores par patient directement réutilisables en covariables PopPK | **Non garantis par V1** |

Le CSV de départ décrit `subject_id,sample_id,assay_id,omic,condition,timepoint,batch,technical_replicate,outcome`, avec deux visites **`T0` et `T12`** chez huit sujets. L'exporteur V1 ne conserve que `collection_time` et `collection_time_unit`, qui étaient **absentes** des fichiers synthétiques d'origine, et ignore `timepoint`. **Les 49 assays perdent donc l'heure/date de visite dans l'interface actuelle.** La sémantique `T12` est inconnue : ne jamais deviner 12 heures, 12 jours ou 12 semaines.

Autres problèmes vérifiés :
- `technical_replicate` (source) n'est pas transformé en `technical_replicate_id` (export) : l'existence de R015A/R015B est déductible des identifiants mais le type de répétition est perdu dans la colonne normalisée.
- L'exporteur fixe `software.commit` à **`unknown`** : la chaîne de provenance source n'est pas figée dans le manifest V1 lui-même.
- La variable d'origine `outcome` n'est pas exportée comme endpoint individualisé ; c'est raisonnable pour un outil omique de ne PAS y confondre endpoints PK, mais **une sidecar PK séparée doit être liée par ID/study/time avec les mêmes sujets**.
- Les paramètres `transcriptomicsValues` (ici raw_counts), `proteomicsValues/metabolomicsValues` parfois `source_values_unspecified`, et les identifiants annotés `unresolved_source_namespace` ne peuvent être utilisés directement comme caractéristiques normalisées de réseau/PopPK sans preprocessing documenté. La tentative d'alimenter SWEET directement avec `raw_counts` est correctement bloquée.

**Important :** l'outil R de MultiOmics sait déjà analyser une *tendance longitudinale omique* (`condition * time`, modèle mixte via `lmerTest`, BH par famille) avec protocoles et tests dédiés. La lacune observée est l'**interopérabilité temporelle des résultats individuels exportés**, pas nécessairement l'incapacité scientifique de MultiOmics à analyser des séries longitudinales.

## 2bis. Test indépendant du ZIP complet, tel que remis à l'utilisateur

Le code du véritable **packager ZIP JS de MultiOmics** a aussi été exécuté avec les fichiers de démonstration synthétiques ; `provenance/manifest.json` et ses empreintes SHA-256 ont été confrontés aux membres archivés, et `source/projet_reproductible.json` aux `exchange/samples.csv`.

**Résultat : les 49 temps de visite catégoriels T0/T12 demeurent dans les fichiers sources du ZIP**, avec les huit identités biologiques. En revanche, **0/49 sont reportés dans l'export d'échange, et aucun temps chiffré utilisable** n'est ajouté. Le ZIP complet n'est donc PAS irrémédiablement dépourvu de temporalité : cette information est seulement manquante **dans le contrat standardisé V1 utilisé pour l'interopérabilité**. On peut récupérer les libellés du projet original sans inventer une unité ; toute conversion h/j exige le protocole authentique.

Ce test d'intégration examine le **générateur de dossiers original sur une simulation technique** et non un vrai résultat clinique ou un parcours utilisateur Playwright connecté à R. Il ne prouve aucun lien Omics–C(t).

Preuve CI : https://github.com/rberrah/Omics2PK/actions/runs/38954285830.

## 3. Modifications minimales demandées à l'agent MultiOmics

### Priorité P0 — conserver la vraie chronologie, sans attribuer d'unités fictives

Dans `exchange/samples.csv`, étendre de manière additive le format V1 ou créer V1.1 documenté :
- `timepoint` **tel que dans le fichier original** (ex. T0/T12), obligatoire si la colonne originale existe ; ne pas transformer arbitrairement ;
- `collection_time`, `collection_time_unit` seulement si valeurs numériques et unité définies réellement par le protocole ;
- `time_origin` / `reference_event` (ex. `first_dose`, `enrollment`, `transplantation`) et `time_axis_scope` (global, sujet, occasion), explicitement indisponibles sinon ;
- `tissue`, `specimen_type` (sang total/plasma/biopsie/etc.), `matrix` si distinct ; `batch`, `technical_replicate_id`, `biological_replicate_id`, `source_dataset_id` ;
- `drug_exposure_window` et PK temporality **ne doivent pas être inventés** par un export omique.

Tests à ajouter : entrée T0/T12 donne sortie T0/T12, `collection_time` reste vide tant que l'unité ne peut être prouvée ; vrai h/j entre correctement et résultat conserve le temps; 1 répétition technique ne génère pas un nouveau sujet.

### Priorité P1 — scores biologiques individuels, non statistiques de groupe

Exporter **optionnellement**, après une méthode qualifiée, un `exchange/individual_signatures.csv` : `subject_id,sample_id,assay_id,omic,signature_id,score,score_scale,fit_scope,reference_model_sha256,feature_definition_sha256,source_exchange_analysis_id`.

Il peut contenir une expression ou un score CYP/pathway préspécifié, un facteur MOFA ou une représentation réseau SWEET **uniquement s'ils sont véritablement calculés par échantillon**. Un `log2FC`, une `p-value` de groupe ou un réseau SWEET transductif ajusté aussi sur les sujets test **ne sont pas des scores individuels inductifs admissibles**.

Omics2PK a déjà un validateur optionnel `qualify_individual_signatures` : les projections des sujets test doivent utiliser **un modèle de référence figé sur les seuls sujets d'entraînement** et documenter le fit scope ; pas d'extrapolation silencieuse, pas de source SHA absent.

### Priorité P2 — provenance et accès

Mettre le **vrai SHA Git** producteur, identifier l'échelle/normalisation par modalité, préserver le lien aux fichiers d'annotation/namespace/gene mapping, contexte tissulaire et version de R/python, puis tests ZIP browser → téléchargement → `Omics2PK inspect_exchange`.

Conserver les sorties RNA seules, deux omiques ou trois omiques : aucune obligation multi-omique pour prédire la PK.

## 4. Répartition des responsabilités recommandée

**MultiOmics** calcule, qualifie et exporte les mesures / représentations biologiques individuelles avec identités et dates. **Omics2PK** accepte cet export et, **séparément et localement**, reçoit les doses/concentrations/paramètres NCA observés de la même étude. Le second outil vérifie unités, résolutions temporelles, appariement exact et données d'entraînement/test. **Pilier B** estime les effets sur CL/F, V/F, Ka, F ou paramètres structurels selon identifiabilité et données disponibles; la PopPK est ajustée sur les vrais C(t), pas sur les effets de groupe q-values. **PBPK** ne redevient pertinent que pour des voies biologiques et organes adaptés, avec transport physiologique et incertitude de translation. **Pilier C** valide ou refuse les mécanismes. **NewPK** peut rester un comparateur exploratoire plus tard, après référence standard solide.

## 5. Critères de réussite par stade

- **E2** moléculaire→fonction : preuve actuelle transcript ciblé→activité CYP (100 donneurs) et protéines→CLint HLM (21 donneurs), puis test futur sur vrais paramètres **NCA** avec Cmax, Tmax, AUC0–t **sans confondre durée de suivi entre patients**.
- **E3** Omics→PopPK : dose/time/C(t) réel, protocole patient-disjoint, référence clinique+PGx versus covariables RNA seules puis scores bio qualifiés. Effet sur CL/F, V/F, absorption, etc. évalué **en prédiction de concentrations et AUC sur sujets externes**, pas seulement via BIC de cohorte.
- **E4** omiques **elles-mêmes longitudinales** et PK longitudinales : variation état biologique aux landmarks et paramètres PK dépendants du temps, sans RNA obtenu après le début d'une fenêtre cible ni transport futur dans le passé.
- Aucune affirmation d'amélioration prédictive de **vraie PK clinique par MultiOmics** avant d'avoir un jeu source réellement apparié. Tests actuels V1 natifs = synthétiques et disponibilité logicielle seulement.

Les prochains changements du site `rberrah.github.io` dépendront des preuves et de l'accord des agents MultiOmics et PMx Explain, sous conventions FR/EN, sombre/clair, Resources et Other Disciplines.
