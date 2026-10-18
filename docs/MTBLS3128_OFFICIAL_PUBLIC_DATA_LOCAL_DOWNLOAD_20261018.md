# Télécharger directement les sources publiques MTBLS3128 (18 octobre 2026)

**Étude originale :** Silveira et al., rosuvastatine 20 mg, pharmacométabolomique pré-dose chez 40 volontaires humains, DOI [10.3389/fphar.2021.752960](https://doi.org/10.3389/fphar.2021.752960), dépôt public [EBI MetaboLights MTBLS3128](https://www.ebi.ac.uk/metabolights/editor/MTBLS3128/overview).

Le téléchargement public *ne nécessite pas d'identifiants Inserm ou universitaires*. Les données brutes et ISA sont déposées sous le [FTP officiel](https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS3128/). Le navigateur de l'éditeur de l'EBI n'est pas une interface garantissant les matrices de concentrations individuelles. Aucune authentification EBI n'est détournée.

## Commandes depuis un ordinateur local ayant accès à l'EBI

Depuis la racine du dépôt Omics2PK avec Python >=3.11 :

\`\`\`bash
# 1. Vérifier le nombre de fichiers disponibles, sans les télécharger
python scripts/download_mtbls3128_original_public.py --list-only

# 2. Télécharger la totalité des métadonnées ISA originales (étape courte)
python scripts/download_mtbls3128_original_public.py

# 3. Essai contrôlé de deux vrais fichiers mzML, avec nouveau dossier privé
python scripts/download_mtbls3128_original_public.py \
  --output-dir private/MTBLS3128_mzml_pilot \
  --mzml-max-files 2

# 4. Si espace disque/réseau suffisant, récupérer l'ensemble des mzML
#    N doit être le nombre effectivement indiqué par --list-only.
#    Exemple si la source officielle liste 145 mzML :
python scripts/download_mtbls3128_original_public.py \
  --output-dir private/MTBLS3128_mzml_complete \
  --mzml-max-files 145 \
  --max-total-gib 200 \
  --max-file-gib 10
\`\`\`

**Règles :** le script ne mélange jamais les fichiers de différentes études, utilise la connexion FTP **anonyme officielle**, écrit sous \`private/\` (déjà ignoré par Git), applique des limites par fichier/globales, et garde chaque fichier source tel quel avec un SHA256 dans un manifeste **local uniquement** \`private_source_manifest.json\`. En cas de téléchargement interrompu, le fichier partiel est supprimé ; les copies existantes ne sont pas écrasées sans \`--overwrite\`. Les ZIP des acquisitions instrumentales ne sont *pas* téléchargés par défaut ; \`--zip-max-files N\` nécessite une sélection explicite et une capacité disque suffisante. Les fichiers sont sélectionnés par nom de source trié (pas par ID clinique). Aucun import/mappage aux 40 volontaires n'est inféré.

## Quelle conclusion scientifique ?

Le précédent recensement officiel source a identifié **6 fichiers de métadonnées ISA**, **145 mzML** et **156 ZIP** dans la hiérarchie courante observée de \`FILES\` au 17 octobre 2026 ; ces nombres sont un recensement ponctuel de 301 fichiers, **pas 301 individus** ni la totalité garantie de l'index historique OmicsDI. Deux acquisitions ionisation + / − et contrôles qualité doivent être distingués des 40 personnes d'origine. Les ISA décrivent les échantillons/assays, mais pas une table validée de concentrations de rosuvastatine ou AUC/Cmax individuelles.

Les \`mzML\` contiennent des spectres et chromatogrammes, **pas une matrice feature×volontaire immédiatement réutilisable** : il faut réaliser la détection de pics, alignement, normalisation QC/blank et regroupement des échantillons, en respectant le moment pré-dose et les unités/étalonnages. Le dépôt \`rberrah/MultiOmics\` fournit une interface pour des matrices analytiques, mais ne réalise pas automatiquement toutes les étapes de prétraitement LC-MS des spectres dans Omics2PK.

**Verrou E3 persistant :** les 5 CSV traités référencés dans le [notebook original des auteurs](https://github.com/GustavoHBDuarte/Metabolomics_PK_proj/blob/main/Metabolomics_Elastic_net.ipynb) n'ont pas été récupérés par l'accès automatique public précédent. Aucun fichier source officiel vérifié ne fournit encore ensemble les **40 AUC et Cmax quantitatives par volontaire avec un ID identique/une table de correspondance exacte vers les intensités pré-dose**. Il faut ces données pour utiliser le benchmark local préalablement préparé dans \`scripts/benchmark_private_original_rosuvastatin.py\`. Ne pas estimer AUC à partir d'une figure ou fabriquer des correspondances en triant les lignes.

## Confidentialité et validation

- Ne pas ajouter les mzML individuels, les CSV auteurs, les concentrations ni le manifeste détaillé dans GitHub, un workflow public ou les artefacts GitHub Actions.
- Les tests du téléchargeur sont **synthétiques**, et ne démontrent pas un nouveau téléchargement réel dans l'environnement de ChatGPT. La récupération effective requiert le réseau FTP EBI depuis l'ordinateur exécutant le script.
- Seuls les rapports non identifiants et agrégés peuvent être partagés ensuite ; validation de la licence des données et de leur provenance avant redistribution.
