# Rosuvastatine n40 — audit officiel des 301 fichiers de données source MetaboLights MTBLS3128 (17 octobre 2026)

**Source scientifique** : Silveira et al. 2021 DOI [10.3389/fphar.2021.752960](https://doi.org/10.3389/fphar.2021.752960), n40 vrais volontaires, prédose métabolomique et AUC/Cmax de rosuvastatine mesurées.

**Question spécifique** : la page de l'utilisateur `https://www.ebi.ac.uk/metabolights/editor/MTBLS3128/files` contient-elle des matrices traitées qui permettraient de se passer des cinq liens Google Drive des auteurs ?

## Accès réellement testés et résultats

**Run public GitHub Actions** [39030060547](https://github.com/rberrah/Omics2PK/actions/runs/39030060547) ; scripts `scripts/audit_rosuvastatin40_complete_original_filelist.py` et `scripts/audit_rosuvastatin40_public_ftp_census.py`.

1. Endpoint officiel `https://www.ebi.ac.uk/metabolights/ws/studies/MTBLS3128/files?include_raw_data=true` (SHA-256 manifest `4faf93eab9887b23de56074b2f737703b5d31406bf93f1bebd41b77bf6570014`) : **8 entrées**, c.-à-d. **6 fichiers de métadonnées** ISA/MAF et **2 dossiers**. `include_raw_data=false` ne retourne que 6 objets. **Ce n'est PAS la liste récursive des fichiers de spectrométrie** ; affirmer qu'il n'y a que huit fichiers au total était erroné.
2. Guide EBI officiel [Files – MetaboLights](https://ebi-metabolights.github.io/guides/Files/) : données instrumentales brutes et dérivées sont dans le dossier `FILES` du FTP public, et non automatiquement téléchargeables par cette API de métadonnées.
3. Connexion **FTP anonyme réellement réussie** à `ftp.ebi.ac.uk`, dossier `/pub/databases/metabolights/studies/public/MTBLS3128/FILES`. L'inventaire récursif borné a visité **3 dossiers** : **301 objets fichiers** visibles, dont **145 fichiers `.mzML` et 156 archives `.zip`**.
4. **Aucun fichier `.csv`, `.tsv`, `.xlsx`, `.xls`, `.parquet`, `.rds` ou `.h5` visible directement dans cette arborescence FTP**, en dehors des six fichiers ISA/MAF au niveau racine. Les **156 ZIP** n'ont pas été téléchargés/ouverts ; aucune preuve n'autorise d'affirmer qu'ils contiennent ou qu'ils ne contiennent pas des matrices traitées. Toutefois aucun nom d'archive ne suggère explicitement `processed`, `quant`, `intens`, `AUC`, `Cmax` ou `clinical` dans l'inventaire restreint examiné.
5. [OmicsDI MTBLS3128](https://www.omicsdi.org/dataset/metabolights_dataset/MTBLS3128) annonçait **environ 1833 fichiers** au moment de son indexation. **Ne pas confondre ce total historique/tiers avec le nombre des objets directement visibles dans le dossier FTP à la date de l'audit** : il peut y avoir des objets internes d'archives ZIP ou une différence d'indexation. Les 301 sont précisément ce que la source officielle FTP actuelle a exposé, non une preuve d'exhaustivité de tous les spectres au monde.
6. **Aucune intensité ni concentration originale n'a été téléchargée ou publiée dans la CI**. Aucun nom de fichier exposant un échantillon/volontaire n'a été journalisé ; seuls nombres d'objets, types de fichiers et statuts source ont été sortis.

## Conclusion analytique

- **Métabolomique brute source vérifiée :** oui, FTP `mzML` + archives ZIP, et identifiants ISA n142 / assays 69+73 déjà inspectés par PR #49. Cette preuve ne transforme pas 142 acquisitions en 142 personnes (article : n40).
- **Table métabolomique traitée directement récupérable depuis fichiers source publics :** non retrouvée sous une extension tabulaire visible dans le FTP ; **aucune lecture des ZIP originaux**.
- **Cinq matrices traitées référencées par les auteurs :** original notebook GitHub `GustavoHBDuarte/Metabolomics_PK_proj` avec cinq URL Drive, direct HTTPS et `gdown` indisponibles depuis CI (PRs #52/#53).
- **Table PK individuelle vraie et jointe à la matrice prédose :** non qualifiée. Les AUC/Cmax de n40 sont des observations réelles de l'étude mais leur fichier originel et les identifiants d'appariement doivent être acquis, pas reconstruits depuis les graphiques.
- **Statut E3 :** `ORIGINAL_RAW_OMICS_FTP_VERIFIED_301_VISIBLE_OBJECTS; NO_DIRECT_TABULAR_PROCESSED_CSV; AUTHOR_FIVE_CSV_LINKS_UNAVAILABLE; PK_X_OMICS_CROSSWALK_UNVERIFIED; E3_ABSTAIN`.

**Meilleure voie opérationnelle :** le téléchargement humain via Google Drive des cinq CSV traités référencés dans le notebook des auteurs, puis le validateur privé ajouté dans [PR #54](https://github.com/rberrah/Omics2PK/pull/54). Le recalcule complet des 145 `mzML` et/ou des 156 ZIP n'est à envisager qu'en cas d'accès impossible aux matrices traitées ET après acquisition des PK individuelles véritables ; sinon la reconstitution du métabolome ne permettrait toujours pas la validation clinique.

