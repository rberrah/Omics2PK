# E3 / PBPK — audit public tacrolimus intestinal 93 transplantés et APAP trois analytes

Mise à jour : 16 octobre 2026. **Aucune déclaration d'amélioration par omique** n'est autorisée sans données moléculaires *et* PK attribuables au même sujet avec les dates et l'origine des mesures.

## Cohorte Degraeve et al., CPT 2024, 93 transplantés rénaux

Article primaire : [Degraeve et al., *Clinical Pharmacology & Therapeutics*](https://doi.org/10.1002/cpt.3077), 115:104–115, [ClinicalTrials.gov NCT04360031](https://clinicaltrials.gov/study/NCT04360031). Dépôt 16S **NCBI BioProject PRJNA929904** : https://www.ncbi.nlm.nih.gov/bioproject/PRJNA929904 .

- **93 patients humains indépendants** stables, uniquement **CYP3A5 non-exprimeurs**. Un prélèvement fécal / patient, séquençage 16S avec 93 BioSamples selon le registre NCBI ; aucune série de transcriptomes humains associée par défaut.
- **C0 tacrolimus réellement dosée dans le sang** le jour de l'inclusion, doses quotidiennes originales, poids et hématocrite/âge. Sur une partie des patients, concentrations supplémentaires à C1 et/ou C3.
- **AUC0–24 non NCA complète** : calculée par un algorithme de limited sampling machine-learning publié par Woillard et collaborateurs. Les valeurs de prélèvements manquants ont été **remplacées par la médiane de la cohorte**. Cette AUC ne peut donc servir d'étalon or indépendant pour reproduire Cmax, Tmax, CL/F ou un PopPK sur profils denses.
- L'article montre une association taxonomique après ajustement sur hématocrite/âge ; la contribution microbiome issue de régressions **ajustées sur cohorte entière** n'est PAS une preuve prospective de performance sur individus tenus hors entraînement. Il faut requalifier sur le même groupe de patients et les mêmes identités originales.
- Supplement de l'éditeur : `cpt3077-sup-0001-Supinfo.zip`, environ 3,7 Mo, cité par Wiley. L'audit essaie seulement les liens **originaux HTTPS officiels** et conserve les statuts d'accès ; en cas de refus HTTP 403, ne pas contourner le service ni prétendre qu'une table clinique source est publique.
- Code source `scripts/audit_degraeve_93_kidney_tac_public_sources.py` interroge les métadonnées de **93 échantillons SRA**. Ni les reads bruts ni les identifiants individuels ne sont journalisés. **Même nombre 93 côté ENA et article ne prouve pas la correspondance sujet↔bibliothèque**, pas plus que 36 séquences SRA et 36 patients dans Wang2026.
- Il faut les champs originaux `patient_id ↔ run_accession/sample_alias ↔ C0/dose/kg ↔ AUC_estimated` et les *times* C0/C1/C3, l'information d'imputation et les ASV individuels avant tout modèle. Le NCBI PRJNA929904 ne dépose pas automatiquement les concentrations PK.
- Avec cette source vérifiée, le **premier objectif E3 raisonnable** serait un modèle *C0/(dose/kg)* comparant âge, hématocrite, PGx constant (CYP3A5 négatif) à des **taxa/ASV présélectionnés**. Au maximum, modèle cross-sectionnel observationnel, avec CV patients et tests par groupe/site ; aucune vraie PKpop dynamique si profils originaux incomplets.
- **Incertitude :** le microbiote peut être une cause de modification de biodisponibilité ou une conséquence des immunosuppresseurs / antibiotiques et du régime. Ne pas traduire corrélation d'abondance taxonomique en effet mécanistique sur CLint/F sans étude de sensibilité.

## Jeu public APAP parent + métabolites pour PBPK multi-analytes

R package `pksensi`, données `APAP.rda` de **nanhung/pksensi**, figées au commit `a775203eab90fa09059fde7190396d730ae04351`. Documentation [pksensi::APAP](https://search.r-project.org/CRAN/refmans/pksensi/html/APAP.html) et étude source de Zurlinden/Reisfeld (2016), modèles de paracétamol populationnels PBPK bayésiens.

- **32 lignes** publiées associant `Study`, `Wt` (poids moyen), `Dose` en mg/kg, `Time` en heures, et trois sorties concentration mg/L : paracétamol, **APAP-glucuronide** et **APAP-sulfate**.
- Il s'agit de profils provenant de publications anciennes et possiblement agrégés par étude ; **aucun champ d'identifiant individuel de personne dans l'objet R**. 32 lignes ≠ 32 volontaires indépendants et ce jeu ne permet pas une CV interindividuelle.
- Les deux métabolites mesurés sont des **analytes PK de produits de biotransformation**, **pas un profil métabolomique global** de 32 sujets. L'exploitation correcte est *test d'une implémentation PBPK et de voies métaboliques glucuronidation/sulfatation* sur profils rapportés et non une preuve "métabolome→PK".
- À terme, il fournit une référence pour qualifier un modèle multi-sortie `parent → formation AG/AS → distribution/élimination` avant d'ajouter des covariables UGT/SULT **uniquement si mesurées chez les mêmes personnes**.
- Script `scripts/audit_public_apap_three_analyte_pbpk.py` télécharge uniquement le RDA original, compte les lignes et observations réellement présentes, résume les unités et interdit de le qualifier comme PK individuelle.
- Aucun patient, aucune omique, aucune valeur originale par temps n'est exporté depuis le contrôle.

## Positionnement actuel

| Source | Identité patient | Médicament mesuré | Omique mesurée même cohorte | Données temporelles PK exploitables | Décision |
|---|---|---|---|---|---|
| Wang 2026 TAC foie | 36 lignes cliniques | 36 C0 et dose | 36 reads SRA distincts, appariement exact aux cliniques inconnu | Une mesure par patient | E3 **bloqué par crosswalk** |
| Degraeve 2024 TAC rein | 93 patients, 16S | C0 mesurée, AUC **estimée** (sparse) | 93 prélèvements 16S dans le protocole | C0 + parfois C1/C3 ; fichier source individuel non vérifié | **Priorité source/ID** |
| Fannin 2010 APAP humain | 9 volontaires, 88 analyses RNA | C(t) 0–12h sous forme de figure | RNA après traitement | Original numérique PK par patient non vérifié | **Priorité suppléments PDF** |
| Degraeve 2023 TAC souris | 99 animaux, 351 concentrations | PK vraie, modèle SAEM | Microbiote / ABCB1 / transcriptome | Original 351 ID×time×DV non librement vérifié | **Priorité PopPK préclinique** |
| pksensi APAP humain | **Pas d'identifiant patient** | Parent et 2 métabolites (3 C(t)) | Pas d'omique | Lignes d'études agrégées | **PBPK multi-analytes**, pas E3 Omics |
| Zhu 2022 TAC foie | 137 personnes, 693 prélèvements | C0/dose répétée | Métabolomique même protocole | Original apparié individuel non vérifié | **Meilleur candidat E4 clinique** |

**Critère de franchissement scientifique :** une URL accessible et le nombre de lignes ne suffisent pas. Il faut établir l'identité du participant entre assay, original dosage et matrice moléculaire, l'origine temporelle relative aux doses, les unités, la censure et la provenance. **Résultats des vérifications originales exécutées :** la requête officielle ENA confirme **93 16S runs et 93 BioSamples distincts** ; les deux URLs officielles du supplément Wiley renvoient **HTTP 403** à l'audit non interactif (le fichier n'a PAS été récupéré). Pour APAP, l'objet `APAP.rda` original (SHA256 `08b1d33b96d41cdd12fd48b2f09ad28c68728bf047ef4c93772a97af3dca00e0`) comporte **32 lignes sur 3 études** dont **30** ont le parent et ses deux métabolites quantifiés. CI d'origine [GitHub Actions](https://github.com/rberrah/Omics2PK/actions/runs/38960552795), tous les contrôles ont réussi, mais le succès du contrôle ENA n'équivaut pas à l'obtention des 93 valeurs PK individuelles.

Le contrôle des formats et fichiers **ne remplace pas** une validation biologique à l'aveugle.

La priorité est de récupérer ou faire confirmer les tables numériques originales et leurs concordances d'identité ; uniquement ensuite réaliser la comparaison clinique+PGx vs signatures MultiOmics sur E3/E4. Aucun couplage fictif des cohortes Theoph, pksensi et Zhu n'est autorisé.
