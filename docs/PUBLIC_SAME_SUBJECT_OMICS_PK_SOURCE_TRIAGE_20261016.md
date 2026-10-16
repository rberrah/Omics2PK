# Recherche des jeux publics liant réellement omiques et PK — 16 octobre 2026

**Périmètre :** sources originales d'articles et entrepôts GEO/SRA, Zenodo, Hugging Face, PK-DB. Les sources ne deviennent des preuves Omics→PK qu'avec **une correspondance biologique exacte**, la date de mesure moléculaire, la dose et un endpoint PK réellement mesuré. Aucun raccordement par âge/sexe/centre ou identifiants reconstruits après coup.

## Acquisition réelle dans le dépôt (CI et provenance)

### 1. Wang et al. 2026, transplantation hépatique / tacrolimus — prioritaire pour un premier endpoint PK clinique public

Article CC BY : https://doi.org/10.1371/journal.pone.0343817 (Wang G et al., PLOS ONE, 26 février 2026). Supplément clinique original **S1 Raw data XLSX** : https://journals.plos.org/plosone/article/file?id=10.1371/journal.pone.0343817.s001&type=supplementary . Données microbiome brute NCBI SRA : **PRJNA1416570**. Les métabolites fécaux sont présentés dans l'article, mais un **tableau individuel de métabolites portant des identifiants joints à S1** n'est pas présent dans l'export clinique S1 testé et n'a pas encore été vérifié ailleurs.

- Original **XLSX validé**, **15 275 octets**, empreinte **SHA-256 32c8a2d10c7c0448aebcb3ce12c32085098e24139b70c5b4c9fcd27ed3b37e49** ; 1 feuille et 37 lignes XML (16 champs originaux, 1 en-tête + **36 patients**).
- Colonnes originales : `No.`, âge, sexe, stade post-greffe, dose quotidienne mg, taille, poids, IMC, dose mg/kg, concentration sanguine de tacrolimus **ng/mL**, `C/D` en ng·kg/(mL·mg), antibiotiques, immunosuppresseurs et chirurgie. L'article mentionne des groupes d'exposition tacrolimus et microbiome. Le contrôle de la formule C/D sur les 36 lignes est **PASS**.
- **3 des 36 patients ont moins de 18 ans** : découverte à l'audit initial, qui avait correctement fait échouer une règle erronée supposant la cohorte exclusivement adulte. L'audit conserve désormais ce cas réel et rapporte le mélange d'âges ; comparaison pédiatrique dédiée **non permise avec n=3**.
- **Premier modèle de référence clinique seule sur les observations originales** : 36 plis leave-one-patient-out ; prédiction du C/D avec moyenne géométrique de l'entraînement (MAE **92,516919**) vs modèle ridge alpha=1 **âge, sexe, IMC** (MAE **92,392803**, amélioration **0,124116 unités C/D**, soit ~0,13%). **Ce résultat quasiment nul** est strictement un *contrôle de faisabilité de référence*, pas une preuve que ces covariables sont sans effet, ni une preuve de bénéfice omique. Faible n, mélange d'âges, covariables limitées, mesures observées une seule fois.
- **Important :** C0/D n'est ni AUC ni clairance/F sans conditions supplémentaires (horaire prélèvement, état d'équilibre, adhérence, schéma posologique) ; il ne s'agit PAS de NCA Cmax/Tmax ni de PopPK longitudinale. Aucun profil C(t) individuel ni PGx n'est établi dans ce supplément.
- **Goulot E3 :** démontrer la présence d'un champ source commun exact entre `No.` du supplément clinique et les 36 échantillons SRA et/ou métabolomiques, ainsi que les matrices individuelles quantifiées et les heures de collecte. L'accès à SRA ne prouve pas automatiquement ce lien. Tant que ce point n'est pas prouvé : **0 modèle microbiome/omics→C0 exécuté**.
- Exécution reproductible dans GitHub Actions : `python -m scripts.audit_public_tacrolimus_paired_omics_2026` puis `python -m scripts.validate_public_tacrolimus_36_c0_baseline` ; toutes les données originales passent **en mémoire uniquement**, aucun tableau clinique ni identifiant patient n'est commité/loggé. Les tests utilisent des valeurs explicitement synthétiques.

### 2. Fannin et al. 2010, paracétamol — candidat Omics longitudinales + PK clinique

Fannin RD et al., *Hepatology*, DOI https://doi.org/10.1002/hep.23330 ; original **GEO GSE15158**, https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15158 . Six volontaires ont reçu une dose unique de 4 g de paracétamol, trois ont reçu un placebo. RNA sang en prétraitement et 6/18/24/48/72/96 h après dose, selon le protocole. L'article rapporte une figure de concentration sanguine de paracétamol (échelle log10 µg/mL) durant 0–12 h ; *les valeurs numériques individuelles n'ont pas été retrouvées dans un tableau original vérifié*.

- L'archive officielle GEO Series Matrix compressée a été **réellement interrogée** ; les métadonnées de prélèvement, identifiants GEO et un tableau de valeurs d'expression traitées sont présents, contrairement au blocage GEO Seurat de l'étude Peng. Ceci ne démontre pas l'existence d'un tableau original des concentrations PK pour les mêmes sujets : **E4 RNA time-course oui, PK individuelle numérique non**.
- Un RNA mesuré à 18/24/48/72/96 h **ne peut PAS prédire rétrospectivement** l'AUC ou la Cmax de 0–12 h. Le signal après exposition peut soutenir l'axe PK→réponse biologique ou une prévision prospective ultérieure. Sous n=6 réellement exposés, toute estimation RNA→Cmax aurait une puissance extrêmement faible même si les concentrations brutes étaient retrouvées.
- Une digitisation de la figure papier = **donnée secondaire de faible fidélité**, non une vraie confirmation NCA ou une validation PopPK. L'original PDF complet fourni par l'utilisateur aiderait à vérifier les suppléments/identifiants et la présence éventuelle de tableaux numériques, mais **ne remplace pas le fichier PK source**.

### 3. Bushel et al. 2016, paracétamol oral répété — série RNA humaine

Article https://doi.org/10.1038/tpj.2016.8 ; GEO original **GSE70784** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE70784 , **432 échantillons d'expression** sur un protocole de 4 g/j pendant 7 jours comparé au placebo et suivi biologique jusqu'à 14 jours. GEO Series Matrix réelle accessible et vérifiée ; aucune série C(t) individuelle ni AUC publiée à fréquence comparable identifiée. **E4 longitudinal transcriptomique / DILI-PD**, pas E3 PopPK directe. Ne pas interpréter les 432 échantillons comme 432 sujets indépendants.

### 4. Degraeve et al. 2023, tacrolimus – microbiote – ABCB1 (souris) — candidat PBPK/PopPK mécanistique

*Microbiome* https://doi.org/10.1186/s40168-023-01578-y . NCBI GEO **GSE224034** transcriptome intestinal et SRA **PRJNA877868** microbiote. **351 vraies concentrations sanguines de tacrolimus**, environ 3,5 par souris, après administration orale répétée (3 mg/kg/j durant 5 jours), modèle populationnel Monolix SAEM avec absorption 1er ordre ; expression ABCB1 intestinale, variation microbiote et AUC. Des effets d'une déplétion microbiotique et de l'inhibition ABCB1 ont été mesurés.

**Important :** l'article ne fournit pas, dans les sources actuellement qualifiées, le fichier brut des 351 paires `(mouse_id,dose_t,C(t))` ni une correspondance numérique individuelle avec `GSE224034`. Les données RNA ne sont donc pas encore une matrice de covariables PopPK admissible. La publication est **très pertinente comme preuve mécanistique conceptuelle RNA/ABCB1→absorption→C(t)** et comme demande de données originales aux auteurs, mais ce n'est PAS un benchmark multicentrique humain. Contrôler en priorité les suppléments et la correspondance entre souris sacrifiées pour le transcriptome et souris suivies pour la PK.

### 5. Zhu et al. 2022, métabolomique plasmatique et tacrolimus après transplantation hépatique

*Frontiers in Pharmacology* DOI https://doi.org/10.3389/fphar.2022.1045843 . Cohorte longitudinale **137 transplantés / 693 prélèvements** : metabolomics non ciblée et **C0/D tacrolimus** répété, modèle mixte publié. Cette étude est **plus proche du couplage E4 clinique** que la cohorte PLOS 2026 de 36 sujets, mais l'**export individuel matrice LC-MS + C0/D + dose/time + patient ID** n'a pas été identifié comme librement téléchargeable dans le matériel public. L'article dit que les contributions sont dans les suppléments et invite à contacter les auteurs. Demander le fichier original au titulaire, pas dériver 693 lignes à partir d'une moyenne ou des coefficients publiés.

### 6. Zenodo — vraie biologie longitudinale sans PK

DOI de l'archive originale : https://zenodo.org/records/20699015 (Suchak et al. 2026), **24 personnes**, plasma métabolomique **18 prélèvements tous les 2 h**, suivi 34 h avec repas standardisés. Très utile pour tester une covariable métabolique **dépendante du temps** et déconfondre rythme circadien/repas, ou faire des tests de robustesse de `MultiOmics`; **aucune dose de médicament ni courbe C(t) du médicament** confirmée : ne jamais fusionner cette cohorte avec Theoph.

### 7. Autres ressources — triage négatif explicite

- **Open TG-GATEs** https://dbarchive.biosciencedbc.jp/data/open-tggates/LATEST/README_e.html : 170 substances, dose/temps de transcriptome rat et in-vitro, toxicité, biologie; dans les tables officielles README **aucun jeu de concentrations plasmatiques du médicament pour chaque rat** n'est indexé parmi les composants. Ne pas qualifier comme PK sans fichier TK primaire distinct. Intéressant pour modéliser dynamique toxico-génomique et PBPK *qualitatif*.
- **Hugging Face PMxBench** https://huggingface.co/datasets/pmxbench/PMxBench : benchmark de **question-réponse / raisonnement sur la PopPK**, non profils C(t) individuels associés à l'omique. **Hugging Face PKPD-Dataset** https://huggingface.co/datasets/Khalilbraham/PKPD-Dataset : corpus textuel d'articles, utile pour l'acquisition bibliographique, pas un jeu patient original.
- Autres jeux Hugging Face ou Github semblant comporter `PK_timeseries` et `gene_expression` doivent prouver l'origine primaire et le caractère **observé, non synthétique** avant toute évaluation scientifique. Un banc synthétique est utile pour CI et pédagogie, jamais pour n patients externes.
- **PK-DB beta** https://matthiaskoenig.github.io/pkdb/api/ : observations, interventions, C(t), NCA avec identités individuelles et groupes. Navigable sans compte ; **exports ZIP authentifiés** et chaque jeu peut contenir des courbes de groupe plutôt que des données individuelles. Très utile pour monter les comparateurs cliniques/PBPK et reproduire des articles, mais ne démontre pas une omique appariée si le même patient n'a pas été prélevé biologiquement.

## Priorités réellement exécutables sans cohorte privée

1. **P0 : Wang 2026** — source C0/dose humaine qualifiée, modèle référence clinique seule n=36 exécuté. Chercher les **IDs biologiques SRA PRJNA1416570** et le dépôt complet métabolomique de l'article. Ne lancer un modèle microbiome/omics→C0 que si source ID exacte + horodatage validés. Ce premier signal ne sera ni une AUC ni de la PopPK longitudinale.
2. **P0 : Fannin 2010** — GEO original RNA longitudinal accessible, manuscrit original + éventuelles données PK numériques complémentaires ou correspondance figure par sujet. Décrire prospective `RNA0`→Cmax/AUC0–12 vs `PK0–12`→`RNA24–96h` distinctement.
3. **P1 : Degraeve 2023** — vrai modèle tacrolimus oral mouse SAEM, transcriptome intestinal et microbiome ; obtenir original `ID×time×DV×dose` et mapping animaux RNA avant modèle mécanistique ABCB1 ou PBPK.
4. **P1 : Zhu 2022** — modèle clinique le plus pertinent : 137×693 C0/D-métabolites ; solliciter les données originales avec dates pour tester effets biologiques dépendants du temps en tenant compte du traitement et du changement post-greffe.
5. **En parallèle :** modèles PK cliniques de référence sur Theoph/mavoglurant déjà qualifiés, et validation des exports MultiOmics v1.1 lorsqu'ils seront livrés ; il est **interdit** d'utiliser les profils PK de Theoph et les métabolites Zenodo sur des sujets différents comme un jeu Omics-PK.

**État probatoire :** E2 *in vitro* déjà démontré ; PK clinique isolée (tacrolimus C0 et Theoph/mavoglurant C(t)) et omiques longitudinales isolées désormais qualifiées en source. **Aucun modèle prédictif humain Omics→Cmax/Tmax/AUC/CL n'est scientifiquement validé** sur des données individuelles appariées. Le projet doit conserver ce verdict même lorsque tous les jobs sont verts.
