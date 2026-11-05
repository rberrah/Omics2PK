# Omics2PK — choisir les omiques selon le mécanisme, le variant et **la méthode de mesure**

**Version expérimentale du 5 novembre 2026.** Prototype *pré-modélisation* intégré à la [PR #137](https://github.com/rberrah/Omics2PK/pull/137). Ce module prépare le **choix des techniques** et les hypothèses interprétables ; il **ne prédit ni activité enzymatique, ni CL, ni gain de MAE, ni dose**. La hiérarchie est une **priorité d'interprétabilité mécanistique**, pas un classement de performances statistiques.

## 1. Question scientifique

**La meilleure couche dépend de la nature du mécanisme :** une mutation qui modifie l'épissage, une mutation faux-sens et une variante de régulation ne modifient pas les relations ADN→ARN→protéine→activité→PK de la même manière.

Le processus physio-pharmacologique attendu pour le médicament est le **premier filtre**. Le deuxième filtre porte sur les conséquences moléculaires documentées du variant. Le troisième décrit la **capacité réelle de la plateforme analytique à révéler ce mécanisme**, dans un tissu adapté et au bon instant.

Il y a une différence fondamentale entre :

- **« Intéressant à mesurer biologiquement »** : déduit de connaissances **sourcées** et de l'identité de l'assay ; cela produit une hypothèse expérimentale et des avertissements.
- **« Améliore effectivement la prédiction »** : nécessite des valeurs omiques réellement mesurées, un véritable critère PK, des **patients identiques**, un modèle clinique+PGx/PopPK, une sélection intra-entraînement et des patients de test non utilisés pour décider du panel.

La seconde affirmation **ne peut jamais** être déduite du seul mécanisme, d'un score SHAP, d'un réseau STRING, d'une mutation ou d'une modalité avancée.

## 2. La chaîne physio-pharmaco contrôlée

**Médicament et endpoint → mécanisme pharmacologique → organe → variant authentique → conséquence sur le transcrit → méthode omique → signal observable → comparateur PK indépendant.**

| Situation observée ou hypothèse | Attention biologique | Modalités à envisager pour tester la fonction |
|---|---|---|
| Mutation **régulatrice** authentique | L'expression ARN peut changer, sans changement de structure protéique | ARN quantifié et protéine dans le bon tissu ; activité si disponible |
| Mutation **faux-sens** | Quantité d'ARN et de protéine potentiellement normales malgré fonction altérée | Protéine variant-spécifique et surtout **activité enzymatique ou transport mesurée** |
| **Codon stop précoce** avec position du transcrit/jonction connue | NMD classique **possible** si stop ≥55 nt avant dernière jonction exon–exon ; exceptions nombreuses | ARN d'épissage / isoforme longue + peptide C-terminal ou protéforme intacte + activité spécifique |
| Stop **dernier exon** ou à proximité de la dernière jonction | NMD classique non attendu, mais transcrit ou protéine **potentiellement tronqués** et non fonctionnels | ARN total n'établit pas l'activité ; vérifier la protéine intacte et l'activité |
| Variant **d'épissage** (*CYP3A5*3, *6) | Le transcrit aberrant peut introduire un stop ; sans structure complète, **NMD non concluant** | ARN jonction / long-read, puis peptides informatifs et activité CYP3A5-sélective |
| Frameshift (*CYP3A5*7) | Le codon stop traduit dépend de la séquence et de l'isoforme | Reconstituer la conséquence transcript/protéine ; ne jamais inventer l'activité |
| Protéine présente mais **localisation/PTM** problématique | Un peptide commun ou une quantité totale peut être détecté malgré une fonction altérée | Protéome ciblé / protéforme PTM qualifiée + mesure de transport/activité appropriée |
| Concentration d'un **métabolite** | Peut refléter un phénotype intégré (enzyme, inflammation, environnement), mais souvent **non spécifique** | MS/MS avec étalon/identité qualifiée, cinétique et, si pertinent, flux traceur ; jamais « activité d'une seule CYP » par défaut |

### Exemple du codon stop : règle canonique ≠ diagnostic biologique

Un stop situé **à plus de ~50–55 nucléotides en amont de la dernière jonction exon–exon** peut déclencher la dégradation de l'ARNm par NMD (*nonsense-mediated decay*). Le programme accepte des positions **dans un même transcrit versionné** et sépare :
- **≥55 nt** : NMD classique *plausible*, **pas certain** ;
- **51–54 nt** : **frontière indéterminée** ;
- **≤50 nt** : NMD classique **non prédit par cette règle**, pas preuve de stabilité ;
- **stop après la dernière jonction** : NMD classique EJC **non prédit**, mais protéine potentiellement tronquée ;
- **positions manquantes** : **UNKNOWN**, jamais un arrêt de fonction ou une dégradation de l'ARN calculés par supposition.

Même en présence de NMD, une mesure d'ARN total peut capter de l'expression résiduelle, un transcrit d'une autre isoforme ou l'allèle sauvage. Un protocole 3'-tag, une mesure short-read d'un exon et un long-read permettant d'observer l'isoforme complète **ne répondent pas à la même question**. De même, en protéomique bottom-up, un peptide N-terminal conservé peut provenir d'une **protéine tronquée** ; un peptide C-terminal adapté est potentiellement plus informatif sur l'intégrité, **sans prouver à lui seul l'activité fonctionnelle**.

**Important CYP3A5 :** *3 correspond au gain d'un site accepteur cryptique dans l'intron 3 et à un transcrit aberrant ; *6 perturbe l'épissage ; *7 crée un frameshift. Le pipeline n'invente pas une position de stop à partir de l'étiquette star-allele et refuse de déduire *1 d'un dépistage incomplet. Les règles et les exceptions de NMD restent propres au transcrit et au contexte cellulaire.

Sources primaires : [consensus AMP/CAP/CPIC CYP3A4/5](https://pmc.ncbi.nlm.nih.gov/articles/PMC10565868/) ; [PharmVar GeneFocus CYP3A5](https://doi.org/10.1002/cpt.2563) ; [revue NMD EJC/50–55 nt](https://pmc.ncbi.nlm.nih.gov/articles/PMC10045457/).

## 3. Les méthodes ne sont pas interchangeables

| Couche | Méthode explicitement distinguée | Information possible | Limite incontournable |
|---|---|---|---|
| ARN | Gene bulk, 3'-tag, exon short-read | Quantité globale, présence locale de transcrit | Ne certifie ni épissage complet, ni longueur, ni fonction |
| ARN isoforme | Junction reads, long-read isoforme | Épissage anormal, usage des exons, isoforme | Une isoforme intacte n'est pas une protéine active |
| Protéine | Bottom-up peptides partagés | Abondance de fragments peptidiques | Peptides communs, paralogues et fragments tronqués |
| Protéine ciblée | Peptide N-terminal, C-terminal, variant-spécifique | Indices de région conservée ou manquante | Couverture de séquence et origine du peptide à vérifier |
| Protéforme | Top-down avec véritable identification | Intégrité, isoforme, PTM si mesurée | Protéforme intacte ≠ transport ou catalyse mesurés |
| Métabolites | MS1 global | Signaux potentiellement proches du phénotype global | Pic non identifié ≠ métabolite connu ou enzyme |
| Métabolites | MS/MS ciblé, référence authentique, isotope/flux | Identité chimique, abondance ou flux | Substrats multiples, flux non automatiquement isoforme-spécifique |
| Fonction | Sonde CYP sélective, partagée, test OCT2/MATE de transport | Activité sous conditions expérimentales explicites | Dépend du substrat, du tissu, de la méthode ; in vitro ≠ PK humaine |

Une sonde de midazolam mesure l'activité CYP3A **partagée**, pas une activité CYP3A5 isolée. Le prototype contrôle le code et le DOI de la sonde **T-5 N-oxydation** plus sélective versus **midazolam** selon le registre `cyp3a_activity_phenotype.py`. Une sonde fonctionnelle mesurée **après exposition au médicament** ne peut pas être réutilisée comme « omique prédose » sans risque de **fuite du critère**.

## 4. Le tissu et le paramètre PK priment sur la popularité d'une voie

**Tacrolimus.** Métabolisme CYP3A4/5 **foie et intestin**, efflux ABCB1 **intestinal**, effets sur CL et F, partition sang total dépendant de l'hématocrite. Une mesure CYP3A5 dans le sang périphérique est au mieux une **mesure substitutive** de l'activité hépatique/intestinale ; elle ne peut pas être déclarée équivalente. Une concentration C0 seule ne permet pas d'estimer séparément CL et F sans hypothèses.

**Metformine.** Clairance rénale impliquant OCT2/SLC22A2 et MATE1/MATE2-K (SLC47A1/SLC47A2) dans le rein. Un ARN OCT1 hépatique ne représente **pas** une mesure directe de CLr. Un signal P-gp/ABCB1 n'est pas automatiquement un mécanisme prioritaire de clairance de la metformine.

Le module actuel utilise exclusivement le registre de mécanismes **tacrolimus** et **metformine** déjà sourcé dans `physiopharmaco_cascade.KNOWN`. Les autres médicaments **s'abstiennent** tant qu'une carte organe–substrat–mécanisme n'a pas été curatée.

## 5. Intégration logicielle réalisée

**API** : `src/omics2pk/assay_aware_physio_triage.py::plan_physio_omic_acquisition(drug, endpoint, variants, assays)`.

Entrées : médicament/endpoint explicites ; conséquences de variants source déclarée avec statut de génotypage, identifiant, éventuel star-allele, transcript versionné et positions si connues ; **inventaire des assays** proposés ou réellement mesurés avec technologie, tissu, phase prédose, cible, peptide (région) ou identité chimique et référence.

Sorties **sans données patient** : pour chaque assay, processus PK pertinent, niveau d'interprétabilité, pertinence de l'organe, conséquences NMD interprétables ou inconnues, limites explicites, statut « candidat prospectif » / « candidat seulement pour test » / « abstention ». Aucun classement selon un prétendu **gain attendu de MAE**. Le plan indique les covariables PK indispensables, les bras **zéro omique / mono / multiomique**, et demande un benchmark indépendant.

**Démonstration sans patient** : `python scripts/demo_physio_assay_triage.py`. Deux situations :
1. *CYP3A5*3/*3* et tacrolimus : ARN 3' versus isoforme, peptides N/C, top-down non confirmé, activité T-5 sélective versus sonde CYP3A partagée, biomarqueur plasmatique CYP3A4/5 non spécifique ;
2. métformine–CLr : OCT2/MATE rénaux versus ARN OCT2 sanguin indirect, et OCT1 hépatique **rejeté comme explication directe de CLr**.

Le démonstrateur n'est **pas une cohorte clinique source**, et ne contient ni simulation d'un patient ni activité enzymatique ou concentration inventées. Les références de méthodes sont des **assertions de métadonnées** : une URL seule n'authentifie pas des fichiers de spectres originaux ou une vraie compatibilité échantillon–patient.

**Tests** : `tests/test_assay_aware_physio_triage.py` ; CI dédiée `.github/workflows/physiopharmaco-assay-aware-omic-triage.yml`.

## 6. Comment l'utiliser dans Omics2PK / la thèse ?

**Étape A — connaissance physiopharmaco** : carte contrôlée médicaments → enzymes/transporteurs → organe → paramètre PK, avec référence de preuve et absence d'attribution si le mécanisme manque.

**Étape B — conséquence du variant** : normal/régulateur/faux-sens/épissage/stop/frameshift, zygosité si réellement connue, transcrit de référence, niveau de certitude de la NMD. **Aucun coefficient quantitatif d'activité présumé.**

**Étape C — ce que mesure chaque assay** : gene-count vs splicing, C-terminal peptide vs protéine tronquée, protéforme intacte vs activité, métabolite d'identité incertaine vs métabolite spécifique. Vérifier organisme, organe, temporalité, méthode, référence/étalon, identité moléculaire.

**Étape D — design indépendant** : choisir a priori un petit nombre de voies ou d'assays **hypothétiques**, et laisser la validation croisée interne déterminer les signatures utiles parmi les couches réellement mesurées. Conserver explicitement `k=0` et les contrôles négatifs ; confronter les modèles aux vraies valeurs PK sur les **mêmes patients** au-delà du clinique+PGx complet+PopPK ; mesurer MAE/RMSE, erreur multiplicative, biais, calibration, pertinence de dose et incertitude selon l'endpoint. Déclarer les résultats nuls.

**Étape E — réconciliation mécanisme ↔ prédiction** : les attributs SHAP ou les associations STRING ne deviennent pas des mécanismes. Un résultat prédictif décalé du prior biologique crée une **hypothèse à tester** avec mesures de fonction, pas un gène causal affirmé.

### Ce qui reste non démontré

Ce composant n'a **pas** été validé sur les mêmes patients ayant génotype germinal CYP3A5 complet, RNA/protéines/métabolites, activité isoforme-spécifique et PK tacrolimus. L'exécution des tests logiciels ne prouve pas que la modalité priorisée prédit mieux. Une **cohorte humaine indépendamment mesurée** et des tests d'interaction mutation × technique seront indispensables.
