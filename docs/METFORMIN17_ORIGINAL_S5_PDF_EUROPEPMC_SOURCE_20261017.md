# Metformine n17 — supplément original Table S5.1/S5.2 retrouvé officiellement (17 octobre 2026)

**Étude** : Tee et al., *Pharmaceutics* 2022, DOI [10.3390/pharmaceutics14061268](https://doi.org/10.3390/pharmaceutics14061268), PMC [PMC9231303](https://pmc.ncbi.nlm.nih.gov/articles/PMC9231303/), **licence de l'article CC BY**.

**Critère E3** : disposer simultanément des observations individuelles d'une omique **pré-dose** et des concentrations PK **quantifiées en unités plasmatiques** avec *strictement les mêmes identifiants de 17 volontaires*, dose 1000 mg, temps 0 à 24 h et condition de prélèvement. **Statut actuel : concentrations PK sources retrouvées en PDF, matrice individuelle brute de métabolomique non trouvée.**

## Source réelle originale — avancée de l'audit

Source officielle testée avec succès dans GitHub Actions : [Europe PMC supplementaryFiles, PMC9231303](https://www.ebi.ac.uk/europepmc/webservices/rest/PMC9231303/supplementaryFiles). L'API renvoie une **archive ZIP de paquet article**, de SHA256 `e7f8329c2e7262ae6b4b76e8f903a9aedd64486dc0e9151c1e03fb32648fa34b`, contenant un **ZIP supplémentaire original** `pharmaceutics-14-01268-s001.zip` (452 006 octets non compressés). Cette archive imbriquée ne contient qu'un **PDF** `pharmaceutics-1731469-supplementary.pdf` (522 715 octets), **17 pages** dont **15 avec une couche texte extractible**.

L'audit source réel [GitHub Actions 39022212659](https://github.com/rberrah/Omics2PK/actions/runs/39022212659) et le module `scripts/audit_metformin17_441k_original_supplement.py` ont vérifié **dans le PDF original** l'existence littérale des intitulés **Table S5.1** et **Table S5.2**. Les 17 × 15 concentrations plasmatiques quantifiées et le protocole sont décrits par les auteurs. Cette vérification prouve l'existence des tableaux d'origine, mais **n'établit pas encore la conversion fidèle de toutes les cellules d'une grille de 255 mesures** à partir du PDF ni l'intégralité/les éventuelles valeurs sous LLOQ : une extraction PDF supervisée par sujet et vérification des unités est requise avant modélisation.

**Important : aucun PDF, fichier d'auteur, identifiant volontaire, concentration individuelle, matrice m/z ni image de ses tableaux n'a été republié sur GitHub.** Le script ne sort que des hachages de l'archive originale, noms de membres, dimensions et présences de labels, pas les mesures individuelles.

### Pourquoi les premières recherches ne trouvaient pas le supplément

- Le JSON PMC OA `PMC9231303.1` déclare **zéro** objet supplémentaire Office XLSX/DOCX. Un supplément PDF zippé ne fait **pas** partie de cette sous-liste : ce résultat était techniquement vrai, mais ne justifiait pas l'absence de supplément original.
- L'URL directe PMC article/bin a retourné une erreur HTTP et le media S3 testé également, tandis que le endpoint officiel Europe PMC **a réussi**.
- Le paquet Europe PMC contient les images d'article **et le ZIP de supplément imbriqué** : il a fallu inspecter deux niveaux d'archive pour identifier les vrais fichiers.

### La matrice omique brute reste absente du ZIP

Le seul fichier dans le ZIP supplémentaire original est le PDF de 17 pages, contenant tableaux des données cliniques et du protocole. **Aucun fichier CSV/XLSX/TSV ou mzML/mzXML brut des intensités des métabolites individuels n'a été retrouvé dans cette archive.** La preuve dans l'article de `17 paired samples` en pré-dose/à Cmax et les figures PCA/PLS-DA ou voies de mummichog ne donnent pas les **valeurs d'intensités de chaque individu**. Ce qu'il faudrait : `subject_id × pre-dose metabolite features` en positif et négatif, normalisation, batch et identifiants communs aux tableaux S5.1/S5.2.

Sans cela : **NO_CLINICAL_OMICS_TO_PK_BENCHMARK** ; refuser de fabriquer des concentrations depuis les figures, les voies métaboliques ou les charges PCA.

## Décision scientifique

| Critère | État après audit |
|---|---|
| Étude chez 17 mêmes volontaires avec métabolomique | Oui, protocole primaire original |
| Dose et temps PK mesurés (1000 mg, 15 prélèvements/personne) | Oui, protocole |
| Vraies tables PK individuelles S5.1/S5.2 retrouvées | **Oui, PDF original 17 pages, deux labels identifiés** |
| Conversion structurée exacte 255 concentrations et BLQ | **Non démontrée**, PDF en attente de vérification de l'extraction |
| Originale matrice brute des features LC-MS de pré-dose | **Non trouvée dans le supplément PDF/ZIP** |
| Même identifiant volontaire PK ↔ omiques dans fichiers accessibles | **Non qualifié** |
| Omics→AUC/Cmax CL et amélioration clinique hors-échantillon | **Non**, aucune exécution possible honnête |

**Suites prioritaires :** (1) reconstruire en privé les 17×15 concentrations depuis les tableaux S5 avec vérification PDF et contrôle de LLOQ ; (2) rechercher la matrice individuelle LC-MS native dans d'éventuelles dépôts MetaboLights/Metabolomics Workbench/GNPS ou prise de contact avec les auteurs, (3) uniquement ensuite LODO à 17 donneurs, et évaluer des indicateurs moléculaires pré-dose fixés avant traitement.

