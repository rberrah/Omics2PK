# Banque fédérée « Spectres → Molécule → Fonction → PK » — prototype de recherche
**État : implémenté et testé en mode local, annotation exploratoire, pas d'identification nouvelle des 229 inconnus Rotroff. 2026-11-04.**

## Principe : ne pas inventer une bibliothèque de spectres
Le prototype `src/omics2pk/spectral_knowledge_bank.py` crée une petite base SQLite locale avec quatre tables **distinctes** :
1. `spectral_reference` — référence expérimentale à un ionogramme `(m/z, intensité)` et métadonnées instrumentales, structure candidate déjà attribuée par la source, accession originale et droits de réutilisation;
2. `mechanism_link` — assertion sourcée `CHEBI/InChIKey/HMDB → enzyme/transporteur, type de relation, organe, référence`. Une identité structurale candidate et une littérature sur un transporteur ne prouvent **jamais** une activité dans le patient.
3. `omic_feature` — identifiants originaux de signaux non identifiés et métadonnées de prélèvement, **sans inventer de masse, de composé ni de fonction**.
4. `protein_proteoform_reference` — référence originale de protéines/peptides/protéformes (accession PXD, UniProt, site PTM, niveau de preuve), **séparée des spectra GC/LC de petites molécules**. Le moteur spectral protéomique complet n'est pas encore implémenté.

Le catalogue fédéré décrit **où chercher** et les types de données compatibles. Il **n'aspire pas** automatiquement de grandes bases tierces, ne redistribue pas les spectres soumis à licence et ne présente pas un résultat de recherche comme vérifié sans l'avoir exécuté.

### Sources fiables et type analytique

| Ressource source | Objet / technique correcte | Usage dans Omics2PK |
|---|---|---|
| [BinBase / BinDiscover](https://bindiscover.metabolomics.us) et [FiehnLib](https://www.fiehnlab.ucdavis.edu/projects/FiehnLib/) | **GC-MS EI**, dérivatisation et indice de rétention documentés; BinDiscover catalogue également des inconnus | PRIORITAIRE pour les 229 signaux GC-TOF de Rotroff, **à condition de retrouver le spectre EI et l'indice RI originaux** |
| [MassBank](https://massbank.eu/) | Spectres expérimentaux de petites molécules; une partie GC, autres LC/ESI | Comparaison compatible de spectres sous conditions source et droits vérifiés |
| [GNPS](https://ccms-ucsd.github.io/GNPSDocumentation/librarysearch/) | Spectres tandem **LC-MS/MS**, fragments, précurseurs et polarité | Candidats pour LC-MS/MS. **Pas** un substitut direct à une bibliothèque GC-EI |
| [PRIDE / ProteomeXchange](https://www.ebi.ac.uk/pride/) | **Protéines, peptides, fragments de peptides, PTM** (PXD), MS protéomique | Référence à une protéine, peptide et site PTM; **non compatible** avec le moteur EI de petites molécules |
| [UniProt](https://www.uniprot.org/) / [Rhea](https://www.rhea-db.org/) / [ChEBI](https://www.ebi.ac.uk/chebi/) / [RaMP-DB](https://github.com/ncats/RaMP-DB) | Ontologie chimique et protéique, transport, enzymes, réactions documentées | Interprétation mécanistique **après**, et seulement après, une annotation chimique source vérifiée |

Une correspondance GC-EI exige **spectre de fragments** (minimum 6 ions), **même dérivatisation**, **même système d'indices de rétention**, tolérances de masses/RI explicitement déclarées. Une correspondance MS2 exige **même mode d'ionisation/polarité, masse de précurseur tolérée et spectre de fragments**. Le score calculé ici est un **cosinus exploratoire transparent**, pas la recherche officielle GNPS, MassBank ou un identifiant MSI niveau 1. L'acquisition `PROTEOMICS_MSMS` est **enregistrée au catalogue**, mais le moteur de recherche de petits métabolites **la refuse**; PTM protéique nécessite une voie d'annotation peptide–site–protéforme distincte.

**Refus explicites :** pas de spectre → `ABSTAIN_RAW_SPECTRAL_PEAKS_MISSING`; mauvais RI/dérivatisation/polarité → aucun match; pas de relation moléculaire sourcée → aucun transporteur suggéré; un match plausible ≠ structure confirmée, fonction enzymatique, variation de CL ni bénéfice clinique. Même le partage de bin id entre études ne prouve pas la même molécule sans protocole comparable.

## Références mécanistiques préqualifiées comme **hypothèses**, pas comme spectres

La commande `--seed-literature` importe **neuf liens bibliographiques sans importer un seul spectre** :
- **CHEBI:16797**, N¹-méthylnicotinamide (NMN), et **CHEBI:16020**, N¹-méthyladénosine (m¹A), associées expérimentalement et cliniquement à l'axe rénal **OCT2/SLC22A2 — MATE1/SLC47A1 — MATE2-K/SLC47A2**. Source [Nishii et al. 2025, métformine fedratinib et biomarqueurs](https://pubmed.ncbi.nlm.nih.gov/40445153/). La référence **ne sépare pas à elle seule** la contribution de chaque transporteur individuel.
- **CHEBI:61049**, tacrolimus, littérature **ABCB1** intestinal / **CYP3A4–CYP3A5** métabolisation intestinale et hépatique. Source [PharmGKB, cyclosporine/tacrolimus pathways](https://pmc.ncbi.nlm.nih.gov/articles/PMC4119065/). Une telle relation **n'implique pas** qu'un ion inconnu chez un volontaire metformine soit du tacrolimus.
- Identifiants chimiques `CHEBI:16797`, `CHEBI:16020`, `CHEBI:61049` recoupés sur [PubChem CID 457](https://pubchem.ncbi.nlm.nih.gov/compound/457), [CID 27476](https://pubchem.ncbi.nlm.nih.gov/compound/27476) et [CID 445643](https://pubchem.ncbi.nlm.nih.gov/compound/445643).

D'autres graines (OAT1/OAT3, albumine PTM, POR, CYP, UGT, BCRP, etc.) ne doivent être ajoutées qu'avec identifiants source et littérature correcte, en distinguant les organes.

## Comment utiliser le prototype hors GitHub

```bash
python -m pip install -e .
python -m scripts.query_local_spectral_bank --list-sources

# Références propres/explicitement autorisées (ne pas committer ces fichiers).
python -m scripts.query_local_spectral_bank \
  --db /private/spectral_knowledge.sqlite --seed-literature \
  --reference-json /private/licensed_ei_reference.json \
  --query-json /private/unknown_ei_spectrum.json \
  --mz-tolerance-da 0.3 --ri-tolerance 100 --min-matched-peaks 6
```

La référence `licensed_ei_reference.json` doit comporter `library_source`, `source_record_accession`, `source_url`, `reuse_permission=explicitly_permitted`, `chemical_id` (vrai ChEBI/HMDB/InChIKey), `acquisition=GC_EI`, `ionization=EI`, `derivatization`, `retention_index`, `retention_index_system`, `instrument_method`, et `peaks=[[m/z,intensité],...]`. La requête contient `acquisition`, `ionization`, `derivatization`, `retention_index`, `retention_index_system` et les `peaks`. Les tolérances se **justifient à partir de la méthode et de la résolution**, pas de la PK observée.

L'interface simple reste un protocole de recherche et ne remplace pas les algorithmes expérimentaux de déconv./alignement, annotation MSI, séparation d'isomères, ou validation par standard.

## Qualification réelle déjà exécutée

La [CI sur les suppléments originaux Rotroff](https://github.com/rberrah/Omics2PK/actions/runs/41909996908) est **verte** : le script `scripts/register_rotroff33_unresolved_spectral_features.py` a lu les deux originaux XLSX liés par SHA-256 et enregistré **384 headers d'analyses individuelles, 155 noms chimiques et 229 bins inconnus** en SQLite éphémère **sans aucune valeur participant** ni correspondance spectrale fictive. Statut des 229 = `UNRESOLVED_BIN_ID_ONLY_NO_SPECTRUM`. Pour construire localement la base persistante : `python -m scripts.register_rotroff33_unresolved_spectral_features --db /private/rotroff_spectral_catalogue.sqlite`.

Le module `src/omics2pk/proteoform_reference.py` permet d'enregistrer des métadonnées **protéiques** PXD/UniProt, peptide ou protéforme, et PTM uniquement avec niveau de preuve explicitement déclaré. Les modifs non localisées sont autorisées comme observations `unlocalized`, mais jamais comme sites confirmés ni comme preuve de liaison albumine–médicament.

## Rotroff — ce qui est réellement disponible

- **384 signaux GC-TOF prédose** appariés aux 33 volontaires, dont **155 noms chimiques** et **229 bins numériques sans nom**; 32 personnes disposent d'une mesure de PK.
- Le **tableau Excel qualifié dans les suppléments** contient les **intensités** des bin IDs, **pas des spectres EI déconv. par bin**, de masses m/z de fragments ni d'indices de rétention exploitables. Donc le prototype **n'associe encore aucun des 229 signaux** à une structure ou à un transporteur.
- La prochaine opération de données est un **crosswalk original** des identifiants dans BinBase/BinDiscover ou la récupération des spectres d'origine et de leur RI auprès du laboratoire. Ne pas supposer qu'un numéro (`629905` par exemple) encode une masse de 629.905 Da : c'est **un identifiant**, pas une masse.
- L'analyse par voie métabolique des **155 annotations nommées** peut être réalisée indépendamment de cette récupération. Pour les **229 inconnus**, garder les intensités comme variables prédictives, mais les annotations mécanistiques restent `UNRESOLVED`.

## Niveau de preuve avant d'émettre une interprétation

- **E0** : identifiant de signal et intensité seulement : candidat statistique sans structure.
- **E1** : spectre/RI d'origine, recherche de similaires, familles chimiques hypothétiques.
- **E2** : métabolite putativement identifié par référence spectrale et RI suffisamment comparables, annotation structurale candidate (pas encore standard confirmé).
- **E3** : confirmation par standard / analyse structurale et provenance du composé; association à une voie/transporteur avec **référence expérimentale distincte**.
- **E4** : signal sanguin réellement associé à une **activité fonctionnelle du transporteur/enzymes** dans le tissu pertinent.
- **E5** : valeur ajoutée de ce signal sur une **PK humaine mesurée hors échantillon**, contre clinique + PGx + PopPK, avec validation indépendante.

**MultiOmics** doit conserver l'identité/qualité et la provenance du signal; **Omics2PK** utilise la couche réactionnelle pour construire des candidats fonctionnels et tester leurs valeurs dans un vrai modèle PK. On ne doit pas créer de nouvelle application indépendante pour ce catalogue dans l'immédiat.
