# Omics2PK — du code génétique à l'information omique réellement utile

**5 novembre 2026 · recherche méthodologique · extension de la PR #137.**
Le pipeline ne doit pas présumer que « ARN élevé = enzyme active », que
« protéine présente = transport/catalyse », ou que « pas de mutation appelée =
fonction normale ». Il doit **d'abord qualifier la mutation et la méthode
analytique**, puis tester ce qui prédit la PK, jamais l'inverse.

## Chaîne causale / de mesure à ne pas confondre

**ADN de référence versionné → variant génomique et zygosité si connue →
transcrit exact (exons/introns, orientation, CDS, épissage) → codon →
acide aminé/protéine théorique → expression/localisation/protéforme/activité
mesurée → processus médicament-organe → vraie PK individuelle.**

À **chaque flèche**, distinguer :
1. **Calcul depuis la séquence** (par exemple `GAA → GCA` donc `Glu → Ala`) ;
2. **hypothèse biologique** (la nouvelle alanine pourrait modifier catalyse ou
   stabilité ; le stop précoce pourrait déclencher un NMD) ;
3. **effet fonctionnel validé** par mesure/source indépendante ;
4. **gain de prédiction PK** mesuré sur sujets tenus à l'écart.

Ces quatre niveaux ne sont **jamais automatiquement équivalents**.

## 1. Classification fondée sur la séquence — nouveauté exécutable

`src/omics2pk/genetic_code_variant_effect.py` fournit :

- `annotate_genomic_snv(reference_transcript, variant, genome_assay)` :
  vérifie l'assemblage GRCh37/38, le chromosome, l'orientation **+ / −**,
  les exons (coordonnées génomiques **1-based inclusives**), la longueur
  du transcrit, l'identité de la base REF sur le cDNA splicé, les offsets
  de CDS **0-based**, le cadre de lecture et les codons (table **nucléaire
  standard**, pas mitochondriale).
- **Exon codant :** classe `missense_variant`, `synonymous_variant`,
  `stop_gained`, `stop_lost` et `start_lost`, affiche **REF/ALT codon,
  acides aminés et position AA** lorsqu'ils sont calculables.
- **Intron :** distingue, à partir des exons de **ce transcrit**,
  les positions adjacentes donneuses/acceptatrices (1–2 bases), région
  d'épissage élargie, et intérieur intronique. **Aucune conséquence
  d'épissage n'est déduite d'une simple proximité spatiale** :
  quantification de jonctions RNA, isoforme et validation fonctionnelle
  indispensables.
- `annotate_coding_edit()` : cas de référence **transcrit-orientés**
  pour insertion/délétion, décalage de cadre si longueur non divisible
  par trois, changement in-frame sinon. Cette API **ne remplace pas
  l'alignement, la normalisation et la conversion des indels génomiques
  vers HGVS**.
- **NMD / stop** : pour un nouveau stop dont la position est calculable
  et dont la dernière jonction exon–exon est connue, calcule **plausible
  / frontière / non prédit**, **pas dégradation mesurée**. Un codon
  stop après la dernière jonction, une isoforme alternative, une
  réinitiation de traduction ou une structure particulière peuvent
  différer de la règle canonique.
- `annotate_variant_across_transcripts()` : conserve les annotations de
  **plusieurs isoformes du même gène**. La même substitution génomique
  peut être **faux-sens dans un transcrit** et **introniquement proche
  d'un site donneur dans un autre**. C'est **un seul variant**, pas
  deux variants indépendants. Le code ne choisit pas automatiquement
  l'isoforme la plus inquiétante.
- `crosscheck_source_vep_terms()` : compare l'effet local aux termes
  Sequence Ontology d'une annotation **VEP versionnée déclarée**.
  L'assemblage, le chromosome, la position 1-based, REF/ALT
  génomiques **dans l'orientation du génome**, le gène et le transcrit
  doivent correspondre exactement. En cas de désaccord de conséquence
  ou de locus, le pipeline **s'abstient**. Ce contrôle ne télécharge
  pas VEP et n'authentifie pas automatiquement un VCF original.
- `reconcile_curated_cyp3a5_allele()` : corrige les limites de la
  seule traduction via une **preuve externe d'allèle définissant**
  *3 (intron / accepteur cryptique), *6 (substitution **synonyme**
  `c.624G>A; p.Lys208=` **mais défaut d'épissage**), *7 (duplication /
  frameshift). La concordance avec rsID, HGVS, transcrit, conséquence
  calculée et provenance est obligatoire avant association. Cette
  méthode **ne calcule pas le diplotype complet d'un patient**.

**Important :** la séquence, les intervalles exoniques et leurs sources
sont **fourni(e)s par l'appelant et vérifié(e)s pour cohérence interne**,
mais **pas automatiquement téléchargé(e)s ou authentifié(e)s contre
GRCh38/RefSeq**. Une chaîne future utilisant Ensembl VEP/MANE et
PharmVar doit verrouiller la version de l'assemblage/transcrit,
valider les REF et conserver les conséquences sur **tous les transcrits
pertinents**, avant une affirmation sur des données originales.
Un rsID et un gène sans allèle et référence versionnée ne suffisent
pas à calculer le codon.

### Exemples fictifs de changements de codons

| Modification d'une référence **synthétique** | Code génétique | Lecture biologique permise | Meilleure preuve à rechercher |
|---|---|---|---|
| `GAA → GCA` | Glu → Ala, **missense** | L'ARN peut être normal ; la protéine peut être inactive **ou fonctionnelle** | Protéine variante, localisation, activité enzyme/transport |
| `GAA → TAA` | Glu → STOP, **nonsense** | Arrêt de traduction et NMD possibles, **conditionnels au transcrit** | Jonctions/isoformes RNA, peptides C-terminal, protéforme, activité |
| `GCT → GCC` | Ala → Ala, **synonyme** | Pas de changement AA à ce codon, **épissage non exclu** | RNA isoforme/épissage ; confrontation aux variants déjà étudiés |
| Variant proche d'un donneur intronique | Aucun codon changé **directement** | Épissage potentiellement altéré, pas prouvé | RNA junction/long-read, transcriptome |
| **Aucun variant détecté** (assay ciblé complet) | Aucun événement dans **les régions testées** | Ni absence d'autres variants ni pleine activité de la protéine | RNA/protéine/activité dans le tissu selon mécanisme |

Le scénario utilise **une fausse séquence pour un faux exemple**, pas un
transcrit CYP3A5 authentique et pas un sujet réel. Le cas *6 est **une
annotation indépendante sourcée**, pas une variante observée chez le
patient de cet exemple.

## 2. Comment les conséquences orientent les omiques

| Conséquence génomique calculée / sourcée | ARN total | ARN de jonction / isoforme | Protéine totale | Protéforme/peptide variant-spécifique | Activité fonctionnelle | Métabolites |
|---|---|---|---|---|---|---|
| Faux-sens | Peut être normal **malgré perte de fonction** | Généralement contextualisant | Peut être normal | Potentiellement utile pour protéine anormale | **Indispensable pour confirmer l'effet** | Peut renseigner phénotype, peu spécifique |
| Non-sens avec NMD plausible | Peut décroître, pas certain | Contrôle du transcrit concerné | Un peptide conservé peut persister | C-terminal / intacte distingueraient une troncation | Confirme impact | Complément fonctionnel indirect |
| Synonyme | Peut rester normal | **Important si suspicion ou preuve d'épissage** | Potentiellement normal | Non nécessaire si fonction déjà sourcée mais peut aider | Vérifie l'effet | Indirect |
| Intron proche du site d'épissage | Ambigu | **Prioritaire pour valider l'épissage** | Peut être modifié | Cible en aval de l'isoforme aberrante | Confirme en aval | Indirect |
| Non détecté dans un panel ciblé | Valeur possible, aucune causalité attribuable au variant | Selon mécanisme hors variants connus | Quantité totale informative ou non | Selon question | Peut rester le meilleur phénotype | Peut agréger inflammation/état physiologique |

Le code ne calcule aucun **pourcentage de fonction**, ni effet individuel de
mutation sur CL, ni « meilleure omique » en termes de MAE. Il propose
les couches qui peuvent répondre à une hypothèse précise, sachant que
le bénéfice prédictif dépend également du médicament, de l'organe, du
type d'échantillon, de l'assay, du nombre de patients et de la
performance du modèle clinique+PGx.

### Dangers supplémentaires

- Un même variant peut être **intronic** dans une isoforme et **exonic**
  dans une autre : l'appel doit fournir le **transcrit exact** et garder
  les annotations concurrentes ; un seul transcrit n'est pas une
  annotation universelle.
- Une variante **missense** peut être **bénigne**, altérer la liaison,
  la localisation, la stabilité, l'activité substrat-spécifique ou
  modifier l'épissage. Le code génétique **ne permet pas** de
  calculer les % de perte de fonction.
- Deux substitutions observées dans **le même codon** peuvent produire
  **une autre conséquence que deux annotations indépendantes**. Sans
  phasage et retraduction conjointe, le pipeline refuse de combiner ces
  effets au lieu d'inventer une protéine non fonctionnelle. Les SNV
  appartenant à plusieurs transcrits restent un **seul locus**.
- Les indels à plusieurs bases, variants complexes, SV/CNV, duplications
  d'exons, variants mosaïques, haplotypes et phasage exigent une
  annotation indépendante (VEP/ClinGen/PharmVar, HGVS, signatures source).
  La première version refuse une conversion approximative des indels
  génomiques.
- **CYP3A5*6** est un test négatif clé : « synonyme » au niveau de
  l'acide aminé **n'est pas synonyme de fonction normale**.
- **Tissu / plateforme / temporalité** sont indispensables pour ne pas
  confondre RNA sang total, hépatocyte, métabolite prédose,
  médicament postdose et activité de sonde mesurée in vitro.

## 3. Intégration dans la chaîne Omics2PK

L'API déjà existante
`assay_aware_physio_triage.plan_physio_omic_acquisition(...)` accepte
désormais **`sequence_variants`** (référence transcrite et exons, SNV
génomique, méthode de génotypage) et **`genome_assay_background`**
(variants absents des régions couvertes ou génotypage non réalisé).
Elle renvoie :
- `reference_checked_DNA_transcript_codon_consequences` ;
- `variant_source_consequence_interpretations` ;
- `genomic_assay_background_without_or_with_called_variants` ;
- les limites et priorités interprétatives des assays (RNA,
  peptides/protéomes, métabolites, **fonction mesurée**) ;
- un **plan de validation prédictive**, non une sélection gagnante
  déjà effectuée.

Les interactions **variant × modalité** ne doivent être apprises que
lorsqu'il y a assez de porteurs/non-porteurs et sur les **mêmes
personnes**, avec `k=0` explicitement autorisé. Une vraie PK doit
être comparée à **clinique + génotype réellement complet + PopPK**,
sur un critère prespécifié (AUC/CL/C0, log-erreur, RMSE, biais et
validation externe).

**Black-box et mécanistique restent deux circuits distincts** :
un métabolite GC-EI-TOF non identifié peut rester un prédicteur
sans jamais recevoir une fausse identité chimique, tandis que
l'annotation génétique apporte une hypothèse biologique qu'il faut
tester expérimentalement. Les deux circuits peuvent se contredire :
ce désaccord est une **question scientifique**.

## 4. Exécuter et reproduire (aucun sujet réel)

`python -m unittest discover -s tests -p 'test_genetic_code_variant_effect.py' -v`

`python -m unittest discover -s tests -p 'test_assay_aware_physio_triage.py' -v`

`python scripts/demo_genetic_code_omic_interpretation.py`

Le job GitHub [Genetic code reference checked variant consequences and assay
interpretability (offline)](https://github.com/rberrah/Omics2PK/actions)
contrôle les codons, les insertions/délétions transcrit-orientées, les
deux brins, l'assemblage et le REF, stop/missense/synonyme,
l'intron donneur/receveur, la fenêtre NMD, l'exception CYP3A5*6,
le cas sans variant identifié et l'absence de conclusions cliniques
inventées.

**Statut : contrat logiciel, références synthétiques.**
La preuve future exige l'annotation VCF/VEP authentifiée, les transcrits
MANE/RefSeq versionnés, l'inventaire des méthodes MultiOmics,
une cohorte avec génotype/phases/omiques/PK, puis validation externe.

### Documentation primaire

- [Ensembl VEP — variant effects on protein, transcripts, regulatory
  regions](https://www.ensembl.org/info/docs/tools/vep/index.html).
- [AMP/CAP/CPIC CYP3A4/5 — *3/*6/*7 and splicing](https://pmc.ncbi.nlm.nih.gov/articles/PMC10565868/).
- [PharmVar CYP3A5 nomenclature](https://www.pharmvar.org/gene/cyp3a5).
- [ClinGen variant classification PVS1 and NMD caveats](https://clinicalgenome.org/working-groups/sequence-variant-interpretation/).
