# CYP3A5 functional-evidence reconciliation pilot (2026-11-03)

**Scope:** research-only, qualitative rule/uncertainty evaluation. No clinically actionable dosing, no clearance model, no fitted biological score, and no claim of superiority versus PGx. Tested fixtures are *synthetic*, not human omics/PK outcomes.

## Scientific assertion versus measured quantity

| Layer | What can actually be established | What must NOT be inferred automatically |
|---|---|---|
| Qualified PGx diplotype | CYP3A5*1 normal function, *3/*6/*7 no-function allele categories (CPIC); *1/*1, *1/no-function and no-function/no-function are expected expressers/intermediate/nonexpressers | Fixed 0%/75%/95% activity losses, absolute CL/F, *1 from an untyped or partial panel |
| Total-gene RNA | Transcript sequences or reads are detected by a specified assay | Normal splicing, full-length CYP3A5 protein, enzymatic activity |
| Splice-junction RNA | Canonical or aberrant splicing signal, subject to mapping/isoform-specific validation | Normal protein or activity solely from canonical RNA (residual *3 transcripts can exist) |
| Proteomics | Pan-CYP3A, protein-specific peptide or full-length signal according to peptide coverage | Functional protein from a non-specific peptide, proportional catalytic activity |
| Enzyme activity | A substrate-specific assay can support functional phenotyping with controls and substrate specificity | CYP3A5-specific activity from shared CYP3A4/5 activity |
| Metabolomics | Abundance reflects formation, further metabolism, transport and elimination; isotope tracer can estimate some fluxes | Direct CYP3A5 activity from a steady-state concentration or a metabolite with competing routes |

CYP3A5*3 is a splice defect (rs776746); *6 is a synonymous splice-altering change (rs10264272); **CYP3A5*7 is an insertion causing a frameshift**, not a deletion (rs41303343). CPIC/PharmVar function categories are **not quantitative percentage priors**. Full-length transcript or residual activity may exist in some settings; triage flags contradictions but never rewrites the raw observations.

The module `omics2pk.functional_knowledge.assess_cyp3a5_evidence` consumes explicitly called, optionally molecular observations and produces rule codes. `call_status="complete"` means an independently qualified and sufficiently covered diplotype; the function does not establish completeness by itself. An unresolved/partial genotype abstains. Relevant enzyme activity requires an assay-specific **validated external reference**. This is a triage contract, not a molecular classifier or automatic annotation pipeline. No subject-level data should be stored in the repository.

### Example (synthetic)

```python
from omics2pk.functional_knowledge import assess_cyp3a5_evidence
result = assess_cyp3a5_evidence(
    diplotype="*3/*3", call_status="complete",
    observations=[{"omic": "transcriptomics", "target": "total_gene",
                   "signal": "detected", "tissue": "liver"}])
assert result["status"] == "INTERPRETATION_ONLY"
assert "TOTAL_RNA_NOT_FUNCTIONAL_TRANSCRIPT" in result["evidence_warnings"]
assert result["estimated_enzyme_activity"] is None
```

### Integration boundary

- **MultiOmics (producer):** retain DOWNSTREAM_EXCHANGE_V1 and the Research v1 frozen statistical core unchanged. It can report feature/assay provenance without assigning causal function. In future, add a separate versioned optional *annotation sidecar* with exact feature IDs, tissue, assay specificity, primary evidence, ontology version and pre-specified mapping. Do **not** modify historical exchange manifests or invent isoform-specific assays from a gene label.
- **Omics2PK (consumer):** join original assays to curated annotations only with source-attested biological person/specimen/collection landmarks, per-fold preprocessing, gene/variant panel coverage, held-out participants and explicit PK endpoint. Keep uncertain/ambiguous annotations and missing genotype separate rather than imputing *1.
- **Future v0.2:** expand with evidence-curated CYP3A4 (including *22) and ABCB1 and reaction/isoenzyme mapping; drug- and tissue-specific contributions; versioned curation review, negative controls and a genuinely calibrated probabilistic activity model. A gene-level *total expression* value must never be silently recoded as functional expression.

### Prespecified quantitative benchmark (NOT yet run)

On **identical original patients with qualified clinical PK and authentic complete PGx**, compare (A) clinical, (B) clinical+PGx, (C) B+omics (single/multi-omic), (D) C+biological priors/reconciliation. Use nested patient-held-out folds or external validation, pre-dose molecular measurement and consistent outcomes. Primary contrast D vs B; also D vs C for incremental value of curation. Predeclare MAE/RMSE/calibration and paired uncertainty. Include a deliberately misleading total-RNA negative-control feature and test knowledge-ablation, tissue/source/allele-coverage sensitivity. No clinical predictive gain follows from these synthetic unit tests.

### Sources, consulted 2026-11-03

- CPIC tacrolimus CYP3A5: https://cpicpgx.org/guidelines/guideline-for-tacrolimus-and-CYP3A5/
- PharmVar CYP3A5: https://www.pharmvar.org/gene/cyp3a5
- PharmVar GeneFocus, doi:10.1002/cpt.2563: https://pmc.ncbi.nlm.nih.gov/articles/PMC9399309/
- Joint AMP/CPIC/CAP CYP3A4/CYP3A5 genotyping consensus: https://pmc.ncbi.nlm.nih.gov/articles/PMC10565868/
