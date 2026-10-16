# Native MultiOmics producer → Omics2PK: public synthetic integration gate

This gate calls the **actual** MultiOmics JavaScript exporter,
\`src/lib/multiomics/downstream-exchange.js\`, pinned at
\`rberrah/MultiOmics@5a71db3bce5990b96cc487c329e36b7cc055cda9\`.

Input is the three-omic **bundled synthetic demo**:
- \`demo_metadata.csv\`
- \`demo_transcriptomics.csv\`
- \`demo_proteomics.csv\`
- \`demo_metabolomics.csv\`

Workflow steps: execute the native producer in Node 24 → verify its real files
and SHA-256 using \`omics2pk.exchange.inspect_exchange\` → assert 16 distinct
biological subjects and 49 assays (17 RNA, 16 proteomics, 16 metabolomics) →
assert one subject has two RNA assays → confirm raw RNA counts are **rejected**
by \`prepare_sweet\`. The produced package is deleted after the checks,
and cannot accidentally become a GitHub artifact. No real patient data
are uploaded.

Scientific interpretation: passing confirms the two repositories exchange
a syntactically and referentially valid three-omic contract. It **does not**
prove that the unsupervised network, genotype, PK model or downstream
estimates are valid.

Important outstanding integration tasks:
- The native browser producer reports \`software.commit: unknown\`. This is
  honest but insufficient for fully frozen research traceability.
- Public MultiOmics demo RNA consists of **raw counts** and includes a
  repeated biological subject. The user must apply documented appropriate
  expression transformation and decide explicit repeat handling **outside**
  SWEET before inference. Do not silently relabel, normalize, aggregate,
  discard or treat repeated measurements as independent patients.
- The producer conservatively sets \`data_access: private\`, even for a
  publicly bundled demo. Honor that setting when operating on real studies.
- Before any clinical project, require patient-level authorization,
  contemporary sample/time matching and outcome-split safety checks.

Run only in a local, disposable checkout with the pinned MultiOmics repository
at \`upstream/multiomics\`, or use the public-only CI workflow.
