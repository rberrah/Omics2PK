# Software IP, clinical confidentiality, local-only compute, and publication

**Decision 2026-10-16:** Research computation stays LOCAL. No Oracle instance, hosted Python/R service, or other server deployment is planned or authorized. Public GitHub Pages content, if later needed, is limited to static documentation and predefined public-data examples. 

- The repository was public at initial inspection. Public code does not imply public patient data and does not remove institutional IP/coauthor obligations. **No license chosen automatically.**
- Never commit BIOMARGIN, ABIS/ISBA, AADAPT, PCCP, individual PK concentration profiles, controlled-access transcriptomes, credentials or private Pillar A data.
- SWEET is upstream MIT; reuse it with original license/attribution if and when the scientific algorithm is adopted. Do not copy it into this repository yet.
- MultiOmics exchange contains pseudonymous but still potentially sensitive omic/clinical links: execute locally, use synthetic examples in CI, do not send raw matrices to external cloud by default.
- GitHub Pages cannot run Python SWEET, R/Bioconductor, XGBoost or Julia SAEM. The current rberrah.github.io MultiOmics public URLs may display public predefined examples when qualification completes.
- Run all Python SWEET/XGBoost, R/Bioconductor and Julia/SAEM jobs on authorized local or institutional machines only. Any future hosting decision requires separate explicit authorization, benchmarking and privacy review.
- Proposed license options and compatibility limits are recorded in [LICENSE_STRATEGY.md](LICENSE_STRATEGY.md). No licensing change is made automatically.
