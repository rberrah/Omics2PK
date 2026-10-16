# Software IP, clinical confidentiality, local compute, and publication

- The repository was public at initial inspection. Public code does not imply public patient data and does not remove institutional IP/coauthor obligations. **No license chosen automatically.**
- Never commit BIOMARGIN, ABIS/ISBA, AADAPT, PCCP, individual PK concentration profiles, controlled-access transcriptomes, credentials or private Pillar A data.
- SWEET is upstream MIT; reuse it with original license/attribution if and when the scientific algorithm is adopted. Do not copy it into this repository yet.
- MultiOmics exchange contains pseudonymous but still potentially sensitive omic/clinical links: execute locally, use synthetic examples in CI, do not send raw matrices to external cloud by default.
- GitHub Pages cannot run Python SWEET, R/Bioconductor, XGBoost or Julia SAEM. The current rberrah.github.io MultiOmics public URLs may display public predefined examples when qualification completes.
- Oracle Cloud Always Free is potentially viable **only for small public compute demos**. Official 2026 limits: Ampere A1 1,500 OCPU-hours and 9,000 GB-hours per month (equivalent to 2 OCPU / 12 GB continuous); quotas depend on tenancy, home region and capacity. Idle instances can be reclaimed. Require ARM package support, non-public input blocking, quotas, authentication, request limits and monitoring before hosting even a demo. **Do not use it for clinical patient uploads absent authorized institutional controls.**
- Prefer local Windows/Linux CPU runs first. If remote compute becomes needed, first measure peak RSS, time per dataset, queue requirements and re-identification/privacy risk, then consider approved infrastructure.
