# Adaptive Gauss-Hermite — original n120 mavoglurant patient-disjoint NLME

## Result checked in original source CI (16 October 2026)

Original 120 biological humans, 196 qualified original finite-IV-infusion occasions, 2,403 positive observed blood concentrations; original publicly sourced RDA pinned by SHA-256 ba1d75e20bcc09f20a55a05a9e4424244efa15b2dfcef6079cd27385d5876174. Both occasions for a biological person remain in one of five fully held-human-out folds.

See https://github.com/rberrah/Omics2PK/actions/runs/38973121688 , 2026-10-16: all synthetic tests and original source benchmark succeeded.

| Out-of-fold endpoint | Fixed pooled 2C | 2C marginal NLME with CL IIV and adaptive quadrature |
|---|---:|---:|
| Root mean log concentration MSE, equal weight per human | **0.383383** | **0.387256** |

The new adaptive per-human posterior proposal integrates an exact finite-infusion 2C model and a single human-level lognormal clearance effect. It uses numerical mode/curvature, normalized prior and quadrature integration without ever creating per-occasion pseudodonors. The held-out **point prediction** integrates the prior without using held-out concentrations or EBEs. The **joint held-out log predictive density** mathematically integrates over each held-out person's effect using all that person's observed concentrations for scoring, not for model fitting or point prediction.

Original held-out **negative joint log predictive density per independent human** using one identical fitted parameter set per fold:

| Adaptive Gauss-Hermite nodes | Mean NLL per human |
|---|---:|
| 7 | 3.908502 |
| 11 | 3.908502 |
| 15 | 3.908502 |
| 21 | 3.908502 |

Displayed absolute range **0.000000** at the 6-decimal output resolution. This is a numerical improvement over the previous *unadapted* Gauss-Hermite NLL 7→31 variation reported in docs/MAVOGLURANT_2C_NLME_IIV_CL_PILOT.md. It is **not** a demonstration of population PK model correctness, parameter identifiability, or perfect predictive calibration.

## Scientific scope

- Single subject-shared ETA on log(CL); no independent IOV, no random Vc/Q/Vp, no residual variability by sampling time.
- Conditional Gaussian residual on log observed positive concentration, no BLQ likelihood, no experimental LLOQ or analytical censoring correction.
- Approximate marginal maximum likelihood optimized with L-BFGS-B and bounded multi-initializations, not SAEM or FOCEi; parameter recovery and uncertainty intervals are not yet validated against an independent mixed-effects estimator.
- The original mavoglurant dataset has **no omics**. No SWEET, multiomics, pillier-A causal score, in-vivo tacrolimus CL/F or benefit from dynamic omics exists here.
- Population point prediction is not necessarily improved by adding random effects. The small higher held-out RMSLE (+0.003873) is a transparent negative finding, not a material clinical deterioration claim.

Next: independent SAEM/FOCEi model implementation on the *same original events* (separate dose/infusion and MDV handling), random CL/Vc and IOV sensitivity, posterior-predictive calibration, identifiability, BLQ handling, nested clinical covariates only if observed, then E3 molecule-specific biological covariates once original same-person omic+PK arrays are verified.

Source code: src/omics2pk/adaptive_iv_nlme.py, scripts/benchmark_mavoglurant_iv_adaptive_nlme.py, tests/test_adaptive_iv_nlme.py. Run with python -m pip install -e ".[pk-structural]".
