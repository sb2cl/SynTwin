# SynTwin -- Lib5 L1O Reduced-Model Estimation Workflow

## Overview

This folder contains the full SynTwin workflow for parameter estimation
and post-processing of the Lib5 sublibrary under L1O (leave-one-out)
cross-validation, using the reduced modeling setting.

In Lib5:

-   1 plasmid origin (pGreen; fixed context)
-   1 promoter (J23100; parameter to estimate)
-   5 RBS variants (indices \[1 2 3 4 5\])

Only the promoter strength parameter Ω (Omega) is estimated.\
All RBS-related parameters and gene copy number are inherited from
previous inference steps.

------------------------------------------------------------------------

## Scientific Context

Lib5 isolates a single promoter (J23100) in a fixed Ori context (pGreen)
and evaluates it across five RBS backgrounds.

Under the L1O scheme, one RBS-defined transcriptional unit (TU) is left
out at a time:

-   Ω is inferred using the remaining 4 TUs.
-   The left-out TU is excluded from the cost function.
-   This enables cross-validation and predictive assessment.

Inherited parameters:

-   RBS k0_sigma0
    -   RBSs present in Lib24 → inherited from Lib24 L1O reduced
    -   RBS B0034 (index 3) → inherited from Lib6 L1O reduced
-   Gene copy number (Gene_cn) → inherited from Lib6 L1O reduced
-   inv_sigma0 → fixed in the scripts

------------------------------------------------------------------------

# Workflow Structure

This folder contains three main scripts:

## 1. Estimate_Lib5_L1O_reduced_model.m

Runs leave-one-out cross-validation across the 5 RBS-defined constructs.

For each left-out RBS:

1.  Builds inherited parameter samples.
2.  Runs num_runs BADS optimizations.
3.  Estimates promoter strength Ω.
4.  Saves one file per fold in ./Estimated_results/

Output naming pattern:

Results_BADS_Lib5_L1O_reduced\_`<Ori>`{=html}`<Prom>`{=html}`<RBS>`{=html}\_`<Use_mean>`{=html}.mat

Each file contains:

Results_BADS_J23100_L1O_approx{num_run}.results = \[Omega, Jmin\]

Uncertainty propagation:

-   A Monte Carlo index vector is sampled.
-   Each optimization run uses one consistent inherited parameter
    sample.

------------------------------------------------------------------------

## 2. Generate_Results_Lib5_L1O_reduced_model.m

Post-processes all per-fold estimation outputs and compiles them into:

Results_Tensor_Lib5_L1O_reduced_model\_`<Use_mean>`{=html}.mat

The generated structure includes:

-   Fold-specific ("local") Ω statistics\
    (Parameters_local_raw / mean / std)
-   Global pooled Ω statistics across folds
-   Inherited RBS and Gene_cn statistics
-   Experimental Mu(t) and Pi(t) data
-   Digital-twin predictions
-   Absolute sensitivities
-   Monte Carlo propagated predictions
-   Kernel density estimates and prediction quantiles

------------------------------------------------------------------------

## 3. Show_Results\_\* (optional)

Visualization scripts load the generated tensor and produce:

-   Pi--μ characteristic curves
-   Confidence envelopes
-   Sensitivity plots
-   Cross-validation diagnostics

------------------------------------------------------------------------

# Cost Function

J4_LogPI_Lib5_L1O_reduced.m minimizes a weighted logarithmic mismatch:

| log10(Pi_pred) \* log10(Pi_exp) − log10(Pi_exp)\^2 \|

The cost is:

-   Aggregated across all training TUs (4 per fold)
-   Averaged over time points
-   Deterministic for fixed (params, inherited sample)

------------------------------------------------------------------------

# Dependencies

Required components:

-   SynTwin initialization (init_SynTwin)
-   HEM surrogate: Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
-   Experimental data: Experimental_Data/ExpData_Tensor_lib5_micro.mat
-   Inherited parameter sources:
    -   Lib24 L1O reduced
    -   Lib6 L1O reduced
-   BADS optimizer (bundled with SynTwin)
-   Parallel Computing Toolbox (optional)

------------------------------------------------------------------------

# Configuration Parameters

Inside Estimate_Lib5_L1O_reduced_model.m:

Use_mean = 'Global' \| 'Instances' \| 'Wells'\
num_runs = number of BADS runs per fold

Recommended:

-   Use_mean = 'Wells'
-   num_runs ≈ 50--100

------------------------------------------------------------------------

# Outputs

Final compiled results:

Results_Tensor_Lib5_L1O_reduced_model\_`<Use_mean>`{=html}.mat

Main variable:

Results_Tensor_Lib5_L1O_reduced

Contains:

-   Local L1O Ω statistics per fold
-   Global Ω distribution
-   Inherited parameter statistics
-   Experimental trajectories
-   Deterministic predictions
-   Monte Carlo propagated predictive envelopes

------------------------------------------------------------------------

# Interpretation

This workflow enables:

-   Promoter strength estimation under controlled Ori context
-   Cross-validation across RBS backgrounds
-   Quantification of predictive capability
-   Explicit uncertainty propagation from upstream inference steps

Lib5 L1O acts as a validation layer in the hierarchical SynTwin
inference pipeline:

Lib24 → Lib6 → Lib5

------------------------------------------------------------------------

# Reproducibility Notes

-   All paths are portable (no absolute paths required).
-   Results depend on Monte Carlo sampling of inherited parameters.
-   For deterministic reproducibility, fix the random seed before
    sampling.

------------------------------------------------------------------------

© SynTwin Digital-Twin Framework
