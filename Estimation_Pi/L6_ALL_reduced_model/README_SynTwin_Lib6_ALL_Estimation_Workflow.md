# SynTwin -- L6 Estimation Workflow (ALL, Reduced Model)

This folder implements the complete SynTwin workflow for parameter
estimation and results generation of the **L6 sublibrary** under the
**ALL scheme** using the **reduced model**.

------------------------------------------------------------------------

## L6 Context

L6 is a 6‑TU sublibrary defined by:

-   **2 plasmid origins**
-   **3 promoters**
-   **1 shared RBS** (single RBS across all constructs)

Total transcriptional units (TUs): 2 × 3 × 1 = 6

### Key methodological feature

In L6:

-   Only the **shared RBS intrinsic initiation capacity**\
    `k0_sigma0` is estimated.
-   The RBS sensitivity parameter `inv_sigma0` is fixed.
-   All **origin- and promoter-dependent parameters**\
    (e.g., `Omega`, `Gene_cn`) are **inherited from Lib24** inference
    results.
-   Experimental measurements are retrieved from the Lib30-format tensor
    container and restricted to the L6 indices.

Thus, L6 estimation propagates uncertainty from Lib24 while identifying
the RBS parameter specific to this sublibrary.

------------------------------------------------------------------------

## Workflow Structure

The workflow consists of three stages.

------------------------------------------------------------------------

### Stage I --- Parameter Estimation

`Estimate_L6_ALL_reduced_model.m`

-   Runs `num_runs` independent BADS optimizations.
-   For each run, selects the corresponding Monte Carlo sample (`Omega`,
    `Gene_cn`) inherited from Lib24.
-   Estimates the shared RBS parameter `k0_sigma0`.
-   Propagates inherited-parameter uncertainty across runs.

**Output:**

    ./Estimated_results/Results_BADS_Lib6_ALL_reduced_<Use_mean>.mat

Each entry contains:

    Results_BADS_B0034_ALL_approx{num_run}.results = [k0_sigma0, Jmin]

------------------------------------------------------------------------

### Stage II --- Results Compilation

`Generate_Results_L6_ALL_reduced_model.m`

-   Loads L6 estimation outputs.
-   Loads inherited Lib24 parameter statistics.
-   Retrieves experimental data for the L6 indices.
-   Builds `Results_Tensor_Lib6_ALL_reduced`.
-   Computes:
    -   Digital-twin predictions
    -   Local sensitivities
    -   Monte Carlo propagated predictions
    -   Kernel density--based predictive statistics

**Output:**

    Results_Tensor_Lib6_ALL_reduced_model_<Use_mean>.mat

This structure contains:

-   Inherited parameter statistics (Omega, Gene_cn)
-   Estimated RBS statistics (raw, mean, std)
-   Experimental Mu/Pi trajectories
-   Deterministic predictions
-   Monte Carlo predictive distributions

For detailed description of the tensor structure, refer to the *Software
and Data Supplementary* document.

------------------------------------------------------------------------

### Stage III --- Visualization

`Show_Results_L6_ALL_reduced_model.m`

-   Loads the generated tensor.
-   Produces diagnostic plots and summary visualizations.

------------------------------------------------------------------------

## How to Run

From any working directory:

``` matlab
ROOT = init_SynTwin('experimental',true,'bads',true);
addpath(fileparts(mfilename('fullpath')));

Estimate_L6_ALL_reduced_model
Generate_Results_L6_ALL_reduced_model
Show_Results_L6_ALL_reduced_model
```

------------------------------------------------------------------------

## Configuration

Inside the scripts, the user must set:

-   `Use_mean`
    -   `'Global'`\
    -   `'Instances'`\
    -   `'Wells'`
-   `num_runs`\
    Number of Monte Carlo / optimization repetitions (e.g., 100).

------------------------------------------------------------------------

## Dependencies

-   HEM surrogate model\
    `Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat`

-   Experimental tensor container (Lib30 format)

-   Inherited Lib24 parameter tensor\
    `Results_Tensor_Lib24_L1O_reduced_*.mat`

-   BADS optimizer (bundled in SynTwin)

------------------------------------------------------------------------

## Scientific Role within SynTwin

L6 provides a controlled scenario where:

-   The cellular context (Ori + Promoter) is fixed by inheritance.
-   Only the RBS intrinsic initiation capacity is re-identified.
-   Uncertainty from previously characterized parts (Lib24) is
    explicitly propagated into the new inference.

This makes L6 a bridge between independent-library characterization and
context-aware parameter transfer within the SynTwin framework.

------------------------------------------------------------------------

## Reproducibility

All outputs produced in this folder correspond to the computational
workflow described in the associated manuscript and Supplementary
Information.
