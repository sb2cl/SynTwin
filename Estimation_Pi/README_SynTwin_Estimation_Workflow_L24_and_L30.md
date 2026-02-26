# SynTwin -- Estimation Workflow Folder

This folder implements a complete **SynTwin parameter-estimation
workflow** for a specific combination of:

-   **Library:** LibXX (e.g., Lib24, Lib30)
-   **Validation scheme:** ALL (full dataset) or L1O (Leave-One-Out
    cross-validation)
-   **Model type:** reduced or full

The exact configuration is encoded in the script filenames, e.g.:

    Estimate_LibXX_<Validation>_<Model>.m
    Generate_Results_LibXX_<Validation>_<Model>.m
    Show_Results_LibXX_<Validation>_<Model>.m
    J4_LogPI_LibXX_<Validation>_<Model>.m

Refer to the associated manuscript and the *Software and Data
Supplementary* document for the theoretical background and detailed
description of the underlying digital-twin model and data structures.

------------------------------------------------------------------------

## Workflow Structure

The SynTwin estimation workflow is organized into three main stages:

### Stage I --- Parameter Estimation

`Estimate_*.m`

-   Initializes SynTwin and loads the experimental data container.
-   Configures the selected validation scheme (ALL or L1O).
-   Defines the model configuration (reduced or full).
-   Calls the optimization solver (BADS or MEIGO/eSS).
-   Stores one result per estimation run (and per fold, if L1O).

**Outputs:**

    ./Estimated_results/Results_<Solver>_LibXX_<Validation>_<Model>_*.mat

Each file contains: - Estimated parameter vector(s) - Objective-function
value(s) - Metadata describing the estimation context

------------------------------------------------------------------------

### Stage II --- Results Compilation

`Generate_Results_*.m`

-   Loads optimization outputs from `Estimated_results/`
-   Reconstructs model predictions
-   Aggregates outputs at the selected level:
    -   `Global`
    -   `Instances`
    -   `Wells`
-   Builds a structured Results tensor

**Outputs:**

    Results_Tensor_LibXX_<Validation>_<Model>_<Aggregation>.mat

These `Results_Tensor_*` structures contain the complete set of:

-   Model predictions
-   Experimental comparisons
-   Summary statistics
-   Parameter distributions
-   Derived observables

For a detailed description of these data structures, refer to the
*Software and Data Supplementary* document.

------------------------------------------------------------------------

### Stage III --- Visualization

`Show_Results_*.m`

-   Loads the compiled `Results_Tensor_*` file
-   Generates diagnostic plots
-   Displays parameter estimates and goodness-of-fit metrics

------------------------------------------------------------------------

## How to Run

From any working directory:

``` matlab
% Initialize SynTwin
ROOT = init_SynTwin('experimental',true,'bads',true);

% Add this folder if running from SynTwin root
addpath(fileparts(mfilename('fullpath')));

% Run the workflow
Estimate_LibXX_<Validation>_<Model>
Generate_Results_LibXX_<Validation>_<Model>
Show_Results_LibXX_<Validation>_<Model>
```

------------------------------------------------------------------------

## Validation Schemes

-   **ALL**\
    Uses the full dataset for parameter estimation.

-   **L1O (Leave-One-Out)**\
    Repeats the estimation by leaving out one construct at a time to
    assess predictive generalization.

------------------------------------------------------------------------

## Model Variants

-   **Reduced model**\
    Estimates a subset of intrinsic parameters while fixing
    sensitivity-related parameters.

-   **Full model**\
    Estimates the complete set of intrinsic biopart parameters.

------------------------------------------------------------------------

## Notes

-   SynTwin must be initialized using `init_SynTwin()` before running
    scripts.
-   If Parallel Toolbox is unavailable, replace `parfor` with `for` in
    estimation scripts.
-   This folder represents one specific configuration; multiple such
    folders may coexist within the `Estimation_Pi/` directory.

------------------------------------------------------------------------

## Reproducibility

All outputs generated in this folder are consistent with the host-aware
digital-twin methodology described in the associated manuscript and its
Supplementary Information.
