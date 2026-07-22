# SynTwin parameter-estimation workflows

## Overview

`Estimation_Pi/` contains the host-aware parameter-estimation workflows used to
characterize the background signal, plasmid context, promoters, and RBS
bioparts from the SynTwin combinatorial libraries.

The workflows are organized hierarchically. Parameters inferred from one
library are propagated to smaller downstream libraries, where the remaining
unknown biopart parameters are identified in a controlled genetic context.

The main inference sequence is:

```text
pARKA21 background model
        |
        v
Lib24 reduced model
        |
        +----------------------+
        |                      |
        v                      v
Lib6 reduced model       Lib24 robustness analyses
        |
        v
Lib5 reduced model
```

Lib30 provides complementary all-library analyses:

```text
Lib30 full model
Lib30 reduced-model robustness to d_mA
```

The distributed workflows use BADS multistart estimation and the objective
functions named `J5w_*`. Unless explicitly stated otherwise, the canonical
compiled results use the `Wells` aggregation level.

---

## Directory map

The reviewed estimation subfolders are:

```text
Estimation_Pi/
├── pARKA21_background_model/
├── L24_ALL_reduced_model/
├── L24_ALL_reduced_model_dmA/
├── L24_L1O_reduced_model/
├── L30_ALL_full_model/
├── L30_ALL_reduced_model_dmA/
├── L6_ALL_reduced_model/
├── L6_ALL_reduced_model_dmA/
├── L6_L1O_reduced_model/
├── L5_ALL_reduced_model/
└── L5_L1O_reduced_model/
```

Each subfolder contains its own README with the precise script names, variable
names, result-tensor fields, and distribution policy.

---

## Common workflow pattern

Most subfolders follow three stages.

### 1. Parameter estimation

```text
Estimate_*.m
```

The estimation script:

- initializes SynTwin;
- loads the HEM surrogate and processed experimental data;
- selects the library and model configuration;
- configures parameter bounds;
- runs repeated BADS optimizations, usually with `parfor`;
- stores the parameter vector and objective value from each run in
  `Estimated_results/`.

For leave-one-out workflows, one file is generated for each excluded
construct.

### 2. Result compilation

```text
Generate_Results_*.m
```

The generator loads the BADS outputs and constructs a result tensor or result
structure containing:

- raw parameter samples;
- parameter means and standard deviations;
- objective-function values;
- processed experimental trajectories;
- deterministic synthesis-rate predictions;
- local sensitivities;
- Monte Carlo parameter samples;
- propagated prediction intervals and quantiles.

### 3. Visualization

```text
Show_Results_*.m
```

The visualization script loads the compiled result tensor, validates the
expected variable, and calls a shared plotting function from `Scripts_base/`
or from the corresponding sibling workflow.

Figures are normally exported to:

```text
Figures/
```

as PNG and SVG files.

---

## Experimental aggregation levels

The library-scale objective functions support:

```text
Global
Instances
Wells
```

### `Global`

Uses one global mean trajectory per construct.

### `Instances`

Uses one mean trajectory per experimental instance.

### `Wells`

Uses the trajectories from individual culture wells.

The public compiled tensors reviewed here use:

```matlab
Use_mean = 'Wells';
```

`Use_mean` must remain consistent across estimation, result generation, and
visualization.

---

## Common reduced-model parameterization

For the reduced RBS model, the RBS intrinsic initiation capacity is represented
as:

```text
kappa^0 = k0/sigma0
```

and the inverse scaling parameter is:

```text
rho^0 = 1/sigma0
```

The canonical reduced-model workflows fix:

```text
rho^0 = 0.02
```

and estimate `kappa^0`.

The complete Lib30 model estimates both `kappa^0` and `rho^0`.

The pSC101 effective copy number is fixed at:

```text
N_pSC101 = 5
```

whereas the pGreen/pSC101 copy-number multiplier is inferred when included in
the fitted parameter vector.

---

## Objective functions

### Library parameter estimation

The `J5w_*` functions compare predicted and experimental synthesis rates in
base-10 logarithmic space:

```text
r = log10(Pi_pred) - log10(Pi_exp)
```

They apply a pseudo-Huber loss:

```text
rho(r) = delta^2 * (sqrt(1 + (r/delta)^2) - 1)
```

with a smooth synthesis-rate-dependent weight:

```text
weight = sqrt(1 + log10(1 + Pi_exp))
```

The cost is averaged within the selected aggregation layer and then across the
constructs included in the current fit.

For leave-one-out workflows, the denominator is the number of training
constructs, not the total library size.

### pARKA21 background estimation

The background workflow uses an instance-balanced mean squared residual. The
cost is calculated independently for each experiment and then averaged across
valid instances.

---

# Background-model estimation

## `pARKA21_background_model`

The nonexpressing construct pARKA21, also identified as `PARKA` in some Cytation
layouts, is used to estimate the growth-dependent autofluorescence synthesis
rate.

The model is:

```text
x = 100 mu
Pi_bk(mu) = b0 + A(1 - exp(-x/k)) + q1 x + q2 x^2
```

with parameter order:

```text
[b0, A, k, q1, q2]
```

### Distributed experimental input

The raw Cytation Excel workbooks are not distributed. They are converted into:

```text
Experimental_Data/ExpData_pARKA21.mat
```

containing:

```text
Data_pARKA21
```

The stored trajectories are generated with `consider_BK = 0`, so no previously
fitted background has been subtracted.

### Public workflow

```matlab
Estimate_pARKA21_background_model
Generate_Results_pARKA21_background_model
Show_Results_pARKA21_background_model
```

The estimator writes:

```text
Estimated_results/
    Results_BADS_pARKA21_background_Wells.mat
```

The generator writes:

```text
Results_pARKA21_background_Wells.mat
```

and preserves the downstream-compatible vector:

```text
Params_global_BKmodel_mean
```

The default estimation and plotting workflows discard the first four and last
three trajectory points. The same cropping is used in the objective and the
figures.

---

# Lib24 workflows

Lib24 contains:

```text
2 plasmids x 3 promoters x 4 RBSs = 24 constructs
```

RBSs:

```text
B0030
B0032
J61100
J61101
```

B0034 is not part of Lib24.

## `L24_ALL_reduced_model`

This is the main joint reduced-model characterization workflow.

Estimated parameters:

```text
3 promoter transcription rates:
    J23106
    J23102
    J23101

4 RBS kappa^0 values:
    B0030
    B0032
    J61100
    J61101

1 pGreen/pSC101 copy-number multiplier
```

Total:

```text
8 free parameters
```

Canonical scripts:

```matlab
Estimate_Lib24_ALL_reduced_model
Generate_Results_Lib24_ALL_reduced_model
Show_Results_Lib24_ALL_reduced_model
```

Distributed compiled result:

```text
Results_Tensor_Lib24_ALL_reduced_model_Wells.mat
```

This tensor supports downstream plotting, Jacobian analysis, and
cross-library parameter transfer.

### Fixed-rho robustness

The same subfolder contains a complementary scan over:

```text
rho^0 = 0.005, 0.01, 0.02, 0.03, 0.05
```

using:

```matlab
Estimate_Lib24_ALL_reduced_model_set_sigma
Generate_Results_Lib24_ALL_reduced_model_set_sigma
Show_Results_Lib24_ALL_reduced_model_set_sigma
Show_Compared_Results_Lib24_ALL_reduced_model_set_sigma
Show_Robustness_rho_Lib24_ALL_reduced_model
```

The `rho^0 = 0.02` case is the canonical reference.

Large rho-specific tensors are regeneration products and are not all included
in the compact public distribution.

## `L24_L1O_reduced_model`

This workflow repeats the reduced-model estimation 24 times, leaving one
construct out in each fold.

Each fit uses:

```text
23 training constructs
1 validation construct
```

The same eight parameters are estimated as in the ALL workflow.

Scripts:

```matlab
Estimate_Lib24_L1O_reduced_model
Generate_Results_Lib24_L1O_reduced_model
Show_Results_Lib24_L1O_reduced_model
```

Distributed result:

```text
Results_Tensor_Lib24_L1O_reduced_Wells.mat
```

This tensor is the primary upstream source for Lib6. It contains both:

- local fold-specific distributions;
- pooled distributions across all L1O folds.

The pooled distributions summarize the set of L1O fits and must not be
interpreted as an additional ALL fit.

## `L24_ALL_reduced_model_dmA`

This workflow assesses robustness to the fixed non-ribosomal mRNA degradation
rate:

```text
d_mA = 0.15, 0.17, 0.20, 0.22, 0.25 min^-1
```

The reference is:

```text
d_mA = 0.20 min^-1
```

For each fixed value, the same eight reduced-model parameters are estimated.

Scripts include:

```matlab
Estimate_Lib24_ALL_reduced_model_set_dma
Generate_Results_Lib24_ALL_reduced_model_set_dma
Show_Results_Lib24_ALL_reduced_model_set_dma
Show_Compared_Results_Lib24_ALL_reduced_model_set_dma
Show_Robustness_dma_Lib24_ALL_reduced_model
Analyze_dmA_ratio_all_TUs
```

The distribution retains the estimation outputs. Large per-\(d_{mA}\) result
tensors are regenerated locally.

---

# Lib6 workflows

Lib6 contains:

```text
2 plasmids x 3 promoters x 1 RBS = 6 constructs
```

The shared RBS is:

```text
B0034
```

Promoter rates and effective plasmid copy numbers are inherited from the Lib24
leave-one-out reduced-model tensor.

## `L6_ALL_reduced_model`

Only the B0034 `kappa^0` value is estimated.

For optimization run `num_run`, all six constructs use the corresponding
Lib24 inherited samples:

```text
Omega_MC_samples(num_run)
Gene_cn_MC_samples(num_run)
```

Scripts:

```matlab
Estimate_Lib6_ALL_reduced_model
Generate_Results_Lib6_ALL_reduced_model
Show_Results_Lib6_ALL_reduced_model
```

Distributed result:

```text
Results_Tensor_Lib6_ALL_reduced_Wells.mat
```

Expected tensor dimensions:

```text
2 x 3 x 1
```

## `L6_L1O_reduced_model`

This workflow performs six leave-one-construct-out fits.

Each fold uses:

```text
5 training constructs
1 validation construct
```

Only B0034 `kappa^0` is estimated.

Scripts:

```matlab
Estimate_Lib6_L1O_reduced_model
Generate_Results_Lib6_L1O_reduced_model
Show_Results_Lib6_L1O_reduced_model
```

Distributed result:

```text
Results_Tensor_Lib6_L1O_reduced_Wells.mat
```

This tensor is the upstream source for the Lib5 B0034 and pGreen inherited
samples.

## `L6_ALL_reduced_model_dmA`

This workflow repeats the Lib6 ALL estimation at:

```text
d_mA = 0.15, 0.17, 0.20, 0.22, 0.25 min^-1
```

Only B0034 `kappa^0` is estimated.

Scripts:

```matlab
Estimate_Lib6_ALL_reduced_model_set_dma
Generate_Results_Lib6_ALL_reduced_model_set_dma
Show_Results_Lib6_ALL_reduced_model_dma
Show_Robustness_dma_Lib6_ALL_reduced_model
Show_Robustness_dma_Lib24_plus_Lib6_ALL_reduced_model
```

The compiled per-\(d_{mA}\) tensors are approximately 6 MB each and are treated
as regenerable products. The public distribution retains the BADS estimation
outputs.

The combined robustness script joins:

- Lib24 promoter, plasmid, and four RBS estimates;
- the B0034 estimate from Lib6.

---

# Lib5 workflows

Lib5 contains:

```text
1 plasmid: pGreen
1 promoter: J23100
5 RBS contexts
5 constructs
```

The parameter inferred from Lib5 is:

```text
J23100 Omega
```

Inherited quantities:

```text
pGreen effective copy number
B0030, B0032, J61100, J61101 kappa^0 from Lib24
B0034 kappa^0 from Lib6
```

## `L5_ALL_reduced_model`

All five constructs are used simultaneously to estimate J23100 `Omega`.

Scripts:

```matlab
Estimate_Lib5_ALL_reduced_model
Generate_Results_Lib5_ALL_reduced_model
Show_Results_Lib5_ALL_reduced_model
```

Distributed result:

```text
Results_Tensor_Lib5_ALL_reduced_Wells.mat
```

Expected tensor dimensions:

```text
1 x 1 x 5
```

Each BADS run uses one consistent inherited realization across all five
constructs.

## `L5_L1O_reduced_model`

Five leave-one-construct-out folds are performed.

Each fold uses:

```text
4 training constructs
1 validation construct
```

Scripts:

```matlab
Estimate_Lib5_L1O_reduced_model
Generate_Results_Lib5_L1O_reduced_model
Show_Results_Lib5_L1O_reduced_model
```

Distributed result:

```text
Results_Tensor_Lib5_L1O_reduced_Wells.mat
```

The local distributions correspond to fitting J23100 with one RBS context
excluded. The pooled distribution combines all five L1O folds and is not an
independent ALL fit.

---

# Lib30 workflows

Lib30 contains:

```text
2 plasmids x 3 promoters x 5 RBSs = 30 constructs
```

RBSs:

```text
B0030
B0032
B0034
J61100
J61101
```

## `L30_ALL_full_model`

This workflow fits the complete model.

Estimated parameters:

```text
3 promoter transcription rates
5 RBS kappa^0 values
5 RBS rho^0 values
1 pGreen/pSC101 copy-number multiplier
```

Total:

```text
14 free parameters
```

Scripts:

```matlab
Estimate_Lib30_ALL_full_model
Generate_Results_Lib30_ALL_full_model
Show_Results_Lib30_ALL_full_model
```

Distributed result:

```text
Results_Tensor_Lib30_ALL_full_model_Wells.mat
```

The visualization requires:

```text
Scripts_base/Plot_Results_Lib_complete_model.m
```

## `L30_ALL_reduced_model_dmA`

This workflow assesses the reduced Lib30 model at:

```text
d_mA = 0.15, 0.17, 0.20, 0.22, 0.25 min^-1
```

Estimated parameters:

```text
3 promoter transcription rates
5 RBS kappa^0 values
1 pGreen/pSC101 copy-number multiplier
```

Total:

```text
9 free parameters
```

Scripts include:

```matlab
Estimate_Lib30_ALL_reduced_model_set_dma
Generate_Results_Lib30_ALL_reduced_model_set_dma
Show_Robustness_dma_Lib30_ALL_reduced_model
Analyze_dmA_ratio_Lib30_ALL_TUs
```

The public distribution retains the BADS estimation outputs. The large
per-\(d_{mA}\) tensors are regenerated locally.

---

## Parameter-transfer graph

The cross-library dependencies are:

```text
L24_L1O_reduced_model
    |
    +--> L6_ALL_reduced_model
    |
    +--> L6_L1O_reduced_model
            |
            +--> L5_ALL_reduced_model
            |
            +--> L5_L1O_reduced_model
```

More explicitly:

```text
Lib24 -> Lib6
    inherited:
        promoter Omega samples
        pGreen and pSC101 effective copy-number samples

Lib24 + Lib6 -> Lib5
    inherited:
        pGreen effective copy-number samples
        B0030, B0032, J61100, J61101 kappa^0 from Lib24
        B0034 kappa^0 from Lib6
```

The matching inherited Monte Carlo realization is used consistently within
each optimization run.

---

## Recommended execution order

A complete reconstruction should follow:

```matlab
% 1. Background
Estimate_pARKA21_background_model
Generate_Results_pARKA21_background_model
Show_Results_pARKA21_background_model

% 2. Main Lib24 characterization
Estimate_Lib24_ALL_reduced_model
Generate_Results_Lib24_ALL_reduced_model
Show_Results_Lib24_ALL_reduced_model

% 3. Lib24 leave-one-out tensor required downstream
Estimate_Lib24_L1O_reduced_model
Generate_Results_Lib24_L1O_reduced_model
Show_Results_Lib24_L1O_reduced_model

% 4. B0034 characterization
Estimate_Lib6_ALL_reduced_model
Generate_Results_Lib6_ALL_reduced_model
Show_Results_Lib6_ALL_reduced_model

Estimate_Lib6_L1O_reduced_model
Generate_Results_Lib6_L1O_reduced_model
Show_Results_Lib6_L1O_reduced_model

% 5. J23100 characterization
Estimate_Lib5_ALL_reduced_model
Generate_Results_Lib5_ALL_reduced_model
Show_Results_Lib5_ALL_reduced_model

Estimate_Lib5_L1O_reduced_model
Generate_Results_Lib5_L1O_reduced_model
Show_Results_Lib5_L1O_reduced_model

% 6. Independent Lib30 full-model analysis
Estimate_Lib30_ALL_full_model
Generate_Results_Lib30_ALL_full_model
Show_Results_Lib30_ALL_full_model
```

Robustness workflows can be run after the corresponding canonical library
workflow.

---

## Distribution policy

The public MATLAB distribution follows these rules:

### Included when practical

- estimation scripts;
- objective functions;
- result generators;
- visualization scripts;
- subfolder READMEs;
- `Estimated_results/`;
- canonical compiled Wells tensors when their size is moderate.

### Regenerated locally when large

Large collections of robustness tensors are omitted when they can be recreated
from the distributed optimization outputs.

Examples:

```text
Lib24 rho^0-specific tensors
Lib24 d_mA-specific tensors
Lib30 d_mA-specific tensors
Lib6 d_mA-specific tensors
```

### Experimental data

Raw Excel workbooks are not distributed. Processed experimental data are
provided as `.mat` files under:

```text
Experimental_Data/
```

This includes:

```text
ExpData_pARKA21.mat
ExpData_Tensor_lib5_micro.mat
ExpData_Tensor_lib24_micro.mat
ExpData_Tensor_lib30_micro.mat
```

---

## Result-tensor interpretation

Compiled result tensors generally contain:

```text
TU_Ori
TU_Promoter
TU_RBS
TU_Bioparts
TU_Name
TU_color_code

Gene_cn_*
Omega_*
RBS_k0_sigma0_*
RBS_inv_sigma0_*

Mu_mumax_pmax_*
Pi_mumax_pmax_*
Instances
Synthesis_predictions
MC_samples
MC_mu_slices
```

Leave-one-out tensors additionally contain local fold-specific fields such as:

```text
Parameters_local_raw
Parameters_local_mean
Parameters_local_std
J_local_raw
```

Pooled L1O distributions combine all folds and should not be interpreted as an
additional fit using the full dataset.

For complete field definitions, refer to the SynTwin Software and Data
Supplementary documentation.

---

## Shared dependencies

Common data and model resources include:

```text
Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
Experimental_Data/ExpData_pARKA21.mat
Experimental_Data/ExpData_Tensor_lib5_micro.mat
Experimental_Data/ExpData_Tensor_lib24_micro.mat
Experimental_Data/ExpData_Tensor_lib30_micro.mat
```

Common functions include:

```text
init_SynTwin
SynTwin_path
Get_synthesis_predictions
Get_synthesis_predictions_lite
Plot_Results_Lib_reduced_model
Plot_Results_Lib6_reduced_model
Plot_Results_Lib5_reduced_model
Plot_Results_Lib_complete_model
```

The exact plotting dependency is documented in each subfolder README.

---

## MATLAB requirements

The workflows may require:

- MATLAB with support for local functions in scripts;
- BADS;
- Parallel Computing Toolbox for `parfor`;
- Statistics and Machine Learning Toolbox;
- `exportgraphics`;
- the truncated multivariate-normal sampler distributed with SynTwin.

When Parallel Computing Toolbox is unavailable, estimation loops can be
changed from `parfor` to `for`, at the cost of longer execution time.

GNU Octave compatibility has not been tested.

---

## Reproducibility notes

BADS starts from randomized initial points, and several workflows randomly
select inherited Monte Carlo samples. Exact reruns therefore require setting a
random seed before estimation, for example:

```matlab
rng(1,'twister');
```

The distributed `Estimated_results/` and canonical Wells tensors correspond to
the retained analyses and can be used without rerunning the optimizations.

When regenerating a result tensor, verify that:

- the matching `Use_mean` value is selected;
- all required upstream tensors exist;
- the expected estimation filenames are present;
- the number of requested Monte Carlo samples does not exceed the available
  inherited samples;
- the shared plotting functions are on the MATLAB path.

---

## Terminology

```text
ALL
    Estimate parameters using every construct in the selected library.

L1O
    Leave one construct out, estimate from the remaining constructs, and use
    the excluded construct for predictive assessment.

Reduced model
    Estimate kappa^0 while fixing rho^0.

Full model
    Estimate both kappa^0 and rho^0.

Wells
    Use individual well trajectories.

Instances
    Use one mean trajectory per experiment.

Global
    Use one global mean trajectory per construct.
```
