# Lib24 ALL reduced model

## Overview

This directory contains the joint parameter-estimation workflow for all 24 constructs in Lib24 under the reduced RBS model.

The reduced model estimates the RBS intrinsic initiation capacity, `kappa0 = k0/sigma0`, while fixing the inverse scaling parameter, `rho0 = inv_sigma0`. The main workflow uses `rho0 = 0.02`. A complementary robustness workflow repeats the estimation for several fixed `rho0` values used to assess the dependence of the inferred parameters and objective value on this modelling choice.

This subfolder belongs to `Estimation_Pi`. A global README for `Estimation_Pi` will later describe how this workflow connects to the other library-specific estimation stages.

## Library and free parameters

Lib24 contains:

```text
2 plasmid backbones x 3 promoters x 4 RBSs = 24 constructs
```

The fitted parameter vector contains eight entries:

```text
3 promoter transcription rates:
    J23106, J23102, J23101

4 RBS intrinsic initiation capacities:
    B0030, B0032, J61100, J61101

1 pGreen copy-number multiplier
```

The pSC101 effective copy number is fixed to `N = 5`. The pGreen effective copy number is obtained by multiplying this reference value by the estimated copy-number factor.

## Canonical workflow: rho0 = 0.02

### 1. Estimate parameters

Run:

```matlab
Estimate_Lib24_ALL_reduced_model
```

The script uses BADS with multiple randomized initial points and minimizes `J5w_LogPI_Lib24_ALL_reduced`.

Default configuration:

```matlab
Use_mean = 'Wells';
num_runs = 54;
RBS_inv_sigma0 = 0.02;
delta = 0.2;
```

Output:

```text
Estimated_results/Results_BADS_Lib24_ALL_reduced_Wells.mat
```

### 2. Compile the result tensor

Run:

```matlab
Generate_Results_Lib24_ALL_reduced_model
```

The script combines the independent BADS runs with the processed Lib24 experimental tensor. It adds:

- raw, mean, and standard-deviation parameter estimates;
- global, experiment-level, and well-level synthesis-rate sensitivities;
- deterministic predictions over a growth-rate grid;
- Monte Carlo parameter samples;
- prediction quantiles and kernel-density summaries.

Distributed output:

```text
Results_Tensor_Lib24_ALL_reduced_model_Wells.mat
```

This tensor is included because it is used by downstream plotting, Jacobian analysis, and cross-library workflows.

### 3. Plot the canonical results

Run:

```matlab
Show_Results_Lib24_ALL_reduced_model
```

This script loads the distributed Wells tensor and calls:

```matlab
Plot_Results_Lib_reduced_model(...)
```

The plotting function produces parameter histograms and experimental-versus-predicted synthesis-rate trajectories.

## Robustness analysis across fixed rho0 values

The fixed values considered are:

```matlab
set_RBS_inv_sigma0 = [0.005, 0.01, 0.02, 0.03, 0.05];
```

These cases were used to evaluate the sensitivity of the objective and inferred parameters to the selected fixed value of `rho0`. The `rho0 = 0.02` case is the canonical reference.

### 1. Repeat estimation

Run:

```matlab
Estimate_Lib24_ALL_reduced_model_set_sigma
```

For each fixed `rho0`, the script writes:

```text
Estimated_results/Results_BADS_Lib24_ALL_reduced_<Use_mean>_<rho_tag>.mat
```

The optimization files are retained in the local `Estimated_results/` subdirectory. They are not all required in the compact public distribution because of their combined size.

### 2. Regenerate the rho0-specific tensors

Run:

```matlab
Generate_Results_Lib24_ALL_reduced_model_set_sigma
```

This produces one large compiled tensor per fixed value:

```text
Results_Tensor_Lib24_ALL_reduced_model_<Use_mean>_<rho_tag>.mat
```

These derived tensors are not all distributed because each occupies approximately 24 MB. They can be regenerated from the optimization outputs.

### 3. Inspect one fixed-rho0 case

Run:

```matlab
Show_Results_Lib24_ALL_reduced_model_set_sigma
```

Select the aggregation layer and `rho0` value in the configuration section. The script loads the corresponding tensor, calls `Plot_Results_Lib_reduced_model`, and plots optimizer convergence diagnostics.

### 4. Compare all fixed-rho0 cases

Two complementary scripts are available:

```matlab
Show_Compared_Results_Lib24_ALL_reduced_model_set_sigma
Show_Robustness_rho_Lib24_ALL_reduced_model
```

`Show_Compared_Results_Lib24_ALL_reduced_model_set_sigma` compares parameter estimates against objective values across all fixed `rho0` values.

`Show_Robustness_rho_Lib24_ALL_reduced_model` summarizes the relative objective variation, selected parameter distributions, and distribution overlap. It exports:

```text
Figures/FigS_RhoRobustness.pdf
```

## Objective function

`J5w_LogPI_Lib24_ALL_reduced.m` evaluates the library-scale discrepancy between predicted and experimental synthesis rates.

For each construct and selected data layer, it computes residuals in base-10 logarithmic space:

```text
r = log10(Pi_predicted) - log10(Pi_experimental)
```

It then applies a pseudo-Huber loss with scale `delta` and a smooth synthesis-rate-dependent weight. The objective is averaged across trajectories, wells or experiments, and the 24 constructs.

Supported aggregation layers:

```text
Global
Instances
Wells
```

## Input data and dependencies

Primary inputs:

```text
Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
Experimental_Data/ExpData_Tensor_lib30_micro.mat
Experimental_Data/ExpData_Tensor_lib24_micro.mat
```

The estimation script uses the Lib30 experimental tensor as the source from which the Lib24 subset is selected. The result-generation scripts use the dedicated Lib24 tensor.

Required software and functions include:

- MATLAB;
- BADS;
- `init_SynTwin` and `SynTwin_path`;
- `Get_synthesis_predictions_lite`;
- `Get_synthesis_predictions`;
- Parallel Computing Toolbox for `parfor`, unless replaced with `for`;
- Statistics and Machine Learning Toolbox for `fitdist` and related distribution operations;
- `rmvnrnd` from the truncated multivariate-normal implementation distributed with SynTwin.

## File inventory

```text
Estimate_Lib24_ALL_reduced_model.m
Estimate_Lib24_ALL_reduced_model_set_sigma.m
Generate_Results_Lib24_ALL_reduced_model.m
Generate_Results_Lib24_ALL_reduced_model_set_sigma.m
J5w_LogPI_Lib24_ALL_reduced.m
Plot_Results_Lib_reduced_model.m
Show_Results_Lib24_ALL_reduced_model.m
Show_Results_Lib24_ALL_reduced_model_set_sigma.m
Show_Compared_Results_Lib24_ALL_reduced_model_set_sigma.m
Show_Robustness_rho_Lib24_ALL_reduced_model.m
Estimated_results/
Results_Tensor_Lib24_ALL_reduced_model_Wells.mat
```

## Data-distribution policy

The compact distribution includes the canonical compiled Wells tensor at `rho0 = 0.02` and the scripts required to reproduce the other cases. The large rho0-specific tensors are treated as regenerable derived outputs and are therefore omitted to avoid duplicating several files of approximately 24 MB each.
