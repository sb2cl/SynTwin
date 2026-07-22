# Lib6 leave-one-out reduced model

## Overview

This directory performs leave-one-construct-out analysis for the Lib6 reduced
translation model.

Lib6 contains:

```text
2 plasmid backbones
3 promoters
1 shared RBS: B0034
6 constructs
```

For each fold, one construct is excluded and the shared B0034 intrinsic
initiation capacity, \(\kappa^0\), is estimated from the remaining five
constructs.

Promoter transcription rates and effective plasmid copy numbers are inherited
from the Lib24 leave-one-out reduced-model tensor. The RBS inverse scaling
parameter is fixed at:

```text
rho^0 = 0.02
```

## Distributed contents

This subfolder distributes:

```text
Estimated_results/
Results_Tensor_Lib6_L1O_reduced_Wells.mat
```

The compiled Wells tensor is approximately 6 MB and can be plotted directly.
It can also be regenerated from the supplied estimation files.

## Workflow

### 1. Run the six leave-one-out estimations

```matlab
Estimate_Lib6_L1O_reduced_model
```

For each left-out construct, the script performs repeated BADS optimizations.

Run `num_run` uses the corresponding inherited Lib24 samples:

```text
Omega_MC_samples(num_run)
Gene_cn_MC_samples(num_run)
```

The estimator verifies that `num_runs` does not exceed the number of available
inherited samples.

Outputs follow:

```text
Estimated_results/
    Results_BADS_Lib6_L1O_reduced_<plasmid><promoter><RBS>_<Use_mean>.mat
```

Because Lib6 contains only B0034, the RBS digit in every filename is `3`.

Each file stores:

```text
Results_BADS_B0034_LOOCV_reduced{num_run}.results
```

with:

```text
[kappa0_B0034, objective_value]
```

### 2. Compile the result tensor

```matlab
Generate_Results_Lib6_L1O_reduced_model
```

The generator loads the six fold-specific files and creates:

```text
Results_Tensor_Lib6_L1O_reduced_Wells.mat
```

The saved variables are:

```text
Results_Tensor_Lib6_L1O_reduced
Estimated_parameters
```

The tensor has dimensions:

```text
2 x 3 x 1
```

For each left-out construct, local statistics are stored in:

```text
Parameters_local_raw
Parameters_local_mean
Parameters_local_std
J_local_raw
```

These values correspond to training on the other five constructs.

The generator also pools all folds:

```text
Estimated_parameters.ALL_raw
Estimated_parameters.ALL_mean
Estimated_parameters.ALL_std
Estimated_parameters.J_raw
```

The pooled B0034 fields in the result tensor summarize all L1O folds. They are
not an independent fit using all six constructs.

The compiled tensor additionally contains:

- promoter and plasmid distributions inherited from Lib24;
- Lib6 experimental data selected from the Lib30 tensor;
- global, instance-level, and well-level synthesis-rate sensitivities;
- deterministic predictions over a growth-rate grid;
- Monte Carlo prediction intervals.

### 3. Plot the distributed tensor

```matlab
Show_Results_Lib6_L1O_reduced_model
```

The script loads:

```text
Results_Tensor_Lib6_L1O_reduced_Wells.mat
```

and calls:

```matlab
Plot_Results_Lib6_reduced_model(...)
```

Figures are exported to:

```text
Figures/
```

The plotting function retains seven positional arguments followed by
name-value options. The current `Show_Results` call is compatible with this
interface.

## Objective function

`J5w_LogPI_Lib6_L1O_reduced.m` computes a weighted pseudo-Huber loss in
base-10 logarithmic synthesis-rate space.

The left-out construct is represented by local Lib6 tensor indices:

```matlab
Construct_2_LO = [i,j,k];
```

Here, `k=1` because Lib6 has one RBS dimension. Experimental data are retrieved
from RBS index `3` of the Lib30 tensor, corresponding to B0034.

The objective is averaged over the five training constructs.

## Required resources

```text
Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
Experimental_Data/ExpData_Tensor_lib30_micro.mat
Estimation_Pi/L24_L1O_reduced_model/
    Results_Tensor_Lib24_L1O_reduced_Wells.mat
L6_ALL_reduced_model/
    Plot_Results_Lib6_reduced_model.m
```

## MATLAB requirements

- BADS;
- Parallel Computing Toolbox for `parfor`, unless replaced with `for`;
- Statistics and Machine Learning Toolbox;
- `exportgraphics`;
- the truncated multivariate-normal sampler distributed with SynTwin.

GNU Octave compatibility has not been tested.
