# Lib5 leave-one-out reduced model

## Overview

This directory performs leave-one-construct-out analysis for the Lib5 reduced
model.

Lib5 contains:

```text
1 plasmid backbone: pGreen
1 promoter: J23100
5 RBS contexts
5 constructs
```

For each fold, one RBS-defined construct is excluded and the J23100 promoter
parameter, \(\Omega\), is estimated from the remaining four constructs.

The effective pGreen copy number and RBS intrinsic initiation capacities are
inherited from earlier Lib24 and Lib6 leave-one-out workflows. The RBS inverse
scaling parameter is fixed at:

```text
rho^0 = 0.02
```

## Distributed contents

This subfolder distributes:

```text
Estimated_results/
Results_Tensor_Lib5_L1O_reduced_Wells.mat
```

The compiled Wells tensor can be plotted directly and can also be regenerated
from the supplied fold-specific estimation files.

## Parameter inheritance

The inherited quantities are:

```text
Gene_cn
RBS kappa^0 for five RBSs
```

Sources:

```text
Lib24 L1O reduced:
    B0030, B0032, J61100, and J61101

Lib6 L1O reduced:
    B0034
    pGreen effective copy-number samples used during estimation
```

Each optimization run uses one shared Monte Carlo index across all included
constructs. This preserves the joint inherited realization within that run.

## Workflow

### 1. Run the five leave-one-out estimations

```matlab
Estimate_Lib5_L1O_reduced_model
```

The estimator:

- loads the Lib24 and Lib6 inherited parameter tensors;
- determines the common number of available inherited samples;
- randomly selects one inherited realization per BADS run;
- estimates J23100 \(\Omega\) using four training constructs;
- generates one result file per left-out RBS context.

Outputs:

```text
Estimated_results/
    Results_BADS_Lib5_L1O_reduced_141_Wells.mat
    Results_BADS_Lib5_L1O_reduced_142_Wells.mat
    Results_BADS_Lib5_L1O_reduced_143_Wells.mat
    Results_BADS_Lib5_L1O_reduced_144_Wells.mat
    Results_BADS_Lib5_L1O_reduced_145_Wells.mat
```

The code `14<RBS>` records the fixed pGreen/J23100 context and the left-out RBS
index.

Each file stores:

```text
Results_BADS_J23100_L1O_reduced{num_run}.results
```

with:

```text
[Omega_J23100, objective_value]
```

### 2. Compile the result tensor

```matlab
Generate_Results_Lib5_L1O_reduced_model
```

Output:

```text
Results_Tensor_Lib5_L1O_reduced_Wells.mat
```

Saved variables:

```text
Results_Tensor_Lib5_L1O_reduced
Estimated_parameters_J23100
```

The result tensor has dimensions:

```text
1 x 1 x 5
```

Local fold-specific fields:

```text
Parameters_local_raw
Parameters_local_mean
Parameters_local_std
J_local_raw
```

Pooled fields:

```text
Estimated_parameters_J23100.ALL_raw
Estimated_parameters_J23100.ALL_mean
Estimated_parameters_J23100.ALL_std
Estimated_parameters_J23100.J_raw
```

The pooled \(\Omega\) distribution combines all five leave-one-out folds. It is
not an independent fit using all five constructs.

The tensor also contains:

- inherited pGreen and RBS parameter statistics;
- experimental growth and synthesis-rate trajectories;
- sensitivities at experimental growth rates;
- deterministic synthesis-rate predictions;
- Monte Carlo prediction intervals and kernel-density summaries.

### 3. Plot the distributed tensor

```matlab
Show_Results_Lib5_L1O_reduced_model
```

The script validates and loads:

```text
Results_Tensor_Lib5_L1O_reduced_Wells.mat
```

and calls:

```matlab
Plot_Results_Lib5_reduced_model(...)
```

The plotting function is located in:

```text
Scripts_base/
```

Figures are written to:

```text
Figures/
```

## Objective function

`J5w_LogPI_Lib5_L1O_reduced.m` computes a weighted pseudo-Huber loss in
base-10 logarithmic synthesis-rate space.

For each fold:

```text
5 total constructs
1 left-out construct
4 training constructs
```

The objective is therefore normalized by four training constructs.

`Construct_2_LO` is the local Lib5 RBS index from 1 to 5.

## Required resources

```text
Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
Experimental_Data/ExpData_Tensor_lib5_micro.mat
Estimation_Pi/L24_L1O_reduced_model/
    Results_Tensor_Lib24_L1O_reduced_Wells.mat
Estimation_Pi/L6_L1O_reduced_model/
    Results_Tensor_Lib6_L1O_reduced_Wells.mat
Scripts_base/
    Plot_Results_Lib5_reduced_model.m
```

## MATLAB requirements

- BADS;
- Parallel Computing Toolbox for `parfor`, unless replaced with `for`;
- Statistics and Machine Learning Toolbox;
- `exportgraphics`;
- the truncated multivariate-normal sampler distributed with SynTwin.

GNU Octave compatibility has not been tested.
