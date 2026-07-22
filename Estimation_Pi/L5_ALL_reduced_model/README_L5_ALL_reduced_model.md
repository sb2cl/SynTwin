# Lib5 ALL reduced model

## Overview

This directory estimates the J23100 promoter transcription parameter using all
five Lib5 constructs.

Lib5 contains:

```text
1 plasmid backbone: pGreen
1 promoter: J23100
5 RBS variants
5 constructs
```

Only the promoter parameter `Omega` is estimated. Effective copy number and RBS
intrinsic initiation capacities are inherited from previous workflows.

The RBS inverse scaling parameter is fixed at:

```text
rho^0 = 0.02
```

## Parameter inheritance

The workflow uses two upstream result tensors:

```text
L24_L1O_reduced_model/Results_Tensor_Lib24_L1O_reduced_Wells.mat
L6_L1O_reduced_model/Results_Tensor_Lib6_L1O_reduced_Wells.mat
```

The inherited sources are:

- B0030, B0032, J61100, and J61101 from Lib24;
- B0034 from Lib6;
- pGreen effective copy-number Monte Carlo samples from Lib6, which itself
  carries the Lib24-derived plasmid characterization.

Each optimization run uses one consistent inherited realization across all
five constructs.

## Distributed contents

This subfolder distributes:

```text
Estimated_results/
Results_Tensor_Lib5_ALL_reduced_Wells.mat
```

The plotting function is shared from:

```text
Scripts_base/Plot_Results_Lib5_reduced_model.m
```

## Workflow

### 1. Estimate J23100 `Omega`

```matlab
Estimate_Lib5_ALL_reduced_model
```

The estimator determines the common sample pool available across inherited
copy-number and RBS distributions. It then samples inherited realizations with
replacement and performs one BADS optimization per realization.

Output:

```text
Estimated_results/Results_BADS_Lib5_ALL_reduced_Wells.mat
```

The stored variable is:

```text
Results_BADS_J23100_ALL_reduced
```

Each entry contains:

```text
[Omega_J23100, objective_value]
```

### 2. Generate the result tensor

```matlab
Generate_Results_Lib5_ALL_reduced_model
```

The generator combines:

- the estimated J23100 `Omega` distribution;
- inherited pGreen copy-number statistics;
- inherited RBS parameter distributions;
- Lib5 experimental data;
- synthesis-rate sensitivities;
- deterministic predictions;
- Monte Carlo prediction intervals.

Output:

```text
Results_Tensor_Lib5_ALL_reduced_Wells.mat
```

The saved variables are:

```text
Results_Tensor_Lib5_ALL_reduced
Estimated_omega_J23100
```

`Estimated_omega_J23100` includes the raw objective values in:

```text
Estimated_omega_J23100.J_raw
```

### 3. Plot the distributed tensor

```matlab
Show_Results_Lib5_ALL_reduced_model
```

The script validates the loaded tensor and calls:

```matlab
Plot_Results_Lib5_reduced_model(...)
```

The shared plotter expects a `1 x 1 x 5` tensor and produces:

- the J23100 `Omega` histogram;
- five experimental-versus-predicted synthesis-rate panels.

Figures are exported to:

```text
Figures/
```

## Objective function

`J5w_LogPI_Lib5_ALL_reduced.m` evaluates a weighted pseudo-Huber loss in
base-10 logarithmic synthesis-rate space over all five constructs.

The same estimated `Omega` is used for all five constructs, while each RBS uses
its inherited `kappa^0` sample and pGreen uses the inherited copy-number sample
selected for the current run.

## Required resources

```text
Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
Experimental_Data/ExpData_Tensor_lib5_micro.mat
Estimation_Pi/L24_L1O_reduced_model/
    Results_Tensor_Lib24_L1O_reduced_Wells.mat
Estimation_Pi/L6_L1O_reduced_model/
    Results_Tensor_Lib6_L1O_reduced_Wells.mat
Scripts_base/Plot_Results_Lib5_reduced_model.m
```

## MATLAB requirements

- BADS;
- Parallel Computing Toolbox for `parfor`, unless replaced with `for`;
- Statistics and Machine Learning Toolbox;
- `exportgraphics`;
- the truncated multivariate-normal sampler distributed with SynTwin.

GNU Octave compatibility has not been tested.
