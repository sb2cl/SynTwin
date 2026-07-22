# Lib6 ALL reduced model

## Overview

This directory estimates the shared B0034 RBS parameter using all six Lib6
constructs.

Lib6 contains:

```text
2 plasmid backbones
3 promoters
1 shared RBS: B0034
6 constructs
```

Only the B0034 intrinsic initiation capacity, \(\kappa^0\), is estimated in
this workflow. The promoter transcription rates and effective plasmid copy
numbers are inherited from the Lib24 leave-one-out reduced-model results.

The RBS inverse scaling parameter is fixed at:

```text
rho^0 = 0.02
```

## Parameter-transfer design

For each optimization run, the estimator selects the matching Monte Carlo
sample from:

```text
Results_Tensor_Lib24_L1O_reduced_Wells.mat
```

The inherited quantities are:

```text
Omega_MC_samples
Gene_cn_MC_samples
```

These values are held fixed while BADS estimates the shared B0034
\(\kappa^0\). Thus, variability across Lib6 runs includes uncertainty
propagated from the promoter and plasmid parameters characterized in Lib24.

## Distributed contents

This subfolder distributes:

```text
Estimated_results/
Results_Tensor_Lib6_ALL_reduced_Wells.mat
```

The compiled Wells tensor can be plotted directly and can also be regenerated
from the distributed estimation output and the Lib24 dependency.

## Workflow

### 1. Estimate B0034 \(\kappa^0\)

Run:

```matlab
Estimate_Lib6_ALL_reduced_model
```

The script:

- loads the HEM surrogate and Lib30 experimental tensor;
- loads the Lib24 leave-one-out Wells tensor;
- extracts inherited promoter and plasmid Monte Carlo samples;
- verifies that `num_runs` does not exceed the available inherited samples;
- runs one BADS estimation for each inherited sample.

Output:

```text
Estimated_results/
    Results_BADS_Lib6_ALL_reduced_Wells.mat
```

The stored variable is:

```text
Results_BADS_B0034_ALL_reduced
```

Each entry contains:

```text
[kappa0_B0034, objective_value]
```

### 2. Generate the Lib6 result tensor

Run:

```matlab
Generate_Results_Lib6_ALL_reduced_model
```

The script selects the six B0034 constructs from the Lib30-format experimental
tensor using:

```matlab
indices_plasmids_lib6 = [1,2];
indices_promoters_lib6 = [1,2,3];
indices_rbss_lib6 = 3;
```

It combines:

- the B0034 estimates obtained in Lib6;
- promoter and plasmid statistics inherited from Lib24;
- experimental growth and synthesis-rate measurements;
- local sensitivity calculations;
- deterministic growth-rate-dependent predictions;
- Monte Carlo prediction intervals.

Output:

```text
Results_Tensor_Lib6_ALL_reduced_Wells.mat
```

The stored variable is:

```text
Results_Tensor_Lib6_ALL_reduced
```

### 3. Plot the distributed results

Run:

```matlab
Show_Results_Lib6_ALL_reduced_model
```

The script calls:

```matlab
Plot_Results_Lib6_reduced_model(...)
```

and writes PNG and SVG files to:

```text
Figures/
```

The plotting function produces:

- the B0034 \(\kappa^0\) histogram;
- a \(2\times3\) tiled comparison of experimental and predicted
  synthesis-rate trajectories.

## Plotting compatibility

`Plot_Results_Lib6_reduced_model.m` retains seven positional arguments:

```matlab
Plot_Results_Lib6_reduced_model( ...
    Results_Lib, Lower_Bounds, Upper_Bounds, ...
    indices_plasmids, indices_promoters, indices_rbss, title_text, ...
    Name, Value, ...)
```

The revised `Show_Results` call is compatible with this interface and supplies
the output directory, filename prefix, font sizes, and figure dimensions as
name-value options.

The expected result tensor has dimensions:

```text
2 x 3 x 1
```

where the singleton third dimension corresponds to B0034.

## Objective function

`J5w_LogPI_Lib6_ALL_reduced.m` computes a weighted pseudo-Huber loss in
base-10 logarithmic synthesis-rate space.

For run `num_run`, every Lib6 construct uses the corresponding inherited
Lib24 samples:

```matlab
Omega_MC_samples(num_run)
Gene_cn_MC_samples(num_run)
```

The same B0034 \(\kappa^0\) is fitted across all six constructs.

## Required resources

```text
Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
Experimental_Data/ExpData_Tensor_lib30_micro.mat
Estimation_Pi/L24_L1O_reduced_model/
    Results_Tensor_Lib24_L1O_reduced_Wells.mat
Scripts_base/Get_synthesis_predictions.m
Scripts_base/Get_synthesis_predictions_lite.m
```

## MATLAB requirements

- BADS;
- Parallel Computing Toolbox for `parfor`, unless replaced with `for`;
- Statistics and Machine Learning Toolbox;
- `exportgraphics`;
- the truncated multivariate-normal sampler distributed with SynTwin.

GNU Octave compatibility has not been tested.
