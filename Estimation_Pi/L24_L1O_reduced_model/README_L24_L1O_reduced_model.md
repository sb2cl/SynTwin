# Lib24 leave-one-out reduced model

## Overview

This directory performs leave-one-construct-out analysis for the Lib24 reduced
translation model.

Lib24 contains 24 constructs formed from:

- two plasmid backbones;
- three promoters;
- four RBS variants.

For each leave-one-out case, one construct is excluded and the remaining
23 constructs are used to estimate:

```text
3 promoter transcription rates
4 RBS intrinsic initiation capacities, kappa^0
1 pGreen/pSC101 copy-number multiplier
```

The RBS inverse scaling parameter is fixed at:

```text
rho^0 = 0.02
```

## Distributed contents

This subfolder distributes:

```text
Estimated_results/
Results_Tensor_Lib24_L1O_reduced_Wells.mat
```

The compiled Wells tensor can therefore be plotted directly. The corresponding
tensor can also be regenerated from the supplied estimation outputs.

## Workflow

### 1. Run the 24 leave-one-out estimations

```matlab
Estimate_Lib24_L1O_reduced_model
```

The script iterates over all Lib24 constructs. For each case, it excludes one
construct and performs repeated BADS optimizations using the other 23.

Generated files follow:

```text
Estimated_results/
    Results_BADS_Lib24_L1O_reduced_<plasmid><promoter><RBS>_<Use_mean>.mat
```

The three-digit construct tag stores the original Lib30 tensor indices used to
identify the Lib24 subset. Because Lib24 excludes B0034, the RBS digit can be:

```text
1, 2, 4, or 5
```

Each file contains the repeated optimization runs for one left-out construct.

### 2. Compile the leave-one-out results

```matlab
Generate_Results_Lib24_L1O_reduced_model
```

The generator loads all 24 estimation files and builds:

```text
Results_Tensor_Lib24_L1O_reduced_<Use_mean>.mat
```

The result structure contains two complementary parameter summaries:

- `Estimated_parameters.TU{p,q,r}` stores the parameter distribution obtained
  when construct `(p,q,r)` was left out;
- `Estimated_parameters.ALL_*` pools the runs from all 24 leave-one-out cases.

The compiled tensor also contains:

- processed experimental data;
- local parameter distributions for each left-out construct;
- pooled parameter distributions;
- global, instance-level, and well-level synthesis-rate sensitivities;
- growth-rate-dependent predictions;
- Monte Carlo prediction intervals.

### 3. Plot the distributed tensor

```matlab
Show_Results_Lib24_L1O_reduced_model
```

By default, the script uses:

```matlab
wells = true;
instances = false;
```

and loads:

```text
Results_Tensor_Lib24_L1O_reduced_Wells.mat
```

The figures are exported to:

```text
Figures/
```

The plotting call uses:

```matlab
Plot_Results_Lib_reduced_model( ...
    Results_Tensor_Lib24_L1O_reduced, ...
    Lower_Bounds, Upper_Bounds, ...
    indices_plasmids, indices_promoters, indices_rbss, title_text, ...
    Name, Value, ...)
```

This is compatible with the current plotting function, which retains the seven
original positional arguments and accepts name-value options.

## Objective function

`J5w_LogPI_Lib24_L1O_reduced.m` evaluates the weighted pseudo-Huber objective
over the 23 training constructs in each leave-one-out case.

The left-out construct is defined by local Lib24 tensor indices:

```matlab
Construct_2_LO = [i,j,k];
```

The corresponding experimental data are read from the Lib30 processed tensor,
using the Lib24 index subsets:

```matlab
indices_plasmids_lib24 = [1,2];
indices_promoters_lib24 = [1,2,3];
indices_rbss_lib24 = [1,2,4,5];
```

## Requirements

- MATLAB with local functions in scripts;
- BADS;
- Parallel Computing Toolbox for `parfor`, unless replaced with `for`;
- Statistics and Machine Learning Toolbox;
- `exportgraphics`;
- the truncated multivariate-normal sampler used by SynTwin;
- HEM surrogate data;
- processed Lib24 and Lib30 experimental tensors;
- `Plot_Results_Lib_reduced_model.m` on the MATLAB path.

GNU Octave compatibility has not been tested.
