# Lib30 ALL complete model

## Overview

This directory estimates the complete host-aware translation model using all
30 constructs in Lib30.

The complete model estimates two RBS parameters:

```text
kappa^0 = K^0/sigma^0
rho^0   = 1/sigma^0
```

Lib30 contains two plasmid backbones, three promoters, and five RBS variants.
Because the pSC101 copy number is fixed, the model has 14 free parameters:

```text
3 promoter transcription rates
5 kappa^0 parameters
5 rho^0 parameters
1 pGreen/pSC101 copy-number multiplier
```

## Distributed contents

This subfolder distributes:

```text
Estimated_results/
Results_Tensor_Lib30_ALL_full_model_Wells.mat
```

The Wells tensor can therefore be plotted directly without rerunning the
estimation or result-generation stages.

## Workflow

### 1. Estimate the complete model

```matlab
Estimate_Lib30_ALL_full_model
```

The script performs repeated BADS optimizations using all 30 constructs and
the objective function:

```text
J5w_LogPI_Lib30_ALL_full_model.m
```

The output is:

```text
Estimated_results/
    Results_BADS_Lib30_ALL_complete_<Use_mean>.mat
```

Each run stores the 14 estimated parameters, the objective value, and the
initial parameter vector.

### 2. Generate the result tensor

```matlab
Generate_Results_Lib30_ALL_full_model
```

The script loads the matching estimation file and generates:

```text
Results_Tensor_Lib30_ALL_full_model_<Use_mean>.mat
```

The tensor contains:

- processed experimental data;
- raw and summary parameter estimates;
- promoter, plasmid, \(\kappa^0\), and \(\rho^0\) samples;
- global, instance-level, and well-level sensitivities;
- growth-rate-dependent synthesis predictions;
- Monte Carlo parameter samples;
- propagated prediction intervals.

The distributed folder includes the Wells tensor.

### 3. Plot the distributed results

```matlab
Show_Results_Lib30_ALL_full_model
```

The script loads:

```text
Results_Tensor_Lib30_ALL_full_model_Wells.mat
```

and calls:

```matlab
Plot_Results_Lib_complete_model(...)
```

The call uses seven positional arguments followed by name-value options for
the output directory, filename prefix, export control, font sizes, and figure
dimensions.

Generated figures are written to:

```text
Figures/
```

## Plotting compatibility

The current `Show_Results` script expects the plotting function to retain:

```matlab
Plot_Results_Lib_complete_model( ...
    Results_Lib, Lower_Bounds, Upper_Bounds, ...
    indices_plasmids, indices_promoters, indices_rbss, title_text, ...
    Name, Value, ...)
```

The plotting function itself was not present among the files available for
this review, so its internal `inputParser` could not be verified directly.
The call is internally consistent with the interface used in the supplied
`Show_Results` script.

## Requirements

- MATLAB with local functions in scripts;
- BADS;
- Parallel Computing Toolbox for `parfor`, unless replaced with `for`;
- Statistics and Machine Learning Toolbox;
- `exportgraphics`;
- the truncated multivariate-normal sampler used by SynTwin;
- HEM surrogate data;
- processed Lib30 experimental data;
- `Plot_Results_Lib_complete_model.m` on the MATLAB path.

GNU Octave compatibility has not been tested.
