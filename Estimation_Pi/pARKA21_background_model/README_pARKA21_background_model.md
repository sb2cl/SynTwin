# pARKA21 background-model estimation

## Overview

This directory estimates the growth-dependent autofluorescence synthesis-rate
background measured with the nonexpressing construct pARKA21, also identified as
`PARKA` in some Cytation layouts.

The fitted model is

```text
x = 100 mu
Pi_bk(mu) = b0 + A (1 - exp(-x/k)) + q1 x + q2 x^2
```

with parameter order

```text
[b0, A, k, q1, q2]
```

The factor `100` rescales the growth rate before evaluating the model. All
scripts use this convention.

## Data-distribution policy

The distributed experimental input is

```text
Experimental_Data/ExpData_pARKA21.mat
```

containing

```text
Data_pARKA21
```

The stored synthesis-rate trajectories are generated with
`consider_BK = 0`; therefore, no previously fitted background model has been
subtracted.

## Distributed workflow

### 1. Estimate the model

```matlab
Estimate_pARKA21_background_model
```

The script loads

```text
Experimental_Data/ExpData_pARKA21.mat
```

and performs independent BADS runs in parallel.

The default configuration is

```text
num_runs = 100
aggregation = instance-balanced mean squared error
discarded leading points = 4
discarded trailing points = 3
```

The same temporal cropping is used during estimation and plotting.

Output:

```text
Estimated_results/
    Results_BADS_pARKA21_background_Wells.mat
```

The main variable is

```text
Results_BADS_pARKA21_background
```

and each cell stores

```text
results = [b0, A, k, q1, q2, Jmin]
```

### 2. Compile the result structure

```matlab
Generate_Results_pARKA21_background_model
```

Output:

```text
Results_pARKA21_background_Wells.mat
```

Saved variables:

```text
Results_pARKA21_background
Params_global_BKmodel_mean
```

`Params_global_BKmodel_mean` retains the backward-compatible parameter vector
used by downstream background correction.

The compiled structure contains

- all BADS parameter vectors and objective values;
- the best parameter vector;
- the mean of the five lowest-cost solutions;
- parameter means and standard deviations;
- the experimental `Data_pARKA21` structure;
- per-instance predictions, residuals, and RMSE values;
- a dense model curve over the experimental growth-rate range.

### 3. Plot the results

```matlab
Show_Results_pARKA21_background_model
```

The script exports

```text
Figures/
    pARKA21_background_experimental.png/.svg
    pARKA21_background_fit.png/.svg
    pARKA21_background_parameters.png/.svg
    pARKA21_background_residuals.png/.svg
```

The first figure displays the experimental trajectories before adding the
fitted model.

## Objective function

`J_pARKA21_background.m` computes the mean squared residual within each
experimental instance and then averages the instance-level costs. Consequently,
experiments with different trajectory lengths do not receive weights
proportional to their number of time points.

Invalid parameter vectors or predictions return `Inf`.

## Distributed contents

```text
pARKA21_background_model/
├── README.md
├── Estimate_pARKA21_background_model.m
├── Generate_Results_pARKA21_background_model.m
├── J_pARKA21_background.m
├── Show_Results_pARKA21_background_model.m
├── Estimated_results/
└── Results_pARKA21_background_Wells.mat
```

## MATLAB requirements

- BADS;
- Parallel Computing Toolbox for `parfor`, unless replaced with `for`;
- Statistics and Machine Learning Toolbox;
- `exportgraphics`;
- SynTwin `Scripts_base` data-processing functions.
