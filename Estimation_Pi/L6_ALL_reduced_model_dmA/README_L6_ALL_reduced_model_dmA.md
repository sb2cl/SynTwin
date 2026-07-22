# Lib6 ALL reduced model: robustness to \(d_{mA}\)

## Overview

This directory repeats the Lib6 all-construct reduced-model estimation while
fixing the non-ribosomal mRNA degradation rate at:

```text
0.15, 0.17, 0.20, 0.22, and 0.25 min^-1
```

Lib6 contains six constructs formed from two plasmid backbones, three promoters,
and the shared B0034 RBS.

Only the B0034 intrinsic initiation capacity, \(\kappa^0\), is estimated.
Promoter rates and effective plasmid copy numbers are inherited from the Lib24
leave-one-out reduced-model tensor. The RBS inverse scaling parameter is fixed
at:

```text
rho^0 = 0.02
```

The \(d_{mA}=0.20\ \mathrm{min}^{-1}\) case is the reference in the robustness
analyses.

## Distribution policy

The public distribution includes:

```text
Estimated_results/
```

with one BADS output file for each fixed \(d_{mA}\).

The generated result tensors are not distributed because each file is
approximately 6 MB:

```text
Results_Tensor_Lib6_ALL_reduced_model_Wells_dma_<tag>.mat
```

They can be regenerated with:

```matlab
Generate_Results_Lib6_ALL_reduced_model_set_dma
```

## Workflow

### 1. Estimate B0034 across \(d_{mA}\)

```matlab
Estimate_Lib6_ALL_reduced_model_set_dma
```

For each \(d_{mA}\), the script:

- loads the Lib24 leave-one-out Wells tensor;
- extracts inherited `Omega_MC_samples` and `Gene_cn_MC_samples`;
- verifies that `num_runs` does not exceed the available inherited samples;
- estimates the shared B0034 \(\kappa^0\) with BADS.

Outputs:

```text
Estimated_results/
    Results_BADS_Lib6_ALL_reduced_Wells_dma_<tag>.mat
```

The stored variable is:

```text
Results_BADS_B0034_ALL_reduced
```

### 2. Regenerate complete result tensors

```matlab
Generate_Results_Lib6_ALL_reduced_model_set_dma
```

The script generates one tensor per \(d_{mA}\):

```text
Results_Tensor_Lib6_ALL_reduced_model_Wells_dma_<tag>.mat
```

Each tensor contains:

- the B0034 estimate distribution and objective values;
- promoter and plasmid statistics inherited from Lib24;
- Lib6 experimental data selected from the Lib30 tensor;
- synthesis-rate sensitivities;
- deterministic predictions over the growth-rate grid;
- Monte Carlo prediction intervals.

These files are regeneration products and are intentionally omitted from the
distribution.

### 3. Plot each \(d_{mA}\) case

```matlab
Show_Results_Lib6_ALL_reduced_model_dma
```

The script loads all five regenerated tensors and calls:

```matlab
Plot_Results_Lib6_reduced_model(...)
```

One figure set is exported per \(d_{mA}\) to:

```text
Figures/
```

`Plot_Results_Lib6_reduced_model.m` is supplied by the sibling
`L6_ALL_reduced_model` workflow and must be on the MATLAB path.

### 4. Analyze Lib6 robustness

```matlab
Show_Robustness_dma_Lib6_ALL_reduced_model
```

This script analyzes only the parameter estimated in Lib6:

```text
B0034 kappa^0
```

It exports:

- relative objective changes;
- B0034 distributions by \(d_{mA}\);
- log-log sensitivity;
- balanced collapsed and reference distributions.

Outputs are written to:

```text
Figures/
Tables/
```

### 5. Combine Lib24 and Lib6 robustness

```matlab
Show_Robustness_dma_Lib24_plus_Lib6_ALL_reduced_model
```

This optional script combines:

- promoters, pGreen, B0030, B0032, J61100, and J61101 from
  `L24_ALL_reduced_model_dmA`;
- B0034 from this Lib6 workflow.

Before running it, regenerate the five result tensors in both sibling
workflows. The expected directory layout is:

```text
Estimation_Pi/
├── L24_ALL_reduced_model_dmA/
└── L6_ALL_reduced_model_dmA/
```

## Objective function

`J5w_LogPI_Lib6_ALL_reduced.m` evaluates a weighted pseudo-Huber loss in
base-10 logarithmic synthesis-rate space.

The selected \(d_{mA}\) enters through:

```matlab
model_c.dm_c
```

For run `num_run`, all six Lib6 constructs use the corresponding inherited
Lib24 samples of promoter rate and effective copy number.

## Required resources

```text
Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
Experimental_Data/ExpData_Tensor_lib30_micro.mat
Estimation_Pi/L24_L1O_reduced_model/
    Results_Tensor_Lib24_L1O_reduced_Wells.mat
L6_ALL_reduced_model/Plot_Results_Lib6_reduced_model.m
```

The combined robustness analysis additionally requires regenerated tensors from:

```text
L24_ALL_reduced_model_dmA/
```

## MATLAB requirements

- BADS;
- Parallel Computing Toolbox for `parfor`, unless replaced with `for`;
- Statistics and Machine Learning Toolbox;
- `exportgraphics`;
- the truncated multivariate-normal sampler distributed with SynTwin.

GNU Octave compatibility has not been tested.
