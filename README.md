# SynTwin — Host-aware digital twin for biopart characterization

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.18787107.svg)](https://doi.org/10.5281/zenodo.18787107)
![License](https://img.shields.io/badge/license-MIT-blue.svg)

> **Release status.** This repository contains the SynTwin implementation
> associated with the manuscript
> *Host-aware Identification of Intrinsic Gene Expression Biopart Parameters
> using Combinatorial Libraries*, which has received an in-principle publication
> decision from *Nature Communications* and is undergoing final editorial
> revision.

## Overview

SynTwin is a MATLAB framework for quantitative, host-aware characterization of
genetic bioparts embedded in combinatorial transcriptional-unit libraries.

The framework combines:

- processed time-resolved experimental measurements;
- a Host Equivalent Model (HEM) of *E. coli* physiology;
- a host-aware digital twin driven by the measured specific growth rate,
  $\mu(t)$;
- global parameter estimation with BADS;
- hierarchical parameter transfer between combinatorial libraries;
- deterministic and Monte Carlo synthesis-rate predictions;
- local practical-identifiability analysis;
- scripts that reproduce the main experimental and computational figures.

The inferred biopart quantities include:

```text
effective plasmid copy number
promoter transcription rate
RBS intrinsic initiation capacity
RBS inverse scaling parameter in the full model
```

The repository contains processed experimental data and stored estimation
outputs. Raw Cytation Excel files are not part of the public distribution.

---

## Repository structure

The current top-level layout is:

```text
SynTwin/
├── Estimation_Pi/
├── Experimental_data/
├── Jacobian_analysis/
├── Scripts_base/
├── Show_Data_Exp/
├── Show_Results/
├── Generate_HEM/
├── MATLAB/
├── doc/
├── init_SynTwin.m
├── SynTwin_path.m
├── SynTwin_root.m
├── CITATION.cff
├── Data Availability & Attribution Statement.md
├── LICENSE.md
└── README.md
```

### `Estimation_Pi/`

Host-aware parameter-estimation workflows for the background model and the
Lib24, Lib30, Lib6, and Lib5 combinatorial libraries.

Current subfolders:

```text
Estimation_Pi/
├── L5_ALL_reduced_model/
├── L5_L1O_reduced_model/
├── L6_ALL_reduced_model/
├── L6_ALL_reduced_model_dmA/
├── L6_L1O_reduced_model/
├── L24_ALL_reduced_model/
├── L24_ALL_reduced_model_dmA/
├── L24_L1O_reduced_model/
├── L30_ALL_full_model/
├── L30_ALL_reduced_model_dmA/
├── pARKA21_background_model/
└── README_Estimation_Pi.md
```

Each workflow contains estimation, result-generation, and visualization scripts,
together with its own README and, where included, the stored
`Estimated_results/` and compiled result tensors.

See:

```text
Estimation_Pi/README_Estimation_Pi.md
```

### `Experimental_data/`

Processed MATLAB experimental datasets used by the background, estimation, and
visualization workflows.

Current files:

```text
ExpData_pARKA21.mat
ExpData_Tensor_lib5_micro.mat
ExpData_Tensor_lib24_micro.mat
ExpData_Tensor_lib24.mat
ExpData_Tensor_lib30_micro.mat
ExpData_Tensor_lib30.mat
ExpData_Tensor_lib35.mat
```

Files ending in `_micro.mat` contain the minimum growth- and synthesis-rate
data required by parameter-estimation workflows. Files without `_micro` retain
the complete processed hierarchy, including particles, MEFL, growth,
synthesis-rate, burden, experiment, and well-level information.

See:

```text
Experimental_data/README_Experimental_data.md
```

### `Jacobian_analysis/`

Local practical-identifiability analyses based on experimental synthesis-rate
sensitivities.

The module includes:

- numerical rank and effective-rank comparisons;
- absolute and relative Jacobians;
- singular-value and condition-number analysis;
- least-identifiable parameter combinations;
- structural parameter–construct incidence plots.

See:

```text
Jacobian_analysis/README_Jacobian_analysis.md
```

### `Scripts_base/`

Shared MATLAB functions used across the repository.

These include:

- host-aware synthesis prediction;
- reduced and full model utilities;
- Cytation workbook parsing;
- experimental-data processing;
- shared plotting functions;
- model and tensor helper functions.

Runnable scripts add this folder through `init_SynTwin`.

### `Show_Data_Exp/`

Experimental-data visualization for the 35-construct constitutive-expression
collection.

The main workflow:

```matlab
Show_Data_Exp_lib35
```

loads:

```text
Experimental_data/ExpData_Tensor_lib35.mat
```

and generates the experimental summary and growth-resolved synthesis-rate
figures.

See the local README in `Show_Data_Exp/`.

### `Show_Results/`

Scripts that regenerate the main stored-result visualizations without rerunning
parameter estimation.

The module covers:

- host-physiology validation;
- plasmid, promoter, and RBS parameter distributions;
- host-aware RBS translation maps.

See:

```text
Show_Results/README_Show_Results.md
```

### `Generate_HEM/`

Host Equivalent Model resources and scripts.

This directory contains:

- the distributed HEM surrogate used by the digital twin;
- host-model preparation and validation scripts;
- third-party datasets used for HEM fitting;
- associated subworkflows for wild-type, nutrient, and mass relationships.

The parameter-estimation scripts generally load:

```text
Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
```

### `MATLAB/`

Third-party MATLAB software distributed or interfaced with by SynTwin,
including BADS and MEIGO/eSS where present.

These packages retain their original authorship and licenses.

### `doc/`

User and software documentation, including the SynTwin User Manual.

---

## Quick start

### 1. Initialize the repository

Start MATLAB and run:

```matlab
ROOT = init_SynTwin();
```

Use:

```matlab
ROOT = init_SynTwin('experimental',true);
```

for workflows requiring processed experimental data, and:

```matlab
ROOT = init_SynTwin('experimental',true,'bads',true);
```

for parameter-estimation workflows using BADS.

The initialization is portable and does not require hard-coded absolute paths.

### 2. Resolve repository files

Use `SynTwin_path`:

```matlab
data_file = SynTwin_path( ...
    'Experimental_data', ...
    'ExpData_Tensor_lib30_micro.mat');
```

Use the exact capitalization adopted by the distributed repository.

### 3. Run a workflow

A standard estimation subfolder follows:

```matlab
Estimate_<Library>_<Scheme>_<Model>
Generate_Results_<Library>_<Scheme>_<Model>
Show_Results_<Library>_<Scheme>_<Model>
```

For example:

```matlab
Estimate_Lib24_ALL_reduced_model
Generate_Results_Lib24_ALL_reduced_model
Show_Results_Lib24_ALL_reduced_model
```

Stored canonical result tensors can be visualized directly without rerunning
the optimizer.

---

## Main computational workflow

SynTwin separates experimental processing, background characterization,
parameter estimation, result compilation, and visualization.

```text
processed experimental data
          |
          v
pARKA21 background characterization
          |
          v
host-aware synthesis-rate data
          |
          v
hierarchical library parameter estimation
          |
          v
result tensors and Monte Carlo predictions
          |
          +--> manuscript figures
          |
          +--> Jacobian/identifiability analysis
```

### Background model

The pARKA21 workflow estimates:

```text
x = 100 mu
Pi_bk(mu) = b0 + A(1 - exp(-x/k)) + q1 x + q2 x^2
```

using BADS multistart and processed data from:

```text
Experimental_data/ExpData_pARKA21.mat
```

### Hierarchical library characterization

The principal reduced-model inference sequence is:

```text
Lib24
  |
  +--> Lib6
          |
          +--> Lib5
```

- **Lib24** identifies three promoter rates, four RBS intrinsic initiation
  capacities, and the pGreen copy-number multiplier.
- **Lib6** inherits promoter and plasmid distributions from Lib24 and identifies
  B0034.
- **Lib5** inherits plasmid and RBS distributions from Lib24 and Lib6 and
  identifies J23100.

Lib30 provides an independent full-model analysis and reduced-model robustness
analysis over the full 30-construct library.

---

## Estimation schemes and model variants

### `ALL`

Uses every construct in the selected library for estimation.

### `L1O`

Leaves one construct out per fold, estimates from the remaining constructs, and
retains the excluded construct for predictive assessment.

### Reduced model

Estimates:

```text
kappa^0 = k0/sigma0
```

while fixing:

```text
rho^0 = 1/sigma0
```

The canonical reduced-model workflows use:

```text
rho^0 = 0.02
```

### Full model

Estimates both `kappa^0` and `rho^0`.

### Experimental aggregation

The objective functions support:

```text
Global
Instances
Wells
```

The canonical distributed tensors use `Wells` unless stated otherwise.

---

## Objective functions

The library-estimation workflows use `J5w_*` objective functions.

For each prediction residual:

```text
r = log10(Pi_pred) - log10(Pi_exp)
```

the workflows apply a pseudo-Huber loss:

```text
rho(r) = delta^2 * (sqrt(1 + (r/delta)^2) - 1)
```

with synthesis-rate-dependent weighting:

```text
weight = sqrt(1 + log10(1 + Pi_exp))
```

The scalar cost is averaged over the selected trajectory aggregation and the
constructs included in the fit.

The pARKA21 background workflow instead uses an instance-balanced mean squared
residual.

---

## Experimental data and result tensors

### Experimental tensors

Micro tensors contain the minimum estimation inputs:

```text
TU metadata
mu(t)
Pi(t)
Global summaries
Instance summaries
Well trajectories
```

Complete tensors additionally contain:

```text
Particles
MEFL
MEFL/Particles
growth-rate trajectories
synthesis-rate trajectories
strain and host burden
experiment metadata
well-level measurements
```

### Result tensors

Compiled result tensors can contain:

```text
estimated parameter distributions
objective-function values
TU metadata
experimental trajectories
global, instance, and well sensitivities
deterministic synthesis predictions
Monte Carlo parameter samples
growth-resolved prediction intervals
```

L1O tensors additionally retain fold-specific local parameter distributions.

See the local READMEs and the User Manual for the complete schemas.

---

## Reproducing the main outputs

### Experimental characterization figures

```matlab
Show_Data_Exp_lib35
```

### Stored manuscript-result figures

```matlab
Plot_Host_Physiology_Model_Validation
Plot_Biopart_Parameter_Distributions
Plot_RBS_Host_Aware_Translation_Maps
```

### Parameter-estimation results

Run the corresponding `Estimate_*`, `Generate_Results_*`, and `Show_Results_*`
scripts under `Estimation_Pi/`.

### Practical-identifiability analysis

```matlab
Analyze_Libs_Jacobian_ranks
Analyze_Lib24_reduced_Jacobian
Plot_Combinatorial_Library_Jacobian_Structure
```

---

## Distribution policy

The repository includes:

- processed `.mat` experimental data;
- parameter-estimation scripts;
- objective functions;
- stored BADS estimation outputs where distributed;
- canonical compiled Wells tensors where practical;
- scripts for regeneration and visualization;
- HEM surrogate data;
- shared MATLAB utilities.

The repository does not distribute the raw Cytation Excel workbooks used to
produce the processed experimental tensors.

Large robustness result tensors are omitted when they can be regenerated from
the distributed estimation outputs.

Third-party datasets used for HEM fitting are documented separately.

---

## Requirements

Recommended environment:

- MATLAB R2020a or later;
- support for local functions in scripts;
- BADS for parameter estimation;
- Parallel Computing Toolbox for `parfor`;
- Statistics and Machine Learning Toolbox for distribution and percentile
  functions;
- `exportgraphics`;
- the truncated multivariate-normal sampler distributed with SynTwin;
- `violinplot` for selected manuscript figures.

Workflows using `parfor` can be changed to `for` when Parallel Computing
Toolbox is unavailable.

GNU Octave compatibility has not been tested.

---

## Reproducibility

BADS uses randomized initial points. Some hierarchical workflows also select
inherited Monte Carlo samples randomly.

For exact reruns, set a random seed before estimation:

```matlab
rng(1,'twister');
```

Before regenerating a result tensor, verify that:

- the selected `Use_mean` matches the estimation output;
- required upstream result tensors are available;
- the estimation filenames match the generator;
- sufficient inherited Monte Carlo samples exist;
- shared functions are on the MATLAB path.

Stored estimation outputs and canonical result tensors can be used directly
without repeating the optimizations.

---

## Documentation

Repository-level and module-level documentation includes:

```text
README.md
Estimation_Pi/README_Estimation_Pi.md
Experimental_data/README_Experimental_data.md
Jacobian_analysis/README_Jacobian_analysis.md
Show_Data_Exp/README_Lib35_Experimental_Data.md
Show_Results/README_Show_Results.md
doc/SynTwin_User_Manual.pdf
```

The local READMEs describe the current distributed interfaces. Where older
manual text differs from the present scripts, the current scripts and local
READMEs take precedence.

---

## Citation

Citation metadata are provided in:

```text
CITATION.cff
```

Software DOI:

```text
https://doi.org/10.5281/zenodo.18787107
```

Current software version in the citation metadata:

```text
1.1.0
```

**Associated manuscript:**  
Picó, J.; Arboleda-García, A.; Penas, D. R.; Banga, J. R.;
Vignoni, A.; Boada, Y.; Zach, P. (2026).  
*Host-aware Identification of Intrinsic Gene Expression Biopart Parameters
using Combinatorial Libraries.*  
Nature Communications, accepted in principle; final editorial revision in progress.

Use the repository citation metadata for the archived software version.

---

## Third-party datasets and software

SynTwin includes or interfaces with third-party datasets and software.

The repository-level file:

```text
Data Availability & Attribution Statement.md
```

documents:

- the Chure and Cremer Flux Parity datasets;
- the Zenodo archive `10.5281/zenodo.5893799`;
- upstream values from Bremer and Dennis;
- the role of these data in HEM fitting;
- attribution requirements;
- third-party optimization packages.

Third-party data and software retain their original ownership and license
conditions.

---

## License

SynTwin source code is distributed under the MIT License.

See:

```text
LICENSE.md
```

The license applies to SynTwin code and documentation. It does not replace the
licenses or terms of third-party datasets and software.

---

## Contact

For questions, bug reports, or attribution concerns, open an issue in the
repository or contact:

```text
jpico@upv.edu.es
```
