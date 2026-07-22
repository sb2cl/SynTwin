# Lib35 experimental-data visualization

## Overview

This folder contains the MATLAB workflow used to visualize the processed experimental measurements for the 35-construct constitutive-expression library (Lib35).

The main script is:

```matlab
Show_Data_Exp_lib35
```

It reproduces the experimental-data panel reported as Fig. 2a in:

> *Host-aware Identification of Intrinsic Gene Expression Biopart Parameters using Combinatorial Libraries*

The workflow also generates additional summaries of expression, growth, host burden, variability, and growth-resolved synthesis rates across the library.

## Folder contents

```text
Show_Data_Exp_lib35.m    Main script for loading the Lib35 dataset and generating the figures
Figures/                 Output directory created or populated by the plotting workflow
```

The script calls the shared plotting routine:

```text
Scripts_base/Plot_ExpData_Tensor_v10.m
```

## Input data

The workflow requires:

```text
Experimental_Data/ExpData_Tensor_lib35.mat
```

The MAT-file must contain:

```matlab
Data_Exp_Tensor_lib35
```

`Data_Exp_Tensor_lib35` is a `2 x 4 x 5` MATLAB cell array indexed by plasmid backbone, promoter, and RBS. Empty cells represent combinations that are not included in Lib35.

A construct entry can be accessed as follows:

```matlab
S = load(SynTwin_path('Experimental_Data', ...
    'ExpData_Tensor_lib35.mat'));

EC = S.Data_Exp_Tensor_lib35{1,1,1};
D  = EC.Data;
P  = D.Stats_Particles;
```

A typical top-level construct entry contains:

```text
TU_Ori
TU_Promoter
TU_RBS
TU_Bioparts
TU_Name
TU_color_code
Data
```

The `Data` field stores the processed observables and summary statistics, including:

```text
Construct_name
Global_stats_Mumax_Pi_p
Global_stats_Mumax_Phi_h
Global_stats_Mumax_Phi_s
Stats_Particles
Stats_MEFL
Stats_MEFL_Particles
Stats_Mu
Stats_Pi_p
Stats_Phi_s
Stats_Phi_h
Global_data_time
```

Nested structures such as `Stats_Particles`, `Stats_Mu`, and `Stats_Pi_p` retain experiment-level and well-level measurements in addition to global means, standard deviations, and trajectories. The MAT-file therefore contains the processed experimental dataset used by the plotting workflow, rather than only the values displayed in Fig. 2a.

## Generated figures

The script generates the following outputs from the same Lib35 dataset:

- the Fig. 2a variant showing ordered construct synthesis rate, maximum growth rate, and estimated construct-induced host burden;
- ordered synthesis-rate and growth-rate summaries;
- ordered synthesis-rate and host-burden summaries with construct labels;
- Fano-factor and coefficient-of-variation summaries;
- tiled growth-resolved synthesis-rate trajectories, $\Pi_p(\mu)$, for all expression constructs.

Output files are written to:

```text
Figures/
```

with the filename prefix:

```text
Lib35_Experimental
```

## Usage

Run the script from MATLAB:

```matlab
Show_Data_Exp_lib35
```

The script initializes the SynTwin project with:

```matlab
ROOT = init_SynTwin('experimental', true);
```

It can therefore be launched from any working directory, provided that the SynTwin initialization utilities and repository structure are available.

## Requirements

- MATLAB R2020a or later;
- the SynTwin initialization and path-management functions;
- `Scripts_base/Plot_ExpData_Tensor_v10.m`;
- no non-core MATLAB toolbox is required by `Show_Data_Exp_lib35.m`.

GNU Octave compatibility has not been tested.

## Warning handling

The script suppresses three known MATLAB warnings associated with:

- small negative values displayed on logarithmic axes;
- negligible imaginary components passed to plotting functions;
- MATLAB print/export content-type recommendations.

These warnings arise during figure generation and do not modify the underlying experimental data.

## Data provenance and archived release

The processed Lib35 dataset contains the experimental measurements and summary statistics used to generate the reported plots.

- GitHub: <https://github.com/sb2cl/SynTwin>
- Zenodo: <https://doi.org/10.5281/zenodo.18787107>

The complete data schema is documented in Appendix B of the **SynTwin User Manual** distributed with the archived release.
