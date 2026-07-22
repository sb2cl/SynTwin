# SynTwin processed experimental data

## Overview

`Experimental_data/` contains the processed MATLAB datasets used by SynTwin for
background characterization, host-aware parameter estimation, experimental-data
visualization, and downstream analyses.

The public distribution provides processed `.mat` files rather than the raw
Cytation Excel workbooks. This keeps the distributed inputs compact and gives
all estimation workflows a stable MATLAB interface.

The current directory contains:

```text
Experimental_data/
├── ExpData_pARKA21.mat
├── ExpData_Tensor_lib5_micro.mat
├── ExpData_Tensor_lib24_micro.mat
├── ExpData_Tensor_lib24.mat
├── ExpData_Tensor_lib30_micro.mat
├── ExpData_Tensor_lib30.mat
└── ExpData_Tensor_lib35.mat
```

Raw Cytation spreadsheets used to regenerate selected processed datasets are
maintainer inputs and are not part of the public distribution.

---

## Important path note

Several SynTwin scripts refer to this directory as:

```text
Experimental_Data/
```

whereas the folder may appear as:

```text
Experimental_data/
```

on some local copies.

Windows and default macOS file systems usually ignore this capitalization
difference, but Linux and other case-sensitive file systems do not. The
repository should use one spelling consistently. If the distributed folder is
named `Experimental_data`, update `init_SynTwin`, `SynTwin_path`, and any direct
path references accordingly. If the canonical repository spelling is
`Experimental_Data`, rename the folder before release.

---

## Data-distribution policy

The distributed files are processed experimental datasets.

They contain:

- transcriptional-unit metadata;
- specific growth-rate measurements, \(\mu(t)\);
- synthesis-rate measurements, \(\Pi(t)\);
- experiment-level and well-level summaries;
- complete processed particle, fluorescence, burden, and metadata structures
  where applicable.

The raw Excel files are not required to reproduce the parameter-estimation
workflows.

---

# File inventory

## `ExpData_pARKA21.mat`

### Purpose

Processed experimental input for estimation of the growth-dependent
autofluorescence synthesis-rate background using the nonexpressing pARKA21
construct.

Some Cytation layouts identify this construct as:

```text
PARKA
```

The distributed file contains:

```matlab
Data_pARKA21
```

The structure is generated with background correction disabled:

```matlab
consider_BK = 0;
```

Therefore, the stored synthesis-rate trajectories have not been corrected using
a previously fitted background model.

### Used by

```text
Estimation_Pi/pARKA21_background_model/
```

Primary scripts:

```matlab
Estimate_pARKA21_background_model
Generate_Results_pARKA21_background_model
Show_Results_pARKA21_background_model
```

---

## `ExpData_Tensor_lib5_micro.mat`

### Purpose

Minimum processed experimental tensor for the five-construct Lib5 library.

Lib5 contains:

```text
1 plasmid: pGreen
1 promoter: J23100
5 RBS contexts
5 transcriptional units
```

This file provides the \(\mu(t)\) and \(\Pi(t)\) trajectories used to estimate
the J23100 promoter parameter in:

```text
L5_ALL_reduced_model
L5_L1O_reduced_model
```

### Expected tensor dimensions

```text
1 x 1 x 5
```

### Expected variable

The estimation scripts expect:

```matlab
ExpData_Tensor_lib5_micro
```

---

## `ExpData_Tensor_lib24_micro.mat`

### Purpose

Minimum processed experimental tensor for the 24-construct Lib24 library.

Lib24 contains:

```text
2 plasmids
3 promoters
4 RBSs
24 transcriptional units
```

The RBS set excludes B0034.

### Expected tensor dimensions

```text
2 x 3 x 4
```

### Expected variable

```matlab
ExpData_Tensor_lib24_micro
```

### Typical use

This dedicated Lib24 tensor is used during result compilation and
library-specific analyses. Some estimation scripts instead select the Lib24
subset directly from the Lib30-format micro tensor to preserve the original
index convention.

---

## `ExpData_Tensor_lib24.mat`

### Purpose

Complete processed experimental tensor for Lib24.

Unlike the `_micro` version, this file retains the complete processed
experimental hierarchy, including particle and fluorescence trajectories,
derived growth and synthesis rates, burden estimates, global summaries,
experiment metadata, and well-level values.

### Typical use

- detailed experimental inspection;
- experimental-data plotting;
- validation of processed trajectories;
- analyses requiring quantities not included in the compact estimation tensor.

This file is not required by every parameter-estimation script.

---

## `ExpData_Tensor_lib30_micro.mat`

### Purpose

Minimum processed experimental tensor for the 30-construct Lib30 library.

Lib30 contains:

```text
2 plasmids
3 promoters
5 RBSs
30 transcriptional units
```

### Expected tensor dimensions

```text
2 x 3 x 5
```

### Expected variable

```matlab
ExpData_Tensor_lib30_micro
```

### Typical use

This is the principal experimental input for:

- Lib30 full-model estimation;
- Lib30 reduced-model robustness analysis;
- Lib24 estimation scripts that select the Lib24 subset internally;
- Lib6 estimation scripts that select the B0034 subset internally.

The Lib24 subset uses the Lib30 RBS indices:

```matlab
[1,2,4,5]
```

The Lib6 subset uses:

```matlab
RBS index 3
```

which corresponds to B0034.

---

## `ExpData_Tensor_lib30.mat`

### Purpose

Complete processed experimental tensor for Lib30.

It contains the same 30 transcriptional-unit positions as the Lib30 micro
tensor but retains the full processed experimental structures.

### Typical use

- complete experimental-data analysis;
- visualization beyond the minimum estimation inputs;
- examination of particle, MEFL, MEFL/Particles, growth, synthesis, and burden
  trajectories;
- access to experiment-level metadata.

---

## `ExpData_Tensor_lib35.mat`

### Purpose

Complete processed experimental tensor for the 35 constitutive expression
constructs used in the global experimental characterization figures.

The file contains:

```matlab
Data_Exp_Tensor_lib35
```

The tensor is stored in a:

```text
2 x 4 x 5
```

cell array indexed by plasmid, promoter, and RBS. Empty cells correspond to
combinations that are not part of the 35-construct collection.

### Used by

```matlab
Show_Data_Exp_lib35
```

and the shared plotting routine:

```matlab
Plot_ExpData_Tensor_v10
```

The file supports:

- the ordered synthesis-rate, maximum-growth-rate, and host-burden summary;
- supplementary synthesis/growth and synthesis/burden plots;
- variability summaries;
- growth-resolved \(\Pi(\mu)\) trajectories.

---

# Micro experimental tensors

## Purpose

Files ending in:

```text
_micro.mat
```

contain the minimum processed information required by the host-aware
parameter-estimation workflows.

The micro tensors are three-dimensional MATLAB cell arrays:

```text
num_plasmids x num_promoters x num_rbss
```

Each nonempty cell represents one transcriptional unit defined by a unique
combination of:

```text
origin
promoter
RBS
```

The compact schema is shared by Lib5, Lib24, and Lib30. Only the tensor
dimensions differ.

---

## Transcriptional-unit metadata

A typical cell contains:

```text
TU_Ori
TU_Promoter
TU_RBS
TU_Bioparts
TU_Name
TU_color_code
```

### `TU_Ori`

Origin or plasmid identifier, such as `pGreen` or `pSC101`.

### `TU_Promoter`

Promoter identifier.

### `TU_RBS`

RBS identifier.

### `TU_Bioparts`

Combined description of the origin, promoter, and RBS.

### `TU_Name`

Internal construct label.

### `TU_color_code`

Color used by the plotting workflows.

---

## Characterization interval

The estimation data are restricted to the interval between:

1. the time at which the maximum specific growth rate, \(\mu_{\max}\), is
   reached; and
2. the time at which the maximum particle count, \(P_{\max}\), is reached.

This interval defines the post-\(\mu_{\max}\) host-aware characterization
window.

Because the interval length depends on the experiment and construct, the
trajectory vectors need not all have the same number of samples.

---

## Global aggregation

Global fields include:

```text
Mu_mumax_pmax_global_mean
Mu_mumax_pmax_global_std
Pi_mumax_pmax_global_mean
Pi_mumax_pmax_global_std
```

These contain the global mean and standard deviation of the specific growth
rate and synthesis rate over the characterization interval.

---

## Instance aggregation

The field:

```text
Instances
```

is a cell array with one entry for each independent experiment.

Each instance typically includes:

```text
Mu_mumax_pmax_instance_mean
Mu_mumax_pmax_instance_std
Pi_mumax_pmax_instance_mean
Pi_mumax_pmax_instance_std
Wells
```

Instance-level values are means and standard deviations across the retained
wells from one experiment.

---

## Well-level data

Each instance contains a `Wells` cell array. A well entry includes:

```text
Mu_mumax_pmax
Pi_mumax_pmax
```

These are the individual well trajectories used by the `Wells` estimation
mode.

The number of retained wells may vary because of missing data, experimental
quality control, or outlier filtering. Code should not assume a fixed replicate
count.

---

# Complete experimental tensors

Files without the `_micro` suffix contain the complete processed experimental
data.

A typical complete tensor remains a three-dimensional cell array indexed by:

```text
plasmid x promoter x RBS
```

Each nonempty cell contains TU metadata and a nested `Data` structure.

Typical top-level fields are:

```text
TU_Ori
TU_Promoter
TU_RBS
TU_Bioparts
TU_Name
TU_color_code
Data
```

The nested `Data` structure can include:

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

---

## `Stats_Particles`

Contains filtered particle-count trajectories and metadata.

Typical fields include:

```text
List_instances
Global_data_mean
Global_data_std
Global_value_Particlesmax_mean
Global_index_Time_Particlesmax_mean
Mumax_pmax_global_mean
Mumax_pmax_global_std
Data_Mumax_Pi_p_All
Data_Mumax_Phi_s_All
Data_Mumax_Phi_h_All
```

Each instance may include:

```text
Experiment_IDnum
Experiment_date
Colony_number
Status_exp
Data_time
Data_raw
Data_mean
Data_std
Data_Mumax_Pi_p
Data_stats_Mumax_Pi_p
Data_Mumax_Phi_s
Data_stats_Mumax_Phi_s
Data_Mumax_Phi_h
Data_stats_Mumax_Phi_h
Value_Particlesmax_wells
Index_Time_Particlesmax_wells
Data_mumax_pmax_mean
Data_mumax_pmax_std
```

---

## `Stats_MEFL`

Contains fluorescence trajectories in MEFL units:

```text
List_instances
Global_data_mean
Global_data_std
Mumax_pmax_global_mean
Mumax_pmax_global_std
```

Instance entries typically retain:

```text
Data_raw
Data_mean
Data_std
Data_mumax_pmax_mean
Data_mumax_pmax_std
```

---

## `Stats_MEFL_Particles`

Contains cell-normalized fluorescence:

```text
MEFL / Particles
```

Typical fields are:

```text
List_instances
Global_data_mean
Global_data_std
```

---

## `Stats_Mu`

Contains the estimated specific growth rate.

Typical fields include:

```text
List_instances
Global_data_mean
Global_data_std
Global_value_Mumax_mean
Global_index_Time_Mumax_mean
Mumax_pmax_global_mean
Mumax_pmax_global_std
```

Each instance can include:

```text
Data_raw
Data_mean
Data_std
Values_Mumax_wells
Indices_Times_Mumax_wells
Value_Mumax_mean
Index_Time_Mumax_mean
Data_mumax_pmax_mean
Data_mumax_pmax_std
```

---

## `Stats_Pi_p`

Contains the estimated TU synthesis rate.

Typical fields include:

```text
List_instances
Global_data_mean
Global_data_std
Mumax_pmax_global_mean
Mumax_pmax_global_std
```

The synthesis rate in the constitutive-expression libraries is normally
background-corrected using the pARKA21 background model.

---

## `Stats_Phi_s` and `Stats_Phi_h`

These structures contain estimated burden measures:

```text
Phi_s
    strain-level burden

Phi_h
    host-normalized burden
```

They retain experiment-level raw trajectories, means, standard deviations, and
global summaries where available.

---

# Relationship between files

## Lib30 and Lib24

Lib30 is the complete:

```text
2 x 3 x 5
```

combinatorial tensor.

Lib24 is the subset obtained by excluding B0034:

```matlab
indices_rbss_lib24 = [1,2,4,5];
```

Some estimation workflows use `ExpData_Tensor_lib30_micro.mat` directly and
select this subset internally. The dedicated Lib24 files remain useful for
library-specific result generation, plotting, and complete-data analysis.

## Lib6

Lib6 is the B0034 slice of the Lib30 tensor:

```matlab
indices_plasmids_lib6 = [1,2];
indices_promoters_lib6 = [1,2,3];
indices_rbss_lib6 = 3;
```

No separate `ExpData_Tensor_lib6_micro.mat` is required.

## Lib5

Lib5 uses its own dedicated micro tensor because it introduces J23100 in a fixed
pGreen context across five RBS variants.

## Lib35

Lib35 combines the 30 main combinatorial constructs with the five Lib5
constructs for experimental summary and visualization. It is a complete
processed tensor, not an estimation-only micro tensor.

---

# Loading data

Initialize SynTwin before loading data:

```matlab
ROOT = init_SynTwin('experimental',true);
```

Use `SynTwin_path` rather than hard-coded absolute paths:

```matlab
S = load(SynTwin_path( ...
    'Experimental_data', ...
    'ExpData_Tensor_lib30_micro.mat'));
```

Use the exact canonical directory capitalization adopted by the repository.

To inspect the variables stored in a file:

```matlab
whos('-file',SynTwin_path( ...
    'Experimental_data', ...
    'ExpData_Tensor_lib30.mat'))
```

To load only one expected variable:

```matlab
S = load(SynTwin_path( ...
    'Experimental_data', ...
    'ExpData_pARKA21.mat'), ...
    'Data_pARKA21');
```

Validate the field before use:

```matlab
if ~isfield(S,'Data_pARKA21')
    error('Expected Data_pARKA21 in ExpData_pARKA21.mat.');
end
```

---

# Example: inspect a micro tensor

```matlab
S = load(SynTwin_path( ...
    'Experimental_data', ...
    'ExpData_Tensor_lib30_micro.mat'));

R = S.ExpData_Tensor_lib30_micro;

size(R)

TU = R{1,1,1};

disp(TU.TU_Bioparts)
disp(TU.Mu_mumax_pmax_global_mean)
disp(TU.Pi_mumax_pmax_global_mean)
```

Inspect one experimental instance:

```matlab
I = TU.Instances{1};

mu_instance = I.Mu_mumax_pmax_instance_mean;
pi_instance = I.Pi_mumax_pmax_instance_mean;
```

Inspect one well:

```matlab
W = I.Wells{1};

mu_well = W.Mu_mumax_pmax;
pi_well = W.Pi_mumax_pmax;
```

---

# Example: inspect a complete tensor

```matlab
S = load(SynTwin_path( ...
    'Experimental_data', ...
    'ExpData_Tensor_lib35.mat'));

R = S.Data_Exp_Tensor_lib35;
TU = R{1,1,1};
D = TU.Data;

particles = D.Stats_Particles;
growth = D.Stats_Mu;
synthesis = D.Stats_Pi_p;
```

---

# Data integrity checks

Before running an estimation workflow, verify:

1. the expected `.mat` file exists;
2. the expected top-level variable is present;
3. the tensor dimensions match the intended library;
4. nonempty cells contain the required TU metadata;
5. \(\mu\) and \(\Pi\) vectors have compatible lengths;
6. each instance has a valid `Wells` field when `Use_mean='Wells'`;
7. trajectory values used in logarithmic objectives are finite and positive;
8. the correct background correction was applied;
9. no raw Excel dependency remains in a distributed estimation script.

A minimal check is:

```matlab
assert(iscell(R));
assert(~isempty(R{1,1,1}));
assert(isfield(R{1,1,1},'TU_Ori'));
assert(isfield(R{1,1,1},'Instances'));
```

For micro tensors:

```matlab
assert(isfield(R{1,1,1},'Mu_mumax_pmax_global_mean'));
assert(isfield(R{1,1,1},'Pi_mumax_pmax_global_mean'));
```

---

# Naming conventions

```text
ExpData_*
    Processed experimental data.

Tensor
    Multidimensional cell array indexed by biopart context.

lib5, lib24, lib30, lib35
    Library or combined construct collection.

_micro
    Minimum data required by parameter-estimation workflows.

No _micro suffix
    Complete processed experimental structure.
```

The library number denotes the number of transcriptional units represented by
the corresponding library or collection, not necessarily the total number of
nonempty positions in the containing rectangular tensor.

---

# Units and symbols

```text
mu
    Specific growth rate, normally min^-1.

Pi or Pi_p
    TU synthesis rate.

Particles
    Cytation particle-count measurement.

MEFL
    Molecules of equivalent fluorescein.

MEFL/Particles
    Cell-normalized fluorescence.

Phi_s
    Estimated strain burden.

Phi_h
    Estimated host-normalized burden.
```

Refer to the manuscript and Supplementary Information for the precise
biophysical definitions and normalization conventions.

---

# MATLAB requirements

Reading the `.mat` files requires MATLAB.

The downstream preprocessing and visualization workflows may additionally use:

- local functions in scripts;
- `tiledlayout`;
- `exportgraphics`;
- Statistics and Machine Learning Toolbox functions;
- SynTwin functions from `Scripts_base/`.

GNU Octave compatibility has not been tested.

---

# Further documentation

The original SynTwin User Manual documents:

- the micro experimental tensor in Appendix A;
- the complete experimental tensor in Appendix B;
- the result tensor in Appendix C.

The present README reflects the current contents of `Experimental_data/`,
including the later pARKA21 background dataset and the current distribution
policy. Where the old manual and current scripts differ, the current scripts
and per-workflow READMEs define the distributed interface.
