# Jacobian analysis

## Overview

This directory contains MATLAB workflows for local practical-identifiability analysis and for visualising the combinatorial parameter--construct structure of the SynTwin libraries.

The numerical analyses build experimental Jacobians from stored synthesis-rate sensitivities. The structural plot is complementary: it displays which global parameters are associated with each expression construct, but it does not use numerical sensitivity values.

## Main scripts

### `Analyze_Libs_Jacobian_ranks.m`

Compares the experimental Jacobian ranks of Lib30 and two subsets derived from it: Lib24 and Lib6.

The script loads the Lib30 sensitivity tensor and builds one Jacobian for each selected library using the same experimental aggregation layer. It reports the number of free parameters, the numerical rank, and an SVD-based effective rank with

```text
tol = 1e-2 * s(1)
```

Default input:

```text
Estimation_Pi/L30_L1O_reduced_model/
    Results_Tensor_Lib30_L1O_reduced_Wells.mat
```

Expected variable:

```text
Results_Tensor_Lib30_L1O_reduced
```

Run:

```matlab
Analyze_Libs_Jacobian_ranks
```

### `Analyze_Lib24_reduced_Jacobian.m`

Performs a detailed practical-identifiability analysis for Lib24 under the reduced translation model.

It builds both

```text
absolute Jacobian: dPi/dtheta
relative Jacobian: (theta/Pi) dPi/dtheta
```

and reports exact rank, effective rank, condition number, least identifiable singular-vector direction, approximate parameter standard deviations, parameter correlations, and the least identifiable parameter combination.

Default input:

```text
Estimation_Pi/L24_ALL_reduced_model/
    Results_Tensor_Lib24_ALL_reduced_model_Wells.mat
```

Expected variable:

```text
Results_Tensor_Lib24_ALL_reduced
```

Run:

```matlab
Analyze_Lib24_reduced_Jacobian
```

## Experimental aggregation layer

Both analysis scripts support:

```matlab
DATA_LAYER = 'Global';
DATA_LAYER = 'Instances';
DATA_LAYER = 'Wells';
```

`Global` uses one mean trajectory per construct. `Instances` concatenates experiment-level means. `Wells` concatenates all well-level trajectories.

## Combinatorial structure plot

### `Plot_Combinatorial_Library_Jacobian_Structure.m`

Generates a binary parameter--construct incidence matrix for the full 35-construct combination:

```text
L24 + L6 + L5
```

The plot contains 35 expression constructs on the horizontal axis and the global promoter, RBS, and plasmid parameter families on the vertical axis. A filled cell means that the construct depends structurally on that parameter.

The parameter families are promoter transcription rates, RBS intrinsic initiation capacities, RBS inverse scaling parameters, and effective plasmid copy numbers.

Run:

```matlab
Plot_Combinatorial_Library_Jacobian_Structure
```

Generated files:

```text
Combinatorial_Library_Jacobian_Structure.png
Combinatorial_Library_Jacobian_Structure.svg
```

## Additional plotting script

`plot_jacobian_comparison_blocks.m` is retained as an auxiliary visualisation comparing the full incidence structure with the J23100 sublibrary. The transposed complete-library plot is the preferred representation of the combinatorial structure.

## Helper functions

- `build_jacobian.m` builds the sparse Jacobian for the full parameterisation.
- `build_jacobian_reduced.m` builds absolute and relative Jacobians for the reduced parameterisation.
- `get_free_params_for_library.m` returns the free promoter, RBS, and plasmid parameters for L24, L30, L6, L4, and L5.

## Interpretation

A full-column-rank Jacobian indicates local distinguishability of the selected free parameters at the analysed point and data layer. Rank deficiency indicates that at least one local parameter combination cannot be distinguished from the available synthesis-rate sensitivities.

The effective rank is scale dependent and should be interpreted together with the exact rank, condition number, singular vectors, sensitivity scaling, and experimental aggregation layer. This analysis is a local practical-identifiability diagnostic; it does not establish global identifiability.

## Requirements

- MATLAB with local functions in scripts.
- SynTwin initialisation through `init_SynTwin('experimental',true)`.
- Access to the required result tensors through `SynTwin_path`.
- `svd`, `rank`, and standard plotting functions.
- `exportgraphics` for preferred PNG and SVG export.

GNU Octave compatibility has not been tested.
