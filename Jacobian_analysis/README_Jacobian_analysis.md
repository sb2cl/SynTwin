# Jacobian analysis (practical identifiability)

This folder contains a MATLAB script to assess *practical identifiability* of SynTwin library parameters by computing the rank of the experimental Jacobian (sensitivity) matrix.

## What it does

`Analyze_Libs_Jacobian_ranks.m`:

- Loads **Lib30** results that already include *absolute sensitivities* of the synthesis rate `Pi` w.r.t. model parameters.
- Builds the experimental Jacobian `J` for **Lib30**, and for the **Lib24** and **Lib6** subsets.
- Reports:
  - `rank(J)` (exact numerical rank)
  - an SVD-based **effective rank** using `tol = 1e-2 * s(1)`.

## Requirements

- SynTwin initialized from MATLAB:

```matlab
ROOT = init_SynTwin('experimental',true);
```

- A Lib30 results `.mat` file containing the sensitivity tensor. By default, the script expects:

- File: `L30_L1O_reduced_model/Results_Tensor_Lib30_L1O_reduced_Wells.mat`
- Variable: `Results_Tensor_Lib30_L1O_reduced`

If your distribution uses a different file name or variable name, edit `RESULTS_FILE` and `RESULTS_VAR_NAME` inside the script.

## Usage

From the `Jacobian_analysis` folder:

```matlab
Analyze_Libs_Jacobian_ranks
```

To change the experimental aggregation layer used to build Jacobian rows, edit:

```matlab
DATA_LAYER = 'Wells';   % 'Global' | 'Instances' | 'Wells'
```

## Notes

- This analysis assumes the Jacobian is constructed from *measured-growth* trajectories `Mu(t)` and corresponding `Pi` sensitivities already computed in the Lib30 workflow.
- Rank deficiency indicates that some parameter combinations are locally unidentifiable given the selected data layer.
