# SynTwin – L6 Estimation Workflow (L1O, Reduced Model)

This folder implements the SynTwin workflow for **Leave-One-Out (L1O) cross-validation**
on the **L6 sublibrary** under the **reduced model**.

For detailed descriptions of the experimental-data and results tensor structures, refer to the
*Software and Data Supplementary* document shipped with SynTwin.

---

## L6 context (6-TU sublibrary)

L6 is defined by:

- **2 plasmid origins**
- **3 promoters**
- **1 shared RBS** (single RBS across all constructs; in Lib30 indexing this is `RBS=3`)

Total transcriptional units (TUs): **2 × 3 × 1 = 6**

### Key methodological feature (parameter inheritance)

In L6 we **only estimate the shared RBS intrinsic initiation capacity** `k0_sigma0`:

- **Estimated:** `k0_sigma0` (shared RBS)
- **Fixed:** `inv_sigma0` (sensitivity-related parameter; fixed inside the cost function)

All **Ori-/promoter-dependent parameters** (e.g., `Omega`, `Gene_cn`) are **inherited from Lib24**
inference results (loaded from `Results_Tensor_Lib24_L1O_reduced_*.mat`).

Experimental measurements for the 6 L6 constructs are retrieved by indexing the **Lib30-format**
experimental tensor container (`ExpData_Tensor_lib30_micro`).

---

## Workflow summary (L1O)

L1O runs **6 folds** (one per construct). In each fold:

1. One TU is **left out** (validation).
2. The RBS parameter `k0_sigma0` is estimated using the remaining **5 training TUs**.
3. This fold-specific estimate is stored as a **local** (fold-specific) result for the left-out TU.

---

## Scripts

### 1) Parameter estimation (per-fold)

`Estimate_Lib6_L1O_reduced_model.m`

- Iterates over the 6 constructs in L6 (each becomes `Construct_2_LO = [Ori, Prom, RBS]`).
- For each fold, runs `num_runs` independent BADS optimizations.
- Each run propagates uncertainty from inherited Lib24 parameters by selecting a consistent
  Monte Carlo sample index `num_run` for:
  - `Omega_MC_samples(num_run)`
  - `Gene_cn_MC_samples(num_run)`

**Outputs (saved to disk):**

One file **per left-out construct**:

```
./Estimated_results/Results_BADS_Lib6_L1O_reduced_<Ori><Prom><RBS>_<Use_mean>.mat
```

Each file stores a cell array of per-run results, e.g.:

```
Results_BADS_...{num_run}.results = [k0_sigma0, Jmin]
```

---

### 2) Results compilation (tensor build + predictions)

`Generate_Results_Lib6_L1O_reduced_model.m`

- Loads the per-fold estimation files from `./Estimated_results/`.
- Builds the tensor `Results_Tensor_Lib6_L1O_reduced` (size 2×3×1).
- Stores **fold-specific ("local")** statistics for each TU:
  - `Parameters_local_raw / mean / std`
  (i.e., inferred when that TU was left out).
- Also pools all folds to compute **global ("ALL")** RBS statistics:
  - `RBS_k0_sigma0_raw / mean / std`
- Copies inherited Lib24 statistics (`Omega_*`, `Gene_cn_*`) into each TU entry.
- Attaches experimental data (Global / Instances / Wells) for each TU.
- Computes:
  - digital-twin synthesis-rate predictions on experimental `Mu(t)`,
  - sensitivities (relative → absolute),
  - Monte Carlo propagated predictions on a `Mu` grid (KDE-based PDFs and quantiles).

**Output (saved to disk):**

```
Results_Tensor_Lib6_L1O_reduced_model_<Use_mean>.mat
```

---

### 3) Visualization

`Show_Results_Lib6_L1O_reduced_model.m`

- Loads the generated L6 L1O tensor.
- Produces diagnostic plots and summary visualizations.

---

## How to run

From any working directory:

```matlab
ROOT = init_SynTwin('experimental',true,'bads',true);

% If running from SynTwin root, ensure this folder is on the path:
addpath(fileparts(mfilename('fullpath')));

Estimate_Lib6_L1O_reduced_model
Generate_Results_Lib6_L1O_reduced_model
Show_Results_Lib6_L1O_reduced_model
```

---

## Configuration

Inside the scripts, set:

- `Use_mean` (aggregation level)
  - `'Global'`
  - `'Instances'`
  - `'Wells'`

- `num_runs` (optimizations per fold; also used as the Monte Carlo index for inherited parameters)

---

## Dependencies

- HEM surrogate model  
  `Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat`

- Experimental data container (Lib30 tensor format)  
  `Experimental_Data/ExpData_Tensor_lib30_micro.mat`

- Inherited Lib24 results (parameter distributions/statistics)  
  `L24_L1O_reduced_model/Results_Tensor_Lib24_L1O_reduced_*.mat`

- BADS optimizer (bundled in SynTwin)

---

## Notes for interpretation

- **Local (fold-specific) RBS estimates** are stored per TU via `Parameters_local_*`.
  These correspond to *training on the other 5 constructs* and can be used to assess
  predictive generalization for the left-out construct.

- **Global ("ALL") RBS statistics** pool samples from all folds and are useful as a
  compact summary, but they are not a single ALL-fit (they aggregate L1O folds).

- For full documentation of the fields stored inside `Results_Tensor_Lib6_L1O_reduced`,
  see the *Software and Data Supplementary*.
