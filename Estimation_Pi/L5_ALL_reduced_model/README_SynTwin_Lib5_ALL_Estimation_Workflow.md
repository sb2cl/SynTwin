# SynTwin — Lib5 (ALL) estimation workflow (reduced model)

This folder contains the **SynTwin workflow** to estimate the **promoter transcription parameter** for the **Lib5** sublibrary using the **ALL** scheme (all TUs are used simultaneously), under the **reduced model**.

Lib5 is a small calibration sublibrary used to infer the promoter strength **Ω (Omega)** for **J23100** in a fixed genetic background, while **all other biopart parameters are inherited** from previous SynTwin estimations.

---

## Context: what is Lib5?

**Lib5** consists of **5 transcriptional units (TUs)** defined by:

- **1 plasmid origin:** pGreen (single Ori context)
- **1 promoter:** J23100 (the one to estimate)
- **5 RBSs:** a panel of five RBS parts (including **B0034**)

So each TU is: **pGreen + J23100 + RBSᵢ** (i = 1…5).

### What is estimated vs inherited?

**Estimated (Lib5, ALL):**
- **Ω (Omega)** — promoter transcription strength for **J23100**  
  (`params(1) = Omega`)

**Inherited (fixed during each cost evaluation):**
- **Gene copy number (Gene_cn)** for pGreen contexts  
  inherited from **Lib24 L1O reduced-model results**
- **RBS intrinsic initiation capacity (k0_sigma0)** for the 5 RBSs:
  - For RBSs present in Lib24: inherited from **Lib24 L1O reduced-model results**
  - For **B0034** (not available in Lib24 in this setup): inherited from **Lib6** estimation results

**Fixed (inside scripts / cost function):**
- `inv_sigma0` (RBS sensitivity-related parameter), typically set to `0.02`
- Other digital-twin constants in `model_c` (e.g., GFP parameters, etc.)

This design lets you **estimate Ω(J23100)** while **propagating uncertainty** from the inherited parameters.

---

## Files in this folder

Typical workflow (run in this order):

1. **`Estimate_Lib5_ALL_reduced_model.m`**  
   Runs repeated BADS optimizations to estimate **Ω(J23100)** using all 5 TUs.

2. **`Generate_Results_Lib5_ALL_reduced_model.m`**  
   Post-processes estimation outputs and compiles a **Results tensor** that includes:
   - inherited parameter statistics,
   - estimated Ω(J23100) statistics,
   - experimental Mu(t) and Pi(t),
   - SynTwin predictions, sensitivities, and Monte Carlo summaries.

3. **`Show_Results_Lib5_ALL_reduced_model.m`** (if present)  
   Loads the Results tensor and generates plots / summaries.

4. **`J4_LogPI_Lib5_ALL_reduced.m`**  
   Objective (cost) function called by BADS/MEIGO.

---

## Before you run: initialize SynTwin and paths

From **SynTwin root**, start MATLAB and run:

```matlab
ROOT = init_SynTwin('experimental',true,'bads',true);
```

If you run scripts **from SynTwin root**, you must also add the script folder to the MATLAB path:

```matlab
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);
```

(Each workflow script in this folder already contains this pattern.)

---

## Step 1 — Run the estimation (ALL)

Open and run:

```matlab
Estimate_Lib5_ALL_reduced_model
```

### Key options to set inside the script

- **Data aggregation level** (`Use_mean`):
  - `'Global'`    — one mean trajectory per TU
  - `'Instances'` — one mean trajectory per experiment (instance)
  - `'Wells'`     — well-level trajectories (most detailed)

- **Number of runs** (`num_runs`):  
  Typical values: `50–200` depending on desired Monte Carlo resolution.

### What gets saved

The estimation script creates (if needed) a local folder:

```
./Estimated_results/
```

and saves:

```
Results_BADS_Lib5_ALL_reduced_<Use_mean>.mat
```

This file contains a cell array of per-run solutions:

```matlab
Results_BADS_J23100_ALL_approx{num_run}.results = [Omega, Jmin]
```

---

## Step 2 — Generate the compiled Results tensor

After estimation, run:

```matlab
Generate_Results_Lib5_ALL_reduced_model
```

This compiles:

- the estimated **Ω(J23100)** distribution (mean/std + MC samples),
- inherited parameter distributions (Gene_cn and RBS k0_sigma0),
- experimental data for the 5 TUs,
- SynTwin predictions + sensitivities,
- Monte Carlo prediction summaries over a μ grid.

### Output

A Results tensor is saved in the same folder as the script:

```
Results_Tensor_Lib5_ALL_reduced_model_<Use_mean>.mat
```

This tensor is typically loaded by `Show_Results_*` scripts and is also the recommended starting point for downstream analysis.

---

## Parameter inheritance and uncertainty propagation

This workflow **inherits parameters from two independent previous estimations**:

1. **Lib24 (L1O, reduced model)**  
   Provides pGreen-related `Gene_cn` and the `k0_sigma0` values for RBSs that exist in Lib24.

2. **Lib6 (reduced model)**  
   Provides `k0_sigma0` for **B0034**.

To propagate uncertainty, the estimation script typically:
- selects random Monte Carlo indices from the inherited parameter samples,
- uses those samples consistently within each run `num_run`,
- estimates Ω(J23100) conditioned on that inherited parameter realization.

---

## Dependencies

Required within the SynTwin distribution:

- **HEM surrogate**  
  `Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat`

- **Experimental data (Lib5)**  
  `Experimental_Data/ExpData_Tensor_lib5_micro.mat`

- **Inherited parameters**  
  - `L24_L1O_reduced_model/Results_Tensor_Lib24_L1O_reduced_*.mat` (pGreen + RBSs)
  - `L6_* / Results_Tensor_Lib6_*_*.mat` (B0034)

Optimizers:

- **BADS** (bundled in SynTwin)
- MATLAB **Parallel Computing Toolbox** is optional (if absent, replace `parfor` with `for`)

---

## Notes and common pitfalls

- **Be consistent with `Use_mean`** across Estimate / Generate / Show scripts.
- If you see errors like “file not found”, verify:
  - SynTwin was initialized via `init_SynTwin(...)`,
  - the script folder is on the MATLAB path,
  - `./Estimated_results/` exists and contains the expected `Results_BADS_*` file.
- The cost function is called many times; it should not generate figures or interactive output.

---

## Where to find data structure documentation

This README focuses on the workflow logic and the file outputs.  
For full details on the **Results tensor fields** and data structures, see the **SynTwin “Software and Data Supplementary”** document shipped with this release.
