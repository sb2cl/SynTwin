# SynTwin — Host-aware Digital Twin for Biopart Characterization

> **Release note (pre-publication / peer-review stage).**  
> This repository corresponds to the SynTwin implementation used to generate the computational results reported in the associated manuscript and Supplementary Information currently under review at *Nature Communications*.  
> Some additional capabilities under active development may be intentionally excluded from this release.

---

## What SynTwin provides (this release)

SynTwin is a MATLAB-based framework for **host-aware characterization** of genetic bioparts (plasmid origins, promoters, RBSs) using:

- processed experimental observables (growth rate μ(t) and synthesis rate Π(t)),
- a host-aware **E. coli digital twin**,
- global optimization (BADS or MEIGO/eSS),
- reproducible **parameter-estimation workflows** and **Results Tensor** outputs.

This release focuses on **parameter estimation and result reproduction** for the libraries analyzed in the manuscript (Lib24, Lib30, Lib6, Lib5). The codebase also contains internal components used to build the **Host Equivalent Model (HEM)** required by the digital twin.

---

## Repository structure

Top-level directory layout (main folders):

- `Scripts_base/`  
  Core helper functions used across workflows (e.g., digital-twin prediction utilities).

- `Experimental_Data/`  
  Processed experimental data containers (`ExpData_Tensor_*.mat`) used for estimation.

- `Generate_HEM/`  
  Host-model preparation utilities and third-party datasets (see attribution below).
  - `HEM_Surrogate/`  
    Core HEM surrogate functions and `HEM_Surrogate.mat`.
  - `Third_party_data/`  
    Verbatim copies / re-encodings of upstream datasets used to fit/validate HEM.
  - `Estimate_WT/`, `Nut_estimate/`, `Mass_estimate/`  
    Scripts used to fit/validate the wild-type host baseline and associated host relationships.

- `Estimation_Pi/`  
  Ready-to-run estimation workflows, organized by library / validation scheme / model variant, e.g.:
  - `Lib30_L1O_full_model/`, `Lib30_ALL_full_model/`, `Lib30_L1O_reduced_model/`, `Lib30_ALL_reduced_model/`
  - `Lib24_ALL_full_model/`, `Lib24_L1O_reduced_model/`, `Lib24_ALL_reduced_model/`
  - `Lib6_L1O_reduced_model/`, `Lib6_ALL_reduced_model/`
  - `Lib5_L1O_reduced_model/`, `Lib5_ALL_reduced_model/`

- `Jacobian_analysis/`  
  Utilities to assess identifiability / Jacobian rank structure.

- `Matlab/`  
  Third-party MATLAB toolboxes (see licensing section):
  - `MEIGO/`
  - `bads-master/`

- `init_SynTwin.m`, `SynTwin_root.m`, `SynTwin_path.m`  
  Portable initialization helpers (no absolute paths).

---

## Quick start (MATLAB)

From any working directory:

```matlab
% Initialize SynTwin paths (portable; no absolute paths)
ROOT = init_SynTwin();

% OPTIONAL: add the folder containing the specific workflow script you want to run
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);
```

Then run one of the workflows in `Estimation_Pi/...`:

```matlab
Estimate_LibXX_<Validation>_<Model>.m
Generate_Results_LibXX_<Validation>_<Model>.m
Show_Results_LibXX_<Validation>_<Model>.m
```

See also: `Estimation_Pi/README_SynTwin_Estimation_Workflow_L24_and_L30.md`.

---

## Reproducing manuscript results

This repository contains:

- the processed experimental tensors used for inference,
- the workflow scripts that generate the estimation outputs and Results Tensors,
- the host-model components (WT + HEM surrogate) used by the digital twin.

If you reproduce or extend the results, please cite the SynTwin manuscript (see *Citation* below) and the third-party tools/datasets as appropriate.

---

## Third-party datasets and attribution

SynTwin includes (or may include) verbatim copies / re-encodings of third-party datasets used to fit/validate the HEM component.

**You must follow the original license/terms of the upstream sources when reusing them.**

A repository-level statement is provided in:

- `Data Availability & Attribution Statement.md`

This statement documents:
- Chure & Cremer (Flux Parity) collated datasets (ribosomal mass fractions and peptide elongation rates) and their Zenodo DOI,
- upstream sources such as Bremer & Dennis (EcoSal Plus, 2008),
- where these files are stored in the repository, and how they are used.

---

## Third-party software

SynTwin interfaces with and/or redistributes third-party optimization packages:

- **MEIGO (eSS – Enhanced Scatter Search)**  
- **BADS (Bayesian Adaptive Direct Search)**

These packages retain their original authorship and licenses. Users are responsible for complying with those terms when using or redistributing them. When results rely on these solvers, users are encouraged to cite the corresponding publications in addition to the SynTwin paper.

---

## Citation

If you use SynTwin in academic work, please cite the associated manuscript:

> **[Authors]** (2026). *Host-aware Identification of Intrinsic Gene Expression Biopart Parameters using Combinatorial Libraries.* **Nature Communications** (submitted).  
> DOI: *to be added upon publication.*

A `CITATION.cff` file will be added once the DOI is available.

---

## License and terms of use (SynTwin code)

**SynTwin code (this repository):**  
See `LICENSE` at the repository root.


## Contact

For questions, bug reports, or attribution concerns, please open an issue on the repository or contact the corresponding author listed in the manuscript.
