# Show_Results

## Overview

This directory contains MATLAB scripts that regenerate the main result visualisations associated with Figures 2b--c, 3, and 4 of:

> **Host-aware Identification of Intrinsic Gene Expression Biopart Parameters using Combinatorial Libraries**

The scripts are named by their scientific function rather than by figure number. The exported files retain the manuscript panel identifiers so that each output can be mapped directly to the paper.

The scripts load processed literature data or stored parameter-estimation results. They do **not** rerun the global optimisation or parameter-identification workflows.

## Contents

### `Plot_Host_Physiology_Model_Validation.m`

Validates two fixed host-physiology relationships against literature data:

- peptide elongation rate, $\nu_t(\mu)$;
- mature active ribosomal protein-mass fraction, $\Phi_R(\mu)$.

**Paper association:** Figures 2b and 2c.

The evaluated models are:

```text
nu_t(mu) = gamma1*mu/(gamma2 + mu)
Phi_R(mu) = a + b*mu
```

with the fixed values reported in Supplementary Section 3.5:

```text
a           = 0.0660
b           = 7.1502 min
gamma1      = 25.9083 aa/s
gamma2      = 0.0092 min^-1
Phi_t*Phi_m = 0.72
```

Inputs are resolved from:

```text
Generate_HEM/Third_party_data/
```

Required spreadsheet stems:

```text
ecoli_peptide_elongation_rates_updated
ecoli_ribosomal_mass_fractions_updated
```

Both `.xlsx` and `.xls` extensions are accepted.

Generated files:

```text
Figures/Figure_2b_peptide_elongation.png
Figures/Figure_2b_peptide_elongation.svg
Figures/Figure_2c_ribosomal_mass_fraction.png
Figures/Figure_2c_ribosomal_mass_fraction.svg
```

### `Plot_Biopart_Parameter_Distributions.m`

Plots the stored distributions of:

- effective plasmid copy number, $N$;
- promoter transcription rate, $\omega$;
- RBS intrinsic initiation capacity, $\kappa^0=K^0/\sigma^0$.

**Paper association:** Figures 3a, 3b, and 3c.

Required result tensors:

```text
Estimation_Pi/L24_L1O_reduced_model/
    Results_Tensor_Lib24_L1O_reduced_Wells.mat

Estimation_Pi/L6_L1O_reduced_model/
    Results_Tensor_Lib6_L1O_reduced_Wells.mat

Estimation_Pi/L5_ALL_reduced_model/
    Results_Tensor_Lib5_ALL_reduced_Wells.mat
```

The script extracts:

- `Gene_cn_raw` from L24 for the pGreen-derived plasmid backbone;
- `Omega_raw` from L24 for J23106, J23102, and J23101;
- `Omega_raw` from L5 for J23100;
- `RBS_k0_sigma0_raw` from L24 for B0030, B0032, J61100, and J61101;
- `RBS_k0_sigma0_raw` from L6 for B0034.

Each violin contains the stored raw samples. The black overlays show the median, interquartile range, and 5th--95th percentile interval.

Generated files:

```text
Figures/Figure_3a_plasmid_copy_number.png
Figures/Figure_3a_plasmid_copy_number.svg
Figures/Figure_3a_plasmid_copy_number.pdf

Figures/Figure_3b_promoter_strengths.png
Figures/Figure_3b_promoter_strengths.svg
Figures/Figure_3b_promoter_strengths.pdf

Figures/Figure_3c_RBS_IIC.png
Figures/Figure_3c_RBS_IIC.svg
Figures/Figure_3c_RBS_IIC.pdf
```

### `Plot_RBS_Host_Aware_Translation_Maps.m`

Uses the stored L24 and L6 digital-twin predictions to plot:

- the host-aware mapping from IIC, $\kappa^0$, to ERAC, $K_A^t/\nu_t$, across growth rates;
- the predicted operational ETR range, $K_A^t$, for each characterised RBS.

**Paper association:** Figures 4a and 4b.

Required result tensors:

```text
Estimation_Pi/L24_L1O_reduced_model/
    Results_Tensor_Lib24_L1O_reduced_Wells.mat

Estimation_Pi/L6_L1O_reduced_model/
    Results_Tensor_Lib6_L1O_reduced_Wells.mat
```

The script searches the complete result tensors for entries containing `Synthesis_predictions`. Representative fields used by this workflow include:

```text
TU_RBS
Gene_cn_mean
Omega_mean
Mu_mumax_pmax_global_mean
Pi_mumax_pmax_global_mean
Synthesis_predictions.Mu_values
Synthesis_predictions.fs_pred_values
Synthesis_predictions.KA_t_pred_values
Synthesis_predictions.varphi_pred_values
Synthesis_predictions.load_pred_values
```

Generated files:

```text
Figures/Figure_4a_ERAC_projection.png
Figures/Figure_4a_ERAC_projection.svg
Figures/Figure_4a_ERAC_projection.pdf

Figures/Figure_4b_ETR_range_by_IIC.png
Figures/Figure_4b_ETR_range_by_IIC.svg
Figures/Figure_4b_ETR_range_by_IIC.pdf
```

## Running the scripts

Initialise SynTwin as required by the repository and run any script from MATLAB:

```matlab
Plot_Host_Physiology_Model_Validation
Plot_Biopart_Parameter_Distributions
Plot_RBS_Host_Aware_Translation_Maps
```

Each script calls `init_SynTwin`, resolves repository paths through `SynTwin_path` or the returned repository root, and writes its output to the local `Figures/` directory.

## Result-tensor scope

The MAT files are complete structured outputs from the host-aware parameter-estimation workflows, not reduced plotting tables. Depending on the library, populated tensor entries can include:

```text
TU_Ori
TU_Promoter
TU_RBS
TU_Bioparts
TU_Name
Parameters_local_raw
Parameters_local_mean
Parameters_local_std
Gene_cn_raw
Gene_cn_mean
Gene_cn_std
Omega_raw
Omega_mean
Omega_std
RBS_k0_sigma0_raw
RBS_k0_sigma0_mean
RBS_k0_sigma0_std
S_Pi_NA_global_values
S_Pi_Omega_global_values
S_Pi_RBS_k0_sigma0_global_values
Synthesis_predictions
Gene_cn_MC_samples
Omega_MC_samples
RBS_k0_sigma0_MC_samples
MC_samples
MC_mu_slices
```

The complete schema and parameter-estimation workflows are documented in the SynTwin User Manual.

## Requirements

- MATLAB R2020a or later is recommended.
- `Plot_Host_Physiology_Model_Validation.m` requires spreadsheet import and plotting functions available in core MATLAB.
- `Plot_Biopart_Parameter_Distributions.m` requires `violinplot`, `exportgraphics`, and `prctile`; `prctile` is provided by Statistics and Machine Learning Toolbox.
- `Plot_RBS_Host_Aware_Translation_Maps.m` requires `scatteredInterpolant`, `exportgraphics`, and support for local functions in scripts.
- No optimisation toolbox is required because the scripts load stored results.
- GNU Octave compatibility has not been tested.

## Software and archived release

- GitHub: `https://github.com/sb2cl/SynTwin`
- Zenodo: `https://doi.org/10.5281/zenodo.18787107`
