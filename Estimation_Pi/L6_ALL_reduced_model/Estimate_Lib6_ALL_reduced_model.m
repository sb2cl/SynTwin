%% Estimate_Lib6_ALL_reduced_model
% SynTwin workflow script: parameter estimation for the L6 sublibrary using ALL TUs
% under the reduced model (estimate a single shared-RBS k0_sigma0; fix inv_sigma0).
%
% CONTEXT (L6 sublibrary)
%   L6 is a 6-TU sublibrary defined by 2 plasmid origins and 3 promoters,
%   sharing a common RBS (6 combinations in total). In this workflow, the
%   origin- and promoter-dependent parameters (e.g., Omega and effective gene
%   copy number) are NOT re-estimated; instead, they are inherited from Lib24
%   inference results.
%
% DESCRIPTION
%   Runs <num_runs> independent estimations of the shared-RBS intrinsic initiation
%   capacity (k0_sigma0) using ALL TUs in the L6 sublibrary. Each run uses BADS to
%   minimize the library-scale mismatch between experimental synthesis-rate data
%   and SynTwin digital-twin predictions (Pi).
%
%   Importantly, the workflow propagates uncertainty from inherited Lib24
%   parameters by selecting, for each run index num_run, the corresponding Monte
%   Carlo sample of:
%       - Omega (context-dependent promoter/transcription term), and
%       - Gene_cn (effective gene copy number),
%   and keeping these inherited values fixed while estimating the L6 RBS parameter.
%
% MODEL SETUP (reduced model for RBSs)
%   - Estimated: k0_sigma0 (intrinsic initiation capacity of the shared RBS)
%   - Fixed:     inv_sigma0 (sensitivity-related parameter; fixed inside the cost function)
%
% INPUTS
%   None (configuration is set inside the script).
%   The user must set:
%     - Use_mean: aggregation level ('Global', 'Instances', or 'Wells')
%     - num_runs: number of Monte Carlo / optimization runs (e.g., 100)
%
% OUTPUTS (saved to disk)
%   ./Estimated_results/Results_BADS_Lib6_ALL_reduced_<Use_mean>.mat
%     Contains a cell array of per-run results:
%       Results_BADS_B0034_ALL_approx{num_run}.results = [k0_sigma0, Jmin]
%
% DEPENDENCIES
%   - SynTwin must be initialized before running this script:
%       ROOT = init_SynTwin(...);
%   - This folder must be on the MATLAB path (if running from SynTwin root):
%       addpath(fileparts(mfilename('fullpath')));
%   - Requires the HEM surrogate:
%       Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
%   - Experimental data container (L6 is extracted by indexing):
%       Experimental_Data/ExpData_Tensor_lib30_micro.mat
%   - Inherited Lib24 parameter samples (Omega and Gene_cn):
%       Results_Tensor_Lib24_*.mat (loaded and mapped into ParamsData_lib24)
%   - Uses: J4_LogPI_L6_ALL_reduced (objective function)
%   - Requires BADS (bundled in SynTwin distribution) and optional Parallel TB.
%
% USAGE
%   Estimate_L6_ALL_reduced
%
% NOTES
%   - The L6 experimental measurements are retrieved from the Lib30 tensor-format
%     container by selecting the appropriate indices (2 oris, 3 promoters, RBS=3).
%   - If Parallel Toolbox is unavailable, replace parfor with for.
%   - This script is designed to quantify the impact of uncertainty in inherited
%     Lib24 parameters on the inferred shared-RBS initiation capacity in L6.

clearvars;
close all;
dbstop if error
warning on
% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin('experimental',true,'bads',true);
% Ensure this script folder is on the path (for local functions in this folder)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);

% --- Load data using portable absolute paths ---
load(SynTwin_path('Generate_HEM','HEM_Surrogate','HEM_Surrogate.mat'));          % loads data of the Host Equivalent Model (HEM) 
% Getting the data of Lib30 WE USE THIS AS BASE TO GET THE EXPERIMENTAL DATA EVEN IF WE ONLY TAKE A
% SUB-LIBRARY OF 24 TUs HERE
load(SynTwin_path('Experimental_Data','ExpData_Tensor_lib30_micro.mat'));        % loads data of Lib30: ExpData_Tensor_lib30_micro

% Getting the data of Lib24 L1O reduced model (we ONLY use the estimated parameters from this file): 
load(SynTwin_path('Estimation_Pi/L24_L1O_reduced_model','Results_Tensor_Lib24_L1O_reduced_Wells.mat'));      

% Consider an exogenous circuit. We use a Transcriptional Unit (TU) expressing GFP:
model_c.lp_c = 240; %Length of GFP protein (aa)     
model_c.le_c = model_c.lp_c^0.097/0.0703; 	%Ribosome occupancy length (aa) %lea=la^0.097/0.0703
model_c.dm_c  =  0.2; %Mean degradation rate of non-ribosomal mRNA (1/min)
model_c.Em_c =  model_c.lp_c/model_c.le_c*(1- (model_c.lp_c/(model_c.lp_c+model_c.le_c))^(model_c.lp_c/model_c.le_c)) ;  
model_c.WEm_c =  1 + 1/model_c.Em_c;  
model_c.N_pSC101 = 5; %known

options = bads('defaults');
options.Display='final';

% Options are: 
% - 'Global' (global mean for each construct),
% - 'Instances' (mean of each experiment for each construct),
% - 'Wells' (use data of all individual culture wells)

%Use_mean = 'Instances';  
%Use_mean = 'Global';  
Use_mean = 'Wells'; 
num_runs = 100;  %length(ParamsData_Tensor_lib24_micro{1,1,1}.Gene_cn_MC_samples) is 1000


indices_plasmids = [1,2];
indices_promoters = [1,2,3];

% NOTICE each combination (i,j,k) below will correspond to a library transcriptional unit (TU)
for i=1:length(indices_plasmids) %Plasmids i=1 high copy (pGreen), i=2-> low copy (pSC101)
  for j=1:length(indices_promoters)
    ParamsData_lib24{i,j}.Gene_cn_MC_samples = Results_Tensor_Lib24_L1O_reduced{i,j,1}.Gene_cn_MC_samples;
    ParamsData_lib24{i,j}.Omega_MC_samples = Results_Tensor_Lib24_L1O_reduced{i,j,1}.Omega_MC_samples;
  end
end

% params: vector of estimated parameters of the model
 % params(1) = RBS_k0_sigma0, RBS IIC of (B0034)
 x0 = 0.13; % Starting point for RBS_k0_sigma0
 lb = 0.05;   %lower expected bounds for RBS_k0_sigma0,      B0034         
 ub = 0.5;  %upper expected bounds for RBS_k0_sigma0,  B0034
 plb = 0.08;  %Plausible lower bounds bounds for RBS_k0_sigma0, 
 pub =  0.16; %Plausible upper bounds bounds for RBS_k0_sigma0

 % Run BADS, which returns the minimum X and its value FVAL.

parfor num_run=1:num_runs 
    J=@(parameters) J4_LogPI_L6_ALL_reduced(parameters,model_c,Use_mean,ExpData_Tensor_lib30_micro,ParamsData_lib24, HEM,num_run);
    [params, Jmin_value] = bads(J,x0,lb,ub,plb,pub,[],options);
     Results_BADS_B0034_ALL_approx{num_run}.results=[params, Jmin_value]; 
end

% --- Save results to local Estimated_results folder (portable) ---
results_dir = fullfile(SCRIPT_DIR,'Estimated_results');
if ~exist(results_dir,'dir')
    mkdir(results_dir);
end
file_name = "Results_BADS_Lib6_ALL_reduced_"  + Use_mean + ".mat";
file_tensor = fullfile(results_dir, file_name);
save(file_tensor, "Results_BADS_B0034_ALL_approx", "-mat");

