%% Estimate_Lib5_L1O_reduced_model
% SynTwin workflow script: parameter estimation for Lib5 using L1O cross-validation
% under a reduced setting (estimate promoter strength Omega; inherit RBS/Gene_cn).
%
% CONTEXT (Lib5 sublibrary)
%   Lib5 is a 5-TU sublibrary defined by:
%     - a single Ori context (pGreen; Ori index fixed to 1),
%     - a single promoter to estimate (e.g., J23100), and
%     - five RBS variants (RBS indices: [1 2 3 4 5]).
%
%   In this workflow, ONLY the promoter strength parameter Omega is estimated.
%   All RBS-related parameters (k0_sigma0) and effective gene copy number (Gene_cn)
%   are inherited from prior inference results:
%     - From Lib24 (L1O reduced): RBS k0_sigma0 for the RBSs that exist in Lib24.
%     - From Lib6  (L1O reduced): RBS B0034 (RBS=3) k0_sigma0 and Gene_cn samples.
%
% DESCRIPTION
%   Runs leave-one-out (L1O) cross-validation across Lib5 constructs, where each
%   L1O iteration leaves out one RBS-defined TU (k = 1..5) from the cost function.
%   For each left-out construct, the script runs <num_runs> independent BADS
%   optimizations to estimate the promoter strength Omega by minimizing the
%   mismatch between experimental synthesis-rate data and SynTwin digital-twin
%   predictions (Pi).
%
%   Uncertainty propagation:
%     For each run index num_run, inherited parameters are sampled via a random
%     Monte Carlo index vector (index_MC_pars). The selected inherited sample is
%     used consistently within that run across all included (non-left-out) TUs:
%       - Inherited_ParamsData.Gene_cn_MC_samples(num_run)
%       - Inherited_ParamsData.RBS{i}.RBS_k0_sigma0_MC_samples(num_run)
%     while Omega is optimized.
%
% MODEL SETUP
%   - Estimated: Omega (promoter strength; e.g., J23100)
%   - Inherited: RBS k0_sigma0 (per RBS; includes B0034 from Lib6)
%   - Fixed:     inv_sigma0 (sensitivity-related parameter; fixed inside the cost function)
%
% INPUTS
%   None (configuration is set inside the script).
%   The user must set:
%     - Use_mean: aggregation level ('Global', 'Instances', or 'Wells')
%     - num_runs: number of Monte Carlo / optimization runs (e.g., 50–100)
%
% OUTPUTS (saved to disk)
%   One file per left-out construct in ./Estimated_results/:
%     Results_BADS_Lib5_L1O_reduced_<LO_RBS>_<Use_mean>.mat
%   where <LO_RBS> identifies the left-out RBS index (1..5).
%
%   Each file contains a cell array of per-run results:
%     Results_BADS_J23100_L1O_approx{num_run}.results = [Omega, Jmin]
%
% DEPENDENCIES
%   - SynTwin must be initialized before running this script:
%       ROOT = init_SynTwin(...);
%   - This folder must be on the MATLAB path (if running from SynTwin root):
%       addpath(fileparts(mfilename('fullpath')));
%   - Requires the HEM surrogate:
%       Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
%   - Experimental data container for Lib5:
%       Experimental_Data/ExpData_Tensor_lib5_micro.mat
%   - Inherited parameter sources:
%       * Lib24 (L1O reduced) for RBSs except B0034:
%           L24_L1O_reduced_model/Results_Tensor_Lib24_L1O_reduced_*.mat
%       * Lib6 (L1O reduced) for B0034 and Gene_cn:
%           L6_L1O_reduced_model/Results_Tensor_Lib6_L1O_reduced_*.mat
%   - Uses: J4_LogPI_Lib5_L1O_reduced (objective function)
%   - Requires BADS (bundled in SynTwin distribution) and optional Parallel TB.
%
% USAGE
%   Estimate_Lib5_L1O_reduced_model
%
% NOTES
%   - L1O here leaves out one RBS-defined TU at a time (Construct_2_LO = 1..5).
%   - Each BADS run uses one inherited Monte Carlo sample selected via index_MC_pars,
%     coupling uncertainty propagation to the run index num_run.
%   - If Parallel Toolbox is unavailable, replace parfor with for.
%   - Optionally, the user may skip leaving out B0034 (RBS=3) if demonstrating
%     prediction capability for its Pi–mu characteristic curve is desired.

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
% Getting the data of Lib5 (we ONLY use the experimental data from this
% file):
load(SynTwin_path('Experimental_Data','ExpData_Tensor_lib5_micro.mat'));        % loads data of Lib5: ExpData_Tensor_lib5_micro

% Getting the data of Lib24 L1O reduced model (we ONLY use the estimated parameters for ori pGreen and RBSs from this file): 
load(SynTwin_path('Estimation_Pi/L24_L1O_reduced_model','Results_Tensor_Lib24_L1O_reduced_Wells.mat'));   

% Getting the estimated IIC (K/sigma) of RBS B0034 from Lib6: 
load(SynTwin_path('Estimation_Pi/L6_L1O_reduced_model','Results_Tensor_Lib6_L1O_reduced_Wells.mat'));   

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
num_runs = 56;  %length(ParamsData_Tensor_libX_micro{1,1,1}.Gene_cn_MC_samples) is 1000

Inherited_ParamsData={};
indices_rbss_Lib24_dummy = [1,2,NaN,3,4];
index_MC_pars = randi([1,length(Results_Tensor_Lib6_L1O_reduced{1,1}.Gene_cn_MC_samples)],num_runs,1);
Inherited_ParamsData.Gene_cn_MC_samples = Results_Tensor_Lib6_L1O_reduced{1,1}.Gene_cn_MC_samples(index_MC_pars,1);
for i=1:5
    if i==3
        Inherited_ParamsData.RBS{i}.RBS_k0_sigma0_MC_samples = Results_Tensor_Lib6_L1O_reduced{1,1}.RBS_k0_sigma0_MC_samples(index_MC_pars,1);
    else
         Inherited_ParamsData.RBS{i}.RBS_k0_sigma0_MC_samples = Results_Tensor_Lib24_L1O_reduced{1,1,indices_rbss_Lib24_dummy(i)}.RBS_k0_sigma0_MC_samples(index_MC_pars,:);
    end
end

% params: vector of estimated parameters of the model
 % params(1) = Omega, Promoter strength of (J23100)
 x0 = 0.15; % Starting point f
 lb = 0.02;   %lower expected bounds   Promoter Omega,     J23100          
 ub = 1.5;  %upper expected bounds  Promoter Omega,     J23100    
 plb = 0.10;  %Plausible lower bounds 
 pub =  0.25; %Plausible upper bounds 


     % Run BADS, which returns the minimum X and its value FVAL.
% Lib5 bioparts indices:
indices_rbss = [1,2,3,4,5];  % WE MAY SKIP RBS B0034 (no. 3) if we want to demonstrate that we will be capable to predict its characteristic curve Pi-mu
num_rbss = length(indices_rbss);
for j=1:num_rbss 
    Construct_2_LO = j %Construct to leave out: RBS (indices_rbss(j) )  in experimental data of Lib5
    parfor num_run=1:num_runs 
        J=@(parameters) J4_LogPI_Lib5_L1O_reduced(parameters,Construct_2_LO,model_c,Use_mean,ExpData_Tensor_lib5_micro,Inherited_ParamsData, HEM,num_run);
        [params, Jmin_value] = bads(J,x0,lb,ub,plb,pub,[],options);
         Results_BADS_J23100_L1O_approx{num_run}.results=[params, Jmin_value] 
    end
    % --- Save results to local Estimated_results folder (portable) ---
    results_dir = fullfile(SCRIPT_DIR,'Estimated_results');
    if ~exist(results_dir,'dir')
        mkdir(results_dir);
    end
    file_name = "Results_BADS_Lib5_L1O_reduced_14" + ...
                    num2str(indices_rbss(Construct_2_LO)) + "_" + ...
                    Use_mean + ".mat";
    file_tensor = fullfile(results_dir, file_name);
    save(file_tensor, "Results_BADS_J23100_L1O_approx", "-mat");

end 



