%% Estimate_Lib6_L1O_reduced_model
% SynTwin workflow script: parameter estimation for the L6 sublibrary using
% Leave-One-Out (L1O) cross-validation under the reduced model
% (estimate a single shared-RBS k0_sigma0; fix inv_sigma0).
%
% CONTEXT (L6 sublibrary)
%   L6 is a 6-TU sublibrary defined by:
%     - 2 plasmid origins (Ori indices: [1 2])
%     - 3 promoters       (Promoter indices: [1 2 3])
%     - 1 shared RBS      (RBS index: [3], e.g. B0034)
%
%   In this workflow, only the shared RBS intrinsic initiation capacity
%   k0_sigma0 is estimated. All Ori- and promoter-related parameters are
%   inherited from Lib24 inference results (Omega and Gene_cn), while known
%   constants (e.g. pSC101 baseline copy number) are provided via model_c.
%
% DESCRIPTION
%   Runs leave-one-out (L1O) cross-validation across the 6 constructs in L6.
%   For each left-out construct (defined by Construct_2_LO = [Ori, Prom, RBS]),
%   the script performs a user-defined number of optimization runs (<num_runs>)
%   using BADS, minimizing the library-scale mismatch between experimental
%   synthesis-rate data and SynTwin digital-twin predictions (Pi) computed on
%   the remaining (training) constructs.
%
%   Importantly, uncertainty from inherited Lib24 parameters is propagated by
%   using the run index num_run to select a corresponding Monte Carlo sample of:
%       - Omega (context-dependent promoter/transcription term), and
%       - Gene_cn (effective gene copy number),
%   which are kept fixed during that optimization run.
%
% MODEL SETUP (reduced model for RBSs)
%   - Estimated: k0_sigma0 (intrinsic initiation capacity of the shared RBS)
%   - Fixed:     inv_sigma0 (sensitivity-related parameter; fixed in the cost function)
%
% INPUTS
%   None (configuration is set inside the script).
%   The user must set:
%     - Use_mean: aggregation level ('Global', 'Instances', or 'Wells')
%     - num_runs: number of Monte Carlo / optimization runs per L1O fold
%
% OUTPUTS (saved to disk)
%   ./Estimated_results/Results_BADS_Lib6_L1O_reduced_<Ori><Prom><RBS>_<Use_mean>.mat
%     One file per left-out construct (L1O fold).
%     Each file stores a cell array of per-run results, e.g.:
%       Results_BADS_*{num_run}.results = [k0_sigma0, Jmin]
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
%   - Inherited Lib24 parameter samples (Omega and Gene_cn), loaded from:
%       Results_Tensor_Lib24_L1O_reduced_*.mat
%   - Uses: J4_LogPI_Lib6_L1O_reduced (objective function)
%   - Requires BADS (bundled in SynTwin distribution) and optional Parallel TB.
%
% USAGE
%   Estimate_Lib6_L1O_reduced_model
%
% NOTES
%   - L6 experimental measurements are retrieved from the Lib30 tensor-format
%     container by selecting the appropriate indices (Ori=[1 2], Prom=[1 2 3], RBS=3).
%   - If Parallel Toolbox is unavailable, replace parfor with for.
%   - This script quantifies the predictive generalization of an inherited
%     (Ori/Promoter) context while identifying a shared-RBS parameter in L6.

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
num_runs = 42; 

indices_plasmids = [1,2];
indices_promoters = [1,2,3];
indices_rbss = [3];%Index of RBS in ExpData_Tensor_lib30_micro

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

for i=1:length(indices_plasmids) %Plasmids i=1 high copy (pGreen), i=2-> low copy (pSC101)
    for j=1:length(indices_promoters)
        for k=1:length(indices_rbss)
            Construct_2_LO = [i,j,k] %Construct to leave out: CN (1 to 2), Promoter (1 to 3),  RBS (3=B0034 )  
                            
             % Run BADS, which returns the minimum X and its value FVAL.
            parfor num_run=1:num_runs 
                J=@(parameters) J4_LogPI_Lib6_L1O_reduced(parameters,Construct_2_LO,model_c,Use_mean,ExpData_Tensor_lib30_micro,ParamsData_lib24, HEM,num_run);
                [params, Jmin_value] = bads(J,x0,lb,ub,plb,pub,[],options);
                 Results_BADS_Lib6_LOOCV_approx{num_run}.results=[params, Jmin_value]; 
            end
         % --- Save results to local Estimated_results folder (portable) ---
                results_dir = fullfile(SCRIPT_DIR,'Estimated_results');
                if ~exist(results_dir,'dir')
                    mkdir(results_dir);
                end
                file_name = "Results_BADS_Lib6_L1O_reduced_" + ...
                                num2str(indices_plasmids(Construct_2_LO(1))) + ...
                                num2str(indices_promoters(Construct_2_LO(2))) + ...
                                num2str(indices_rbss(Construct_2_LO(3))) + "_" + ...
                                Use_mean + ".mat";
                file_tensor = fullfile(results_dir, file_name);
                save(file_tensor, "Results_BADS_B0034_LOOCV_approx", "-mat");
        end %rbss
    end % promoters
end %plasmid





