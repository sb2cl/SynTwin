%% Estimate_Lib30_ALL_full_model
% SynTwin workflow script: parameter estimation for Lib30 using ALL TUs
% under the complete model for RBSs (estimate both k0_sigma0 and  inv_sigma0).
%
% DESCRIPTION
%   Runs a user-defined <num_runs> of estimations using ALL constructs in Lib30.
%   Parameters are estimated by minimizing the mismatch between experimental
%   synthesis-rate data and SynTwin digital-twin predictions (Pi) using 
%   the BADS optimizer.
%
% MODEL SETUP (complete model for RBSs)
%   - Estimated: k0_sigma0 (intrinsic initiation capacity)
%                inv_sigma0 (sensitivity-related parameter)
%
% INPUTS
%   None (configuration is set inside the script).
%   The user must set the desired <Use_mean> option. Options are: 
%       - 'Global' (global mean for each construct),
%       - 'Instances' (mean of each experiment for each construct),
%       - 'Wells' (use data of all individual culture wells)
%   The user must set the desired <num_runs> option.
%
% OUTPUTS (saved to disk)
%   ./Estimated_results/Results_BADS_Lib30_ALL_complete_<Use_mean>.mat
%
% DEPENDENCIES
%   - SynTwin must be initialized before running this script:
%       ROOT = init_SynTwin(...);
%   - This folder must be on the MATLAB path (if running from SynTwin root):
%       addpath(fileparts(mfilename('fullpath')));
%   - Uses: J4_LogPI_Lib24_ALL_full_model (objective function)
%   - This distribution requires BADS (bundled in SynTwin distribution) and optional Parallel TB.
%
% USAGE
%   Estimate_Lib30_ALL_full_model
%
% NOTES
%   - The experimental data is loaded from ExpData_Tensor_Lib30_micro.mat.
%   - If Parallel Toolbox is unavailable, replace parfor with for.

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
load(SynTwin_path('Experimental_Data','ExpData_Tensor_lib30_micro.mat'));        % loads data of Lib30: ExpData_Tensor_lib30_micro

% Consider an exogenous circuit. We use a Transcriptional Unit (TU) expressing GFP:
model_c.lp_c = 240; %Length of GFP protein (aa)     
model_c.le_c = model_c.lp_c^0.097/0.0703; 	%Ribosome occupancy length (aa) %lea=la^0.097/0.0703
model_c.dm_c  =  0.2; %Mean degradation rate of non-ribosomal mRNA (1/min)
model_c.Em_c =  model_c.lp_c/model_c.le_c*(1- (model_c.lp_c/(model_c.lp_c+model_c.le_c))^(model_c.lp_c/model_c.le_c)) ;  
model_c.WEm_c =  1 + 1/model_c.Em_c;  
model_c.N_pSC101 = 5; %known

options = bads('defaults');
options.Display='final';
%options.Display='iter';

% Options are: 
% - 'Global' (global mean for each construct),
% - 'Instances' (mean of each experiment for each construct),
% - 'Wells' (use data of all individual culture wells)

%Use_mean = 'Instances';  
Use_mean = 'Global';  
%Use_mean = 'Wells'; 
num_runs = 84;

indices_plasmids_lib30 = [1,2];
indices_promoters_lib30 = [1,2,3];
indices_rbss_lib30 = [1,2,3,4,5];
                                
% params: vector of estimated parameters of the model
    % params(1:3) = omega, promoter strengths of (J23106, J23102, J23101) 
    % params(4:8) = RBS_k0_sigma0, RBS strengths over sensitivities (translational efficiency)  of (B0030, B0032, B0034, J61100, J61101)
    % params(9:13) = RBS_inv_sigma0, RBS inverse sensitivities of (B0030, B0032, B0034, J61100, J61101)
    % params(14) = high copy number multiplying term. So that Cn_high = x(9)*Cn_low, and we know that Cn_low=5
    
    x0 = [0.05,0.19, 0.16... % Starting point for omega 
          0.8, 0.005, 0.13, 0.004, 0.008 ...   % Starting point for RBS_k0_sigma0
          0.02, 0.02, 0.02, 0.02, 0.02 ...   % Starting point for RBS_inv_sigma0
          4.4 ];                 % Starting point for high copy number multiplying term
    lb = [0.025,0.05,0.05... %lower expected bounds for omega, promoter strengths,
          0.05, 2.5e-3, 0.05, 2.5e-3,2.5e-3 ...  %lower expected bounds for RBS_k0_sigma0, 
          0.0025, 0.0025, 0.0025, 0.0025, 0.0025 ...   % lower expected boundsfor RBS_inv_sigma0
          2.5 ];  %lower expected bound for high copy number multiplying term 

    ub = [0.15,0.35,0.35... %upper expected bounds for omega, promoter strengths,
          2.0, 10e-3, 0.35, 8e-3,2e-2 ...  %upper expected bounds for RBS_k0_sigma0, 
          0.05, 0.05, 0.05, 0.05, 0.05 ...   % upper expected bounds for RBS_inv_sigma0
          8.25 ];  %upper expected bound for high copy number multiplying term 

    plb = [0.04,0.15,0.12... %Plausible lower bounds for omega, promoter strengths,
          0.5, 4e-3, 0.12, 3e-3,4e-3 ...  %Plausible lower bounds bounds for RBS_k0_sigma0, 
           0.015, 0.015, 0.015, 0.015, 0.015 ...   %Plausible lower bounds for RBS_inv_sigma0
          3.5 ];  %Plausible lower bounds for high copy number multiplying term 

    pub = [0.06,0.25,0.20... %Plausible upper bounds for omega, promoter strengths,
          1.3, 6e-3, 0.16, 5e-3,10e-3 ...  %Plausible upper bounds bounds for RBS_k0_sigma0, 
          0.025, 0.025, 0.025, 0.025, 0.025 ...   %Plausible upper bounds for RBS_inv_sigma0
          5.5 ];  %Plausible upper bounds for high copy number multiplying term 
                
    % Run BADS, which returns the minimum X and its value FVAL.

parfor num_run=1:num_runs 
    J=@(parameters) J4_LogPI_Lib30_ALL_full_model(parameters,model_c,Use_mean,ExpData_Tensor_lib30_micro,HEM);
    [params, Jmin_value] = bads(J,x0,lb,ub,plb,pub,[],options)
     Results_BADS_Lib30_ALL_complete_v2{num_run}.results=[params, Jmin_value]; 
end

% --- Save results to local Estimated_results folder (portable) ---
    results_dir = fullfile(SCRIPT_DIR,'Estimated_results');
    if ~exist(results_dir,'dir')
        mkdir(results_dir);
    end
    file_name = "Results_BADS_Lib30_ALL_complete_"  + Use_mean + ".mat";
    file_tensor = fullfile(results_dir, file_name);
    save(file_tensor, "Results_BADS_Lib30_ALL_complete_v2", "-mat");




