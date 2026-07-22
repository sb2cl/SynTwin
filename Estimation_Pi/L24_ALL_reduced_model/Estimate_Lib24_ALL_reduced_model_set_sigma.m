%% Estimate_Lib24_ALL_reduced_model_set_sigma
% Estimate Lib24 reduced-model parameters across fixed rho^0 values.
%
% PURPOSE
%   Repeats the Lib24 ALL-construct parameter estimation while fixing the
%   RBS inverse scaling parameter, rho^0 = inv_sigma0, at each value in:
%
%       [0.005, 0.01, 0.02, 0.03, 0.05].
%
%   For every fixed rho^0, BADS estimates three promoter transcription
%   rates, four RBS intrinsic initiation capacities (kappa^0), and the
%   effective pGreen copy-number multiplier. pSC101 is fixed at N = 5.
%
% INPUT DATA
%   - Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
%   - Experimental_Data/ExpData_Tensor_lib30_micro.mat
%
% CONFIGURATION
%   - Use_mean: 'Global', 'Instances', or 'Wells'
%   - num_runs: number of independent BADS starts per rho^0 value
%
% OUTPUTS
%   Estimated_results/Results_BADS_Lib24_ALL_reduced_<Use_mean>_<rho_tag>.mat
%
%   The rho tag is the decimal part used by the workflow; for example,
%   rho^0 = 0.02 produces the suffix "02".
%
% DEPENDENCIES
%   - init_SynTwin(...,'bads',true)
%   - J5w_LogPI_Lib24_ALL_reduced
%   - BADS
%   - Parallel Computing Toolbox when parfor is used
%
% USAGE
%   Estimate_Lib24_ALL_reduced_model_set_sigma
%
% NOTES
%   These optimization files are large. They may be regenerated with this
%   script and need not all be included in a compact software distribution.

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
num_runs = 54;

indices_plasmids_lib24 = [1,2];
indices_promoters_lib24 = [1,2,3];
indices_rbss_lib24 = [1,2,4,5];
                
% params: vector of estimated parameters of the model
    % params(1:3) = omega, promoter strengths of (J23106, J23102, J23101) 
    % params(4:7) = RBS_k0_sigma0, RBS IIC of (B0030, B0032, J61100, J61101)
    % params(8) = high copy number multiplying term. So that Cn_high = x(9)*Cn_low, and we know that Cn_low=5
  
    lb = [0.025,0.05,0.05... %lower expected bounds for omega, promoter strengths,
          0.05, 2.5e-3,  2.5e-3,2.5e-3 ...  %lower expected bounds for RBS_k0_sigma0, 
          2.5 ];  %lower expected bound for high copy number multiplying term 

    ub = [0.25,0.45,0.45... %upper expected bounds for omega, promoter strengths,
          2.0, 20e-3,  20e-3,5e-2 ...  %upper expected bounds for RBS_k0_sigma0, 
          8.25 ];  %upper expected bound for high copy number multiplying term 

    plb = [0.04,0.15,0.12... %Plausible lower bounds for omega, promoter strengths,
          0.5, 4e-3,  3e-3,4e-3 ...  %Plausible lower bounds bounds for RBS_k0_sigma0, 
          3.5 ];  %Plausible lower bounds for high copy number multiplying term 

    pub = [0.06,0.25,0.20... %Plausible upper bounds for omega, promoter strengths,
          1.3, 6e-3,  5e-3,10e-3 ...  %Plausible upper bounds bounds for RBS_k0_sigma0, 
          5.5 ];  %Plausible upper bounds for high copy number multiplying term 

    % Run BADS, which returns the minimum X and its value FVAL.
    set_RBS_inv_sigma0 = [0.005,0.01,0.02,0.03,0.05];
    delta = 0.2;      
    for i= 1:numel(set_RBS_inv_sigma0)
        RBS_inv_sigma0 = set_RBS_inv_sigma0(i);
        parfor num_run=1:num_runs 
            x0 = plb + (pub-plb).* rand(1,length(lb));
            J=@(parameters) J5w_LogPI_Lib24_ALL_reduced(parameters,model_c,Use_mean,ExpData_Tensor_lib30_micro,HEM,RBS_inv_sigma0,delta);
            [params, Jmin_value] = bads(J,x0,lb,ub,plb,pub,[],options)
            Results_BADS_Lib24_ALL_reduced{num_run}.results=[params, Jmin_value,x0]; 
        end
        % --- Save results to local Estimated_results folder (portable) ---
        results_dir = fullfile(SCRIPT_DIR,'Estimated_results');
        if ~exist(results_dir,'dir')
            mkdir(results_dir);
        end
        file_name = "Results_BADS_Lib24_ALL_reduced_"  + Use_mean + "_" + extractAfter(num2str(RBS_inv_sigma0), '.') + ".mat";
        file_tensor = fullfile(results_dir, file_name);
        save(file_tensor, "Results_BADS_Lib24_ALL_reduced", "-mat");
    end






