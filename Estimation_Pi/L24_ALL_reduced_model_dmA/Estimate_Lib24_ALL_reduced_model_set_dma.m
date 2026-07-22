%% Estimate_Lib24_ALL_reduced_model_set_dma
% Estimate Lib24 reduced-model parameters across fixed mRNA degradation rates.
%
% DESCRIPTION
%   Runs repeated BADS optimizations for all 24 Lib24 constructs while fixing
%   the non-ribosomal mRNA degradation rate, d_mA, at several values. The RBS
%   inverse scaling parameter is fixed at rho^0 = 0.02.
%
% ESTIMATED PARAMETERS
%   - promoter transcription rates for J23106, J23102, and J23101;
%   - RBS intrinsic initiation capacities for B0030, B0032, J61100, and J61101;
%   - the pGreen/pSC101 effective copy-number multiplier.
%
% FIXED VALUES
%   - pSC101 effective copy number: 5;
%   - RBS inverse scaling parameter: rho^0 = 0.02;
%   - d_mA values: 0.15, 0.17, 0.20, 0.22, and 0.25 min^-1.
%
% OUTPUTS
%   Estimated_results/
%       Results_BADS_Lib24_ALL_reduced_<Use_mean>_dma_<tag>.mat
%
%   The distributed repository contains these optimization outputs. The larger
%   compiled Results_Tensor files are regenerated with
%   Generate_Results_Lib24_ALL_reduced_model_set_dma.m.
%
% REQUIREMENTS
%   - SynTwin initialized with experimental data and BADS support;
%   - BADS;
%   - Parallel Computing Toolbox for parfor, or replace parfor with for.
%
% USAGE
%   Estimate_Lib24_ALL_reduced_model_set_dma
%
% See README.md for the complete workflow.

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
num_runs = 14;

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
    set_dma = [0.15,0.17,0.20,0.22,0.25];
    delta = 0.2;      
    RBS_inv_sigma0 = 0.02;
    for i = 1:numel(set_dma)
        model_c.dm_c = set_dma(i); % Mean degradation rate of non-ribosomal mRNA (1/min)
        Results_BADS_Lib24_ALL_reduced = cell(num_runs,1);
        parfor num_run = 1:num_runs 
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
        file_name = "Results_BADS_Lib24_ALL_reduced_"  + Use_mean + "_dma_" + extractAfter(num2str(model_c.dm_c), '.') + ".mat";
        file_tensor = fullfile(results_dir, file_name);
        save(file_tensor, "Results_BADS_Lib24_ALL_reduced", "-mat");
    end






