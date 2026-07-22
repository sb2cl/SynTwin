%% Estimate_Lib6_ALL_reduced_model_set_dma
% Estimate the shared B0034 parameter across fixed d_mA values.
%
% DESCRIPTION
%   For each d_mA value, estimates the shared B0034 intrinsic initiation
%   capacity, kappa^0, using all six Lib6 constructs.
%
%   Promoter rates and effective plasmid copy numbers are inherited from the
%   Lib24 leave-one-out reduced-model tensor. Run num_run uses the matching
%   inherited Monte Carlo samples of Omega and Gene_cn.
%
% FIXED VALUES
%   - d_mA = 0.15, 0.17, 0.20, 0.22, and 0.25 min^-1;
%   - rho^0 = 0.02.
%
% OUTPUTS
%   Estimated_results/
%       Results_BADS_Lib6_ALL_reduced_<Use_mean>_dma_<tag>.mat
%
%   The Estimated_results directory is distributed. Complete Results tensors
%   are regenerated with Generate_Results_Lib6_ALL_reduced_model_set_dma.m.
%
% USAGE
%   Estimate_Lib6_ALL_reduced_model_set_dma
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
% Lib6 experimental data are selected from the Lib30-format tensor.
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
num_runs = 70;


indices_plasmids = [1,2];
indices_promoters = [1,2,3];

% NOTICE each combination (i,j,k) below will correspond to a library transcriptional unit (TU)
for i=1:length(indices_plasmids) %Plasmids i=1 high copy (pGreen), i=2-> low copy (pSC101)
  for j=1:length(indices_promoters)
    ParamsData_lib24{i,j}.Gene_cn_MC_samples = Results_Tensor_Lib24_L1O_reduced{i,j,1}.Gene_cn_MC_samples;
    ParamsData_lib24{i,j}.Omega_MC_samples = Results_Tensor_Lib24_L1O_reduced{i,j,1}.Omega_MC_samples;
  end
end

% Verify that each optimization run has matching inherited samples.
available_samples = inf;
for plasmid_idx = 1:numel(indices_plasmids)
    for promoter_idx = 1:numel(indices_promoters)
        available_samples = min(available_samples, ...
            min(numel(ParamsData_lib24{plasmid_idx,promoter_idx}.Gene_cn_MC_samples), ...
                numel(ParamsData_lib24{plasmid_idx,promoter_idx}.Omega_MC_samples)));
    end
end
if num_runs > available_samples
    error('Estimate_Lib6_ALL_reduced_model_set_dma:TooManyRuns', ...
        'num_runs (%d) exceeds the available inherited Lib24 samples (%d).', ...
        num_runs, available_samples);
end

% params: vector of estimated parameters of the model
 % params(1) = RBS_k0_sigma0, RBS IIC of (B0034)
 %x0 = 0.13; % Starting point for RBS_k0_sigma0
 lb = 0.05;   %lower expected bounds for RBS_k0_sigma0,      B0034         
 ub = 0.5;  %upper expected bounds for RBS_k0_sigma0,  B0034
 plb = 0.08;  %Plausible lower bounds bounds for RBS_k0_sigma0, 
 pub =  0.16; %Plausible upper bounds bounds for RBS_k0_sigma0

 % Run BADS, which returns the minimum X and its value FVAL.
set_dma = [0.15,0.17,0.20,0.22,0.25];
delta = 0.2;      
RBS_inv_sigma0 = 0.02;
for case_idx = 1:numel(set_dma)
    model_c.dm_c = set_dma(case_idx); %Mean degradation rate of non-ribosomal mRNA (1/min)
    Results_BADS_B0034_ALL_reduced = cell(num_runs,1);
    parfor num_run = 1:num_runs 
         x0 = plb + (pub-plb).* rand(1,length(lb));
        J=@(parameters) J5w_LogPI_Lib6_ALL_reduced(parameters,model_c,Use_mean,ExpData_Tensor_lib30_micro,ParamsData_lib24, HEM,num_run,RBS_inv_sigma0,delta);
        [params, Jmin_value] = bads(J,x0,lb,ub,plb,pub,[],options);
         Results_BADS_B0034_ALL_reduced{num_run}.results=[params, Jmin_value]; 
    end
    
    % --- Save results to local Estimated_results folder (portable) ---
    results_dir = fullfile(SCRIPT_DIR,'Estimated_results');
    if ~exist(results_dir,'dir')
        mkdir(results_dir);
    end
    file_name = "Results_BADS_Lib6_ALL_reduced_"  + Use_mean + "_dma_" + extractAfter(num2str(model_c.dm_c), '.') + ".mat";
    file_tensor = fullfile(results_dir, file_name);
    save(file_tensor, "Results_BADS_B0034_ALL_reduced", "-mat");
end
