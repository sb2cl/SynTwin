%% Estimate_Lib5_L1O_reduced_model
% Estimate J23100 promoter strength by leave-one-construct-out analysis.
%
% DESCRIPTION
%   Runs five folds. In each fold, one Lib5 RBS-defined construct is excluded
%   and the J23100 promoter parameter, Omega, is estimated from the remaining
%   four constructs.
%
%   Effective pGreen copy number and RBS kappa^0 values are inherited from the
%   Lib24 and Lib6 leave-one-out reduced-model tensors. Each optimization run
%   uses one consistent inherited Monte Carlo realization across all constructs.
%
% FIXED VALUE
%   rho^0 = 0.02.
%
% OUTPUTS
%   Estimated_results/
%       Results_BADS_Lib5_L1O_reduced_14<RBS>_<Use_mean>.mat
%
%   One file is generated for each left-out RBS-defined construct.
%
% DEPENDENCIES
%   - J5w_LogPI_Lib5_L1O_reduced
%   - Results_Tensor_Lib24_L1O_reduced_Wells.mat
%   - Results_Tensor_Lib6_L1O_reduced_Wells.mat
%   - BADS, HEM surrogate, and processed Lib5 experimental data
%
% USAGE
%   Estimate_Lib5_L1O_reduced_model
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
num_runs = 56;

Inherited_ParamsData = struct();
indices_rbss_Lib24 = [1,2,NaN,3,4];

available_samples = numel(Results_Tensor_Lib6_L1O_reduced{1,1}.Gene_cn_MC_samples);
available_samples = min(available_samples, ...
    numel(Results_Tensor_Lib6_L1O_reduced{1,1}.RBS_k0_sigma0_MC_samples));
for rbs_idx = [1,2,4,5]
    source_idx = indices_rbss_Lib24(rbs_idx);
    available_samples = min(available_samples, ...
        numel(Results_Tensor_Lib24_L1O_reduced{1,1,source_idx}.RBS_k0_sigma0_MC_samples));
end
if num_runs > available_samples
    error('Estimate_Lib5_L1O_reduced_model:TooManyRuns', ...
        'num_runs (%d) exceeds the common inherited sample count (%d).', ...
        num_runs, available_samples);
end

index_MC_pars = randi(available_samples,num_runs,1);
Inherited_ParamsData.Gene_cn_MC_samples = ...
    Results_Tensor_Lib6_L1O_reduced{1,1}.Gene_cn_MC_samples(index_MC_pars,1);

for rbs_idx = 1:5
    if rbs_idx == 3
        source_samples = Results_Tensor_Lib6_L1O_reduced{1,1}.RBS_k0_sigma0_MC_samples;
    else
        source_idx = indices_rbss_Lib24(rbs_idx);
        source_samples = Results_Tensor_Lib24_L1O_reduced{1,1,source_idx}.RBS_k0_sigma0_MC_samples;
    end
    Inherited_ParamsData.RBS{rbs_idx}.RBS_k0_sigma0_MC_samples = ...
        source_samples(index_MC_pars,1);
end

% params: vector of estimated parameters of the model
 % params(1) = Omega, Promoter strength of (J23100)
% x0 = 0.15; % Starting point f
 lb = 0.02;   %lower expected bounds   Promoter Omega,     J23100          
 ub = 1.5;  %upper expected bounds  Promoter Omega,     J23100    
 plb = 0.10;  %Plausible lower bounds 
 pub =  0.25; %Plausible upper bounds 


     % Run BADS, which returns the minimum X and its value FVAL.
% Lib5 bioparts indices:
indices_rbss = [1,2,3,4,5];  % WE MAY SKIP RBS B0034 (no. 3) if we want to demonstrate that we will be capable to predict its characteristic curve Pi-mu
num_rbss = length(indices_rbss);
RBS_inv_sigma0 = 0.02;
delta = 0.2;  
for j=1:num_rbss 
    Construct_2_LO = j; % Local Lib5 RBS index of the left-out construct
    Results_BADS_J23100_L1O_reduced = cell(num_runs,1);
    parfor num_run = 1:num_runs 
         x0 = plb + (pub-plb).* rand(1,length(lb));
        J=@(parameters) J5w_LogPI_Lib5_L1O_reduced(parameters,Construct_2_LO,model_c,Use_mean,ExpData_Tensor_lib5_micro,Inherited_ParamsData, HEM,num_run,RBS_inv_sigma0,delta);
        [params, Jmin_value] = bads(J,x0,lb,ub,plb,pub,[],options);
        Results_BADS_J23100_L1O_reduced{num_run}.results = [params, Jmin_value];
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
    save(file_tensor, "Results_BADS_J23100_L1O_reduced", "-mat");

end 



