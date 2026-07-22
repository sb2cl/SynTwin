%% Estimate_Lib24_L1O_reduced_model
% Estimate Lib24 reduced-model parameters by leave-one-construct-out analysis.
%
% DESCRIPTION
%   Runs 24 leave-one-out cases. In each case, one Lib24 construct is excluded
%   and the remaining 23 constructs are used to estimate:
%
%   - three promoter transcription rates;
%   - four RBS intrinsic initiation capacities, kappa^0;
%   - the pGreen/pSC101 effective copy-number multiplier.
%
%   The RBS inverse scaling parameter is fixed at rho^0 = 0.02.
%
% OUTPUTS
%   Estimated_results/
%       Results_BADS_Lib24_L1O_reduced_<plasmid><promoter><RBS>_<Use_mean>.mat
%
%   One file is generated for each left-out construct.
%
% DEPENDENCIES
%   - J5w_LogPI_Lib24_L1O_reduced
%   - BADS
%   - SynTwin experimental data and HEM surrogate
%   - Parallel Computing Toolbox for parfor, or replace parfor with for
%
% USAGE
%   Estimate_Lib24_L1O_reduced_model
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
model_c.dm_c  =  0.2; %Mean degradation rate of non-ribosomal mRNA (1/min)
model_c.Em_c =  model_c.lp_c/model_c.le_c*(1- (model_c.lp_c/(model_c.lp_c+model_c.le_c))^(model_c.lp_c/model_c.le_c)) ;  
model_c.WEm_c =  1 + 1/model_c.Em_c;  
model_c.N_pSC101 = 5; %known

options = bads('defaults');
options.UncertaintyHandling = true;
options.Display='final';

% Options are: 
% - 'Global' (global mean for each construct),
% - 'Instances' (mean of each experiment for each construct),
% - 'Wells' (use data of all individual culture wells)

%Use_mean = 'Instances';  
%Use_mean = 'Global';  
Use_mean = 'Wells'; 
num_runs = 10;

indices_plasmids_lib24 = [1,2];
indices_promoters_lib24 = [1,2,3];
indices_rbss_lib24 = [1,2,4,5];

% params: vector of estimated parameters of the model
% params(1:3) = omega, promoter strengths of (J23106, J23102, J23101) 
% params(4:7) = RBS_k0_sigma0, RBS IIC of (B0030, B0032, J61100, J61101)
% params(8) = pGreen/pSC101 copy-number multiplier, with N_pSC101 = 5

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

    RBS_inv_sigma0 = 0.02;
    delta = 0.2;      

for i=1:length(indices_plasmids_lib24) %Plasmids i=1 high copy (pGreen), i=2-> low copy (pSC101)
    for j=1:length(indices_promoters_lib24)
        for k=1:length(indices_rbss_lib24)
                Construct_2_LO = [i,j,k]; % Local Lib24 tensor indices of the left-out construct
                
                % Run BADS, which returns the minimum X and its value FVAL.
                Results_BADS_Lib24_LOOCV_reduced = cell(num_runs,1);
                parfor num_run = 1:num_runs 
                    x0 = plb + (pub-plb).* rand(1,length(lb));
                    J=@(parameters) J5w_LogPI_Lib24_L1O_reduced(parameters,Construct_2_LO,model_c,Use_mean,ExpData_Tensor_lib30_micro,HEM,RBS_inv_sigma0,delta);
                    [params, Jmin_value] = bads(J,x0,lb,ub,plb,pub,[],options)
                     Results_BADS_Lib24_LOOCV_reduced{num_run}.results=[params, Jmin_value,x0]; 
                end
                 % --- Save results to local Estimated_results folder (portable) ---
                results_dir = fullfile(SCRIPT_DIR,'Estimated_results');
                if ~exist(results_dir,'dir')
                    mkdir(results_dir);
                end
                file_name = "Results_BADS_Lib24_L1O_reduced_" + ...
                                num2str(indices_plasmids_lib24(Construct_2_LO(1))) + ...
                                num2str(indices_promoters_lib24(Construct_2_LO(2))) + ...
                                num2str(indices_rbss_lib24(Construct_2_LO(3))) + "_" + ...
                                Use_mean + ".mat";
                file_tensor = fullfile(results_dir, file_name);
                save(file_tensor, "Results_BADS_Lib24_LOOCV_reduced", "-mat");
        end %RBS
    end %PROMOTER
end %PLASMID






