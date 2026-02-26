%% Generate_Results_Lib5_ALL_reduced_model
% SynTwin workflow script: post-processing and compilation of Lib5 (ALL)
% reduced-model results (promoter-only estimation with inherited Ori/RBS parameters).
%
% CONTEXT (Lib5 sublibrary)
%   Lib5 is a 5-TU sublibrary defined by:
%     - 1 plasmid origin: pGreen (single Ori context)
%     - 1 promoter:       J23100 (the promoter being estimated)
%     - 5 RBSs:           [1 2 3 4 5] (five translational contexts)
%
%   In this workflow:
%     - Only the promoter transcription parameter Omega (J23100) is estimated
%       from Lib5 data (ALL scheme).
%     - The effective gene copy number (Gene_cn) and four RBS k0_sigma0 values
%       are inherited from Lib24 reduced-model inference (pGreen context).
%     - The k0_sigma0 of RBS B0034 (not present in Lib24) is inherited from the
%       dedicated Lib6 reduced-model inference.
%
% DESCRIPTION
%   Loads the optimization outputs produced by Estimate_Lib5_ALL_reduced_model
%   (BADS runs estimating Omega for J23100) and compiles them into a unified
%   Results tensor for Lib5:
%     (i)  inherited parameter statistics (Gene_cn; RBS k0_sigma0 for all 5 RBSs),
%     (ii) estimated promoter statistics (Omega for J23100),
%     (iii) experimental Mu(t) and Pi(t) trajectories from Lib5,
%     (iv) digital-twin predictions, sensitivities, and Monte Carlo summaries
%          for downstream analysis and plotting.
%
%   Monte Carlo predictions propagate uncertainty from:
%     - inherited parameters (Gene_cn and RBS k0_sigma0),
%     - the estimated promoter parameter Omega (sampled from its inferred distribution).
%
% INPUTS
%   None (loads results from ./Estimated_results/ and uses local configuration).
%
%   Configuration is set inside the script.
%   The user must set the desired <Use_mean> option. Options are:
%       - 'Global'    (global mean for each construct),
%       - 'Instances' (mean of each experiment for each construct),
%       - 'Wells'     (use data of all individual culture wells).
%
% OUTPUTS (saved to disk)
%   Results_Tensor_Lib5_ALL_reduced_model_<Use_mean>.mat
%     Example: Results_Tensor_Lib5_ALL_reduced_model_Wells.mat
%
% DEPENDENCIES
%   - SynTwin initialization recommended:
%       ROOT = init_SynTwin(...);
%   - Requires the HEM surrogate:
%       Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
%   - Requires Lib5 experimental tensor:
%       Experimental_Data/ExpData_Tensor_lib5_micro.mat
%   - Requires inherited Lib24 reduced-model results (pGreen context):
%       .../Results_Tensor_Lib24_L1O_reduced_Wells.mat
%   - Requires inherited Lib6 reduced-model results (for RBS B0034):
%       .../Results_Tensor_Lib6_L1O_reduced_Wells.mat
%   - Requires local estimation outputs (Omega for J23100):
%       ./Estimated_results/Results_BADS_Lib5_ALL_reduced_<Use_mean>.mat
%
% USAGE
%   Generate_Results_Lib5_ALL_reduced_model
%
% NOTES
%   - Use_mean controls the aggregation level: 'Global', 'Instances', or 'Wells'.
%   - The output file is typically loaded by Show_Results_* scripts.
%   - This script assumes the Lib5 tensor already matches the intended indexing
%     (single Ori/promoter context with 5 RBS entries).

% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin('experimental',true);
% Ensure this script folder is on the path (for local functions in this folder)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);

% --- Load data using portable absolute paths ---
load(SynTwin_path('Generate_HEM','HEM_Surrogate','HEM_Surrogate.mat'));          % loads data of the Host Equivalent Model (HEM) 
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

% Options are: 
% - 'Global' (global mean for each construct),
% - 'Instances' (mean of each experiment for each construct),
% - 'Wells' (use data of all individual culture wells)

%Use_mean = 'Instances';  
%Use_mean = 'Global';  
Use_mean = 'Wells'; 

RBS_inv_sigma0 = 0.02;

% Getting the inherited estimated parameters
Inherited_ParamsData={};
indices_rbss_Lib24_dummy = [1,2,NaN,3,4];
Inherited_ParamsData.Gene_cn_MC_samples = Results_Tensor_Lib6_L1O_reduced{1,1}.Gene_cn_MC_samples;
for i=1:5
    if i==3
        Inherited_ParamsData.RBS{i}.RBS_k0_sigma0_MC_samples = Results_Tensor_Lib6_L1O_reduced{1,1}.RBS_k0_sigma0_MC_samples;
        Inherited_ParamsData.RBS{i}.RBS_k0_sigma0_raw = Results_Tensor_Lib6_L1O_reduced{1,1}.RBS_k0_sigma0_raw;
        Inherited_ParamsData.RBS{i}.RBS_k0_sigma0_mean =  Results_Tensor_Lib6_L1O_reduced{1,1}.RBS_k0_sigma0_mean;
        Inherited_ParamsData.RBS{i}.RBS_k0_sigma0_std =   Results_Tensor_Lib6_L1O_reduced{1,1}.RBS_k0_sigma0_std;
    else
        Inherited_ParamsData.RBS{i}.RBS_k0_sigma0_MC_samples = Results_Tensor_Lib24_L1O_reduced{1,1,indices_rbss_Lib24_dummy(i)}.RBS_k0_sigma0_MC_samples;
        Inherited_ParamsData.RBS{i}.RBS_k0_sigma0_raw = Results_Tensor_Lib24_L1O_reduced{1,1,indices_rbss_Lib24_dummy(i)}.RBS_k0_sigma0_raw;
        Inherited_ParamsData.RBS{i}.RBS_k0_sigma0_mean =  Results_Tensor_Lib24_L1O_reduced{1,1,indices_rbss_Lib24_dummy(i)}.RBS_k0_sigma0_mean;
        Inherited_ParamsData.RBS{i}.RBS_k0_sigma0_std =   Results_Tensor_Lib24_L1O_reduced{1,1,indices_rbss_Lib24_dummy(i)}.RBS_k0_sigma0_std;
    end
end

% Getting the estimated Omega for J23100: 
% --- Load estimation results from the local Estimated_results folder (portable) ---
results_dir = fullfile(SCRIPT_DIR,'Estimated_results');
if ~exist(results_dir,'dir')
     error('Estimated_results folder not found: %s', results_dir);
end
file_name = "Results_BADS_Lib5_ALL_reduced_" +   Use_mean + ".mat";
file_tensor = fullfile(results_dir, file_name);
load(file_tensor, "Results_BADS_J23100_ALL_approx", "-mat"); 

num_runs = length(Results_BADS_J23100_ALL_approx);
Matrix_tempo_All = [];
for i=1:num_runs
   Matrix_tempo_All=[Matrix_tempo_All;Results_BADS_J23100_ALL_approx{i}.results(1,1)];
end
Estimated_omega_J23100.ALL_raw = Matrix_tempo_All;
Estimated_omega_J23100.ALL_mean = mean(Matrix_tempo_All,1);
Estimated_omega_J23100.ALL_std = std(Matrix_tempo_All,0,1);

A=-1;
B=-1e-6;
n_samples_MC = 1000;
Estimated_omega_J23100.Omega_MC_samples = rmvnrnd(Estimated_omega_J23100.ALL_mean, Estimated_omega_J23100.ALL_std.^2, n_samples_MC,A,B);
       
indices_plasmids_lib5 = 1;
num_plasmids_Lib5 = length(indices_plasmids_lib5);
indices_promoters_lib5 = 4;
num_promoters_Lib5 = length(indices_promoters_lib5);
indices_rbss_lib5 = [1,2,3,4,5];
num_rbss_Lib5 = length(indices_rbss_lib5);

% Adds the estimated parameters, experimental and predicted synthesis and their statistics to the structure Results_Tensor_lib5_ALL_approx
for r=1:num_rbss_Lib5
    if strcmp(Use_mean,'Wells' )   
         Results_Tensor_Lib5_ALL_reduced{1,1,r}.Use_mean = 'Wells'; 
    elseif strcmp(Use_mean,'Instances')
         Results_Tensor_Lib5_ALL_reduced{1,1,r}.Use_mean = 'Instances'; 
    elseif strcmp(Use_mean,'Global')
         Results_Tensor_Lib5_ALL_reduced{1,1,r}.Use_mean ='Global'; 
    else
          error('Invalid option for Use_mean. Choose either "Wells", "Instances", or "Global".');
    end
    % General TU information:
  Results_Tensor_Lib5_ALL_reduced{1,1,r}.TU_Ori = ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.TU_Ori;
  Results_Tensor_Lib5_ALL_reduced{1,1,r}.TU_Promoter = ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.TU_Promoter;
  Results_Tensor_Lib5_ALL_reduced{1,1,r}.TU_RBS =  ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.TU_RBS;
  Results_Tensor_Lib5_ALL_reduced{1,1,r}.TU_Bioparts  =  ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.TU_Bioparts;
  Results_Tensor_Lib5_ALL_reduced{1,1,r}.TU_Name =  ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.TU_Name;
  Results_Tensor_Lib5_ALL_reduced{1,1,r}.TU_color_code =  ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.TU_color_code;
    
   %Estimated parameters (Ori and RBSs inherited):
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.Gene_cn_raw = Results_Tensor_Lib24_L1O_reduced{1,1,1}.Gene_cn_raw;
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.Gene_cn_mean = Results_Tensor_Lib24_L1O_reduced{1,1,1}.Gene_cn_mean;
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.Gene_cn_std =  Results_Tensor_Lib24_L1O_reduced{1,1,1}.Gene_cn_std;
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.Omega_raw = Estimated_omega_J23100.ALL_raw; 
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.Omega_mean =  Estimated_omega_J23100.ALL_mean; 
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.Omega_std =  Estimated_omega_J23100.ALL_std; 
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.RBS_k0_sigma0_raw = Inherited_ParamsData.RBS{r}.RBS_k0_sigma0_raw; 
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.RBS_k0_sigma0_mean = Inherited_ParamsData.RBS{r}.RBS_k0_sigma0_mean; 
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.RBS_k0_sigma0_std = Inherited_ParamsData.RBS{r}.RBS_k0_sigma0_std; 
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.RBS_inv_sigma0_mean = RBS_inv_sigma0; 
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.RBS_inv_sigma0_std = 0; 

   % TU experimental data:
      % GLOBAL INFORMATION:
   Results_Tensor_Lib5_ALL_reduced{1,1,r}.Mu_mumax_pmax_global_mean = ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.Mu_mumax_pmax_global_mean;
   Results_Tensor_Lib5_ALL_reduced{1,1,r}.Mu_mumax_pmax_global_std = ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.Mu_mumax_pmax_global_std;
   Results_Tensor_Lib5_ALL_reduced{1,1,r}.Pi_mumax_pmax_global_mean = ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.Pi_mumax_pmax_global_mean;
   Results_Tensor_Lib5_ALL_reduced{1,1,r}.Pi_mumax_pmax_global_std = ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.Pi_mumax_pmax_global_std;

    % INSTANCES INFORMATION:
   % go through all the instances:
   instance_count = length(ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.Instances);
   Results_Tensor_Lib5_ALL_reduced{1,1,r}.Instances = cell(instance_count, 1);
   for w = 1:instance_count
     Results_Tensor_Lib5_ALL_reduced{1,1,r}.Instances{w}.Mu_mumax_pmax_instance_mean = ...
            ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.Instances{w}.Mu_mumax_pmax_instance_mean;
     Results_Tensor_Lib5_ALL_reduced{1,1,r}.Instances{w}.Mu_mumax_pmax_instance_std = ...
            ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.Instances{w}.Mu_mumax_pmax_instance_std;
     Results_Tensor_Lib5_ALL_reduced{1,1,r}.Instances{w}.Pi_mumax_pmax_instance_mean = ...
            ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.Instances{w}.Pi_mumax_pmax_instance_mean;
     Results_Tensor_Lib5_ALL_reduced{1,1,r}.Instances{w}.Pi_mumax_pmax_instance_std = ...
            ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.Instances{w}.Pi_mumax_pmax_instance_std;
   end % for each instance

   % INDIVIDUAL WELLS INFORMATION:
   % go through all the instances and then all the wells for each instance:
  for w = 1:instance_count
      no_wells = length(ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.Instances{w}.Wells);
      Results_Tensor_Lib5_ALL_reduced{1,1,r}.Instances{w}.Wells = cell(no_wells, 1);
      for ww = 1:no_wells
          Results_Tensor_Lib5_ALL_reduced{1,1,r}.Instances{w}.Wells{ww}.Mu_mumax_pmax = ...
              ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.Instances{w}.Wells{ww}.Mu_mumax_pmax;
          Results_Tensor_Lib5_ALL_reduced{1,1,r}.Instances{w}.Wells{ww}.Pi_mumax_pmax = ...
              ExpData_Tensor_lib5_micro{1,1,indices_rbss_lib5(r)}.Instances{w}.Wells{ww}.Pi_mumax_pmax;
      end % for each well
  end % for each instance

 %Get the sensitivities of the synthesis rate w.r.t. the estimated parameters at the experimental values of growth rate:
              % GLOBAL 
              Mu_vector = Results_Tensor_Lib5_ALL_reduced{1,1,r}.Mu_mumax_pmax_global_mean;
              Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,Results_Tensor_Lib5_ALL_reduced{1,1,r}.RBS_k0_sigma0_mean,...
                                           RBS_inv_sigma0,Results_Tensor_Lib5_ALL_reduced{1,1,r}.Omega_mean,Results_Tensor_Lib5_ALL_reduced{1,1,r}.Gene_cn_mean);
              Results_Tensor_Lib5_ALL_reduced{1,1,r}.S_Pi_NA_global_values = Results_Tensor_Lib5_ALL_reduced{1,1,r}.Pi_mumax_pmax_global_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_NA_values;
              Results_Tensor_Lib5_ALL_reduced{1,1,r}.S_Pi_Omega_global_values = Results_Tensor_Lib5_ALL_reduced{1,1,r}.Pi_mumax_pmax_global_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_Omega_values;
              Results_Tensor_Lib5_ALL_reduced{1,1,r}.S_Pi_RBS_k0_sigma0_global_values = Results_Tensor_Lib5_ALL_reduced{1,1,r}.Pi_mumax_pmax_global_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_k0_sigma0_values;
              Results_Tensor_Lib5_ALL_reduced{1,1,r}.S_Pi_RBS_inv_sigma0_global_values = Results_Tensor_Lib5_ALL_reduced{1,1,r}.Pi_mumax_pmax_global_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_inv_sigma0_values;
             % INSTANCES
             p=1;
             q=1;
             instance_count = length(Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances);
             for i = 1:instance_count
               if ~(isempty(Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Mu_mumax_pmax_instance_mean)) %some experiments can be empty if all wells of a construct are outliers
                 Mu_vector = Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Mu_mumax_pmax_instance_mean;
                 Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,Results_Tensor_Lib5_ALL_reduced{p,q,r}.RBS_k0_sigma0_mean,...
                                           RBS_inv_sigma0,Results_Tensor_Lib5_ALL_reduced{p,q,r}.Omega_mean,Results_Tensor_Lib5_ALL_reduced{p,q,r}.Gene_cn_mean);
                 Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.S_Pi_NA_instance_values = Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Pi_mumax_pmax_instance_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_NA_values;
                 Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.S_Pi_Omega_instance_values = Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Pi_mumax_pmax_instance_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_Omega_values;
                 Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.S_Pi_RBS_k0_sigma0_instance_values = Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Pi_mumax_pmax_instance_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_k0_sigma0_values;
                 Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.S_Pi_RBS_inv_sigma0_instance_values = Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Pi_mumax_pmax_instance_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_inv_sigma0_values;
               end % if not empty
            end % for each instance
            % WELLS
            for i = 1:instance_count
               if ~(isempty(Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Mu_mumax_pmax_instance_mean)) %some experiments can be empty if all wells of a construct are outliers
                  no_wells = length(Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Wells);
                  for j = 1:no_wells
                      Mu_vector = Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.Mu_mumax_pmax;
                      Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,Results_Tensor_Lib5_ALL_reduced{p,q,r}.RBS_k0_sigma0_mean,...
                                           RBS_inv_sigma0,Results_Tensor_Lib5_ALL_reduced{p,q,r}.Omega_mean,Results_Tensor_Lib5_ALL_reduced{p,q,r}.Gene_cn_mean);
                      Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.S_Pi_NA_well_values = Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.Pi_mumax_pmax.*Synthesis_predictions.Relative_sensitivities.S_Pi_NA_values;
                      Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.S_Pi_Omega_well_values = Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.Pi_mumax_pmax.*Synthesis_predictions.Relative_sensitivities.S_Pi_Omega_values;
                      Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.S_Pi_RBS_k0_sigma0_well_values = Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.Pi_mumax_pmax.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_k0_sigma0_values; 
                      Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.S_Pi_RBS_inv_sigma0_well_values = Results_Tensor_Lib5_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.Pi_mumax_pmax.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_inv_sigma0_values;
                  end % for each well
               end % if not empty
           end % for each instance

     % Get synthesis predictions for a range of values of growth rate:
     Mu_vector =  0.005:0.003:0.032;
     Results_Tensor_Lib5_ALL_reduced{1,1,r}.Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,Results_Tensor_Lib5_ALL_reduced{1,1,r}.RBS_k0_sigma0_mean,...
                                                   Results_Tensor_Lib5_ALL_reduced{1,1,r}.RBS_inv_sigma0_mean,Results_Tensor_Lib5_ALL_reduced{1,1,r}.Omega_mean,Results_Tensor_Lib5_ALL_reduced{1,1,r}.Gene_cn_mean);
     Results_Tensor_Lib5_ALL_reduced{1,1,r}.Mu_slices_values = Mu_vector;
    
    % MONTE CARLO predictions for a range of values of growth rate:
    % First generate distributions for the parameters
    % To avoid values less or equal to zero, uses Tim Benham (2025). Truncated multivariate normal (https://www.mathworks.com/matlabcentral/fileexchange/34402-truncated-multivariate-normal), MATLAB Central File Exchange.
     A=-1;
     B=-1e-6;
     Results_Tensor_Lib5_ALL_reduced{1,1,r}.Gene_cn_MC_samples = Inherited_ParamsData.Gene_cn_MC_samples;
     Results_Tensor_Lib5_ALL_reduced{1,1,r}.Omega_MC_samples = Estimated_omega_J23100.Omega_MC_samples;
     Results_Tensor_Lib5_ALL_reduced{1,1,r}.RBS_k0_sigma0_MC_samples = Inherited_ParamsData.RBS{r}.RBS_k0_sigma0_MC_samples;
     MC_indexes = 1:1:n_samples_MC; 
     for j=1:n_samples_MC
         sample_gen_cn = Results_Tensor_Lib5_ALL_reduced{1,1,r}.Gene_cn_MC_samples(MC_indexes(j));
         sample_Omega = Results_Tensor_Lib5_ALL_reduced{1,1,r}.Omega_MC_samples(MC_indexes(j));
         sample_RBS_k0_sigma0 = Results_Tensor_Lib5_ALL_reduced{1,1,r}.RBS_k0_sigma0_MC_samples(MC_indexes(j));
         RBS_inv_sigma0 = 0.02;
         Results_Tensor_Lib5_ALL_reduced{1,1,r}.MC_samples(j).Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,...
                                                      sample_RBS_k0_sigma0,RBS_inv_sigma0,sample_Omega,sample_gen_cn);
     end % for MC samples

    % Next, for each value of growth rate, we can get the probability density function of the predicted synthesis rate and get statistics from it
   MC_Pi_pred_mu = NaN*ones(n_samples_MC,length(Mu_vector)); %each column is a value of growth rate, each row is a MC sample
   MC_Pi_pred_mu_approx = NaN*ones(n_samples_MC,length(Mu_vector)); 
   for j=1:n_samples_MC
       MC_Pi_pred_mu(j,:) = Results_Tensor_Lib5_ALL_reduced{1,1,r}.MC_samples(j).Synthesis_predictions.Pi_pred_values';
   end
 
   for m=1:length(Mu_vector)
    pd = fitdist(MC_Pi_pred_mu(:,m),'kernel','Kernel','epanechnikov');
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.MC_mu_slices(m).Pi_pdf = pd;
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.MC_mu_slices(m).Mu_slice = Mu_vector(m);
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.MC_mu_slices(m).Pi_pred_mean = mean(pd);
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.MC_mu_slices(m).Pi_pred_std = std(pd);
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.MC_mu_slices(m).Pi_pred_q50 = icdf(pd,0.5);
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.MC_mu_slices(m).Pi_pred_q25 = icdf(pd,0.25);
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.MC_mu_slices(m).Pi_pred_q75 = icdf(pd,0.75);
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.MC_mu_slices(m).Pi_pred_q2p5 = icdf(pd,0.025);
    Results_Tensor_Lib5_ALL_reduced{1,1,r}.MC_mu_slices(m).Pi_pred_q97p5 = icdf(pd,0.975);
   end

end %rbss

% SAVE the results
% --- Save generated results in the same folder as this script (portable) ---
file_name  = "Results_Tensor_Lib5_ALL_reduced_" + Use_mean + ".mat";
file_tensor = fullfile(SCRIPT_DIR, file_name);
save(file_tensor, "Results_Tensor_Lib5_ALL_reduced", "-mat");

 
