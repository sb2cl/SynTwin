%% Generate_Results_Lib6_ALL_reduced_model
% SynTwin workflow script: post-processing and compilation of L6 (ALL) reduced-model results.
%
% CONTEXT (L6 sublibrary)
%   L6 is a 6-TU sublibrary defined by 2 plasmid origins and 3 promoters,
%   sharing a common RBS. In this workflow:
%     - The shared RBS parameter (k0_sigma0) is estimated in L6 (from BADS runs).
%     - Ori-/promoter-dependent parameters (e.g., Omega and Gene copy number)
%       are inherited from Lib24 inference results.
%     - Experimental measurements for the corresponding 6 TUs are retrieved
%       by indexing a larger experimental tensor container (Lib30-format).
%
% DESCRIPTION
%   Loads the optimization outputs produced by Estimate_L6_ALL_reduced and
%   compiles them into a unified Results tensor for the L6 sublibrary.
%   The generated structure aggregates:
%     (i) inherited parameter statistics (Omega, Gene_cn) from Lib24,
%     (ii) estimated RBS statistics (k0_sigma0) from L6 BADS runs,
%     (iii) experimental Mu(t) and Pi(t) trajectories for the L6 TUs, and
%     (iv) digital-twin predictions, sensitivities, and Monte Carlo summaries
%          required for downstream analysis and plotting.
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
%   Results_Tensor_Lib6_ALL_reduced_model_<Use_mean>.mat
%     Example: Results_Tensor_Lib6_ALL_reduced_model_Wells.mat
%
% DEPENDENCIES
%   - SynTwin initialization recommended:
%       ROOT = init_SynTwin(...);
%   - Requires local estimation outputs:
%       ./Estimated_results/Results_BADS_Lib6_ALL_reduced_<Use_mean>.mat
%   - Requires inherited Lib24 parameter tensor (Omega, Gene_cn):
%       Results_Tensor_Lib24_*.mat
%   - Requires the HEM surrogate:
%       Generate_HEM/HEM_Surrogate/HEM_Surrogate.mat
%   - Requires the experimental tensor container from which L6 is extracted
%     by indexing (Lib30-format container).
%
% USAGE
%   Generate_Results_L6_ALL_reduced_model
%
% NOTES
%   - Use_mean controls the aggregation level: 'Global', 'Instances', or 'Wells'.
%   - The output file is typically loaded by Show_Results_* scripts.
%   - Monte Carlo predictions are generated to propagate uncertainty from both
%     inherited (Lib24) and estimated (L6) parameters.

% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin('experimental',true);
% Ensure this script folder is on the path (for local functions in this folder)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);

% --- Load data using portable absolute paths ---
load(SynTwin_path('Generate_HEM','HEM_Surrogate','HEM_Surrogate.mat'));          % loads data of the Host Equivalent Model (HEM) 
load(SynTwin_path('Experimental_Data','ExpData_Tensor_lib30_micro.mat'));        % loads data of Lib30: ExpData_Tensor_lib30_micro

% Getting the data of Lib24 L1O reduced model (we ONLY use the estimated parameters from this file): 
load(SynTwin_path('Estimation_Pi/L24_L1O_reduced_model','Results_Tensor_Lib24_L1O_reduced_Wells.mat'));     

% Options are: 
% - 'Global' (global mean for each construct),
% - 'Instances' (mean of each experiment for each construct),
% - 'Wells' (use data of all individual culture wells)

%Use_mean = 'Instances';  
%Use_mean = 'Global';  
Use_mean = 'Wells'; 


indices_plasmids_lib30 = [1,2];
num_plasmids_Lib30 = length(indices_plasmids_lib30);
indices_promoters_lib30 = [1,2,3];
num_promoters_Lib30 = length(indices_promoters_lib30);
indices_rbss_lib30 = [1,2,3,4,5];
num_rbss_Lib30 = length(indices_rbss_lib30);
indices_plasmids_lib6 = [1,2];
num_plasmids_lib6 = length(indices_plasmids_lib6);
indices_promoters_lib6 = [1,2,3];
num_promoters_lib6 = length(indices_promoters_lib6);
indices_rbss_lib6 = [3];
num_rbss_lib6 = length(indices_rbss_lib6);

color_blue = hex2rgb('#00008B');
color_grey = hex2rgb('#4b4b4b');
color_grey_boira = hex2rgb('#F2F2F2');
color_grey_neutre = hex2rgb('#CCCCCC');
RBS_colors= ['#A8DADC';'#F4A261';'#B5E48C';'#CDB4DB';'#FFE066'];
Promoter_colors = ['#264653';'#E76F51';'#2A9D8F';'#6A0572'];

%  We use a Transcriptional Unit (TU) expressing GFP:
model_c.lp_c = 240; %Length of GFP protein (aa)     
model_c.le_c = model_c.lp_c^0.097/0.0703; 	%Ribosome occupancy length (aa) %lea=la^0.097/0.0703
model_c.dm_c  =  0.2; %Mean degradation rate of non-ribosomal mRNA (1/min)
model_c.Em_c =  model_c.lp_c/model_c.le_c*(1- (model_c.lp_c/(model_c.lp_c+model_c.le_c))^(model_c.lp_c/model_c.le_c)) ;  
model_c.WEm_c =  1 + 1/model_c.Em_c;  
model_c.N_pSC101 = 5; %known

n_samples_MC = 1000;

RBS_inv_sigma0 = 0.02;

% Loads the estimated parameters for RBS B0034 and generates the structure Estimated_parameters
Matrix_tempo_All =[];
% --- Load estimation results from the local Estimated_results folder (portable) ---
results_dir = fullfile(SCRIPT_DIR,'Estimated_results');
if ~exist(results_dir,'dir')
     error('Estimated_results folder not found: %s', results_dir);
end
file_name = "Results_BADS_Lib6_ALL_reduced_" +   Use_mean + ".mat";
file_tensor = fullfile(results_dir, file_name);
load(file_tensor, "Results_BADS_B0034_ALL_approx", "-mat");
num_runs = length(Results_BADS_B0034_ALL_approx);
for i=1:num_runs
   Matrix_tempo_All=[Matrix_tempo_All;Results_BADS_B0034_ALL_approx{i}.results(1,1)];
end
Estimated_parameters_B0034.ALL_raw = Matrix_tempo_All;
Estimated_parameters_B0034.ALL_mean = mean(Matrix_tempo_All,1);
Estimated_parameters_B0034.ALL_std = std(Matrix_tempo_All,0,1);


% Adds the estimated parameters and their statistics to the structure Results_Tensor_lib6_L1O_approx
for p=1:length(indices_plasmids_lib6)
    for q=1:length(indices_promoters_lib6)
        for r=1:length(indices_rbss_lib6)
            if strcmp(Use_mean,'Wells' )   
                 Results_Tensor_Lib6_ALL_reduced{p,q,r}.Use_mean = 'Wells'; 
            elseif strcmp(Use_mean,'Instances')
                 Results_Tensor_Lib6_ALL_reduced{p,q,r}.Use_mean = 'Instances'; 
            elseif strcmp(Use_mean,'Global')
                 Results_Tensor_Lib6_ALL_reduced{p,q,r}.Use_mean ='Global'; 
            else
                  error('Invalid option for Use_mean. Choose either "Wells", "Instances", or "Global".');
            end
            % General TU information:
          Results_Tensor_Lib6_ALL_reduced{p,q,r}.TU_Ori = ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.TU_Ori;
          Results_Tensor_Lib6_ALL_reduced{p,q,r}.TU_Promoter = ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.TU_Promoter;
          Results_Tensor_Lib6_ALL_reduced{p,q,r}.TU_RBS =  ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.TU_RBS;
          Results_Tensor_Lib6_ALL_reduced{p,q,r}.TU_Bioparts  =  ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.TU_Bioparts;
          Results_Tensor_Lib6_ALL_reduced{p,q,r}.TU_Name =  ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.TU_Name;
          Results_Tensor_Lib6_ALL_reduced{p,q,r}.TU_color_code =  ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.TU_color_code;
            
           %Estimated parameters (Ori and promoters inherited from L24):
           %Plasmids p=1 high copy (pGreen), p=2-> low copy (pSC101)
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.Gene_cn_raw = Results_Tensor_Lib24_L1O_reduced{p,q,1}.Gene_cn_raw;
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.Gene_cn_mean = Results_Tensor_Lib24_L1O_reduced{p,q,1}.Gene_cn_mean;
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.Gene_cn_std =  Results_Tensor_Lib24_L1O_reduced{p,q,1}.Gene_cn_std;
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.Omega_raw = Results_Tensor_Lib24_L1O_reduced{p,q,1}.Omega_raw; 
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.Omega_mean =  Results_Tensor_Lib24_L1O_reduced{p,q,1}.Omega_mean; 
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.Omega_std =  Results_Tensor_Lib24_L1O_reduced{p,q,1}.Omega_std; 
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.RBS_k0_sigma0_raw = Estimated_parameters_B0034.ALL_raw(:,1); 
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.RBS_k0_sigma0_mean = Estimated_parameters_B0034.ALL_mean(:,1); 
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.RBS_k0_sigma0_std = Estimated_parameters_B0034.ALL_std(:,1); 
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.RBS_inv_sigma0_mean = RBS_inv_sigma0; 
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.RBS_inv_sigma0_std = 0; 

            % Adds the experimental data to the structure
            % Results_Tensor_lib6_ALL_reduced inherited from ExpData_Tensor_lib30_micro:

              % GLOBAL INFORMATION:
           Results_Tensor_Lib6_ALL_reduced{p,q,r}.Mu_mumax_pmax_global_mean = ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.Mu_mumax_pmax_global_mean;
           Results_Tensor_Lib6_ALL_reduced{p,q,r}.Mu_mumax_pmax_global_std = ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.Mu_mumax_pmax_global_std;
           Results_Tensor_Lib6_ALL_reduced{p,q,r}.Pi_mumax_pmax_global_mean = ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.Pi_mumax_pmax_global_mean;
           Results_Tensor_Lib6_ALL_reduced{p,q,r}.Pi_mumax_pmax_global_std = ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.Pi_mumax_pmax_global_std;

            % INSTANCES INFORMATION:
           % go through all the instances:
           instance_count = length(ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.Instances);
           Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances = cell(instance_count, 1);
           for w = 1:instance_count
             Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{w}.Mu_mumax_pmax_instance_mean = ...
                    ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.Instances{w}.Mu_mumax_pmax_instance_mean;
             Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{w}.Mu_mumax_pmax_instance_std = ...
                    ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.Instances{w}.Mu_mumax_pmax_instance_std;
             Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{w}.Pi_mumax_pmax_instance_mean = ...
                    ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.Instances{w}.Pi_mumax_pmax_instance_mean;
             Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{w}.Pi_mumax_pmax_instance_std = ...
                    ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.Instances{w}.Pi_mumax_pmax_instance_std;
           end % for each instance

           % INDIVIDUAL WELLS INFORMATION:
           % go through all the instances and then all the wells for each instance:
          for w = 1:instance_count
              no_wells = length(ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.Instances{w}.Wells);
              Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{w}.Wells = cell(no_wells, 1);
              for ww = 1:no_wells
                  Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{w}.Wells{ww}.Mu_mumax_pmax = ...
                      ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.Instances{w}.Wells{ww}.Mu_mumax_pmax;
                  Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{w}.Wells{ww}.Pi_mumax_pmax = ...
                      ExpData_Tensor_lib30_micro{p,q,indices_rbss_lib6(r)}.Instances{w}.Wells{ww}.Pi_mumax_pmax;
              end % for each well
          end % for each instance

          %Get the sensitivities of the synthesis rate w.r.t. the estimated parameters at the experimental values of growth rate:
              % GLOBAL 
              Mu_vector = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Mu_mumax_pmax_global_mean;
              Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,Results_Tensor_Lib6_ALL_reduced{p,q,r}.RBS_k0_sigma0_mean,...
                                           RBS_inv_sigma0,Results_Tensor_Lib6_ALL_reduced{p,q,r}.Omega_mean,Results_Tensor_Lib6_ALL_reduced{p,q,r}.Gene_cn_mean);
              Results_Tensor_Lib6_ALL_reduced{p,q,r}.S_Pi_NA_global_values = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Pi_mumax_pmax_global_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_NA_values;
              Results_Tensor_Lib6_ALL_reduced{p,q,r}.S_Pi_Omega_global_values = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Pi_mumax_pmax_global_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_Omega_values;
              Results_Tensor_Lib6_ALL_reduced{p,q,r}.S_Pi_RBS_k0_sigma0_global_values = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Pi_mumax_pmax_global_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_k0_sigma0_values;
              Results_Tensor_Lib6_ALL_reduced{p,q,r}.S_Pi_RBS_inv_sigma0_global_values = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Pi_mumax_pmax_global_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_inv_sigma0_values;
             % INSTANCES
             instance_count = length(Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances);
             for i = 1:instance_count
               if ~(isempty(Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Mu_mumax_pmax_instance_mean)) %some experiments can be empty if all wells of a construct are outliers
                 Mu_vector = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Mu_mumax_pmax_instance_mean;
                 Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,Results_Tensor_Lib6_ALL_reduced{p,q,r}.RBS_k0_sigma0_mean,...
                                           RBS_inv_sigma0,Results_Tensor_Lib6_ALL_reduced{p,q,r}.Omega_mean,Results_Tensor_Lib6_ALL_reduced{p,q,r}.Gene_cn_mean);
                 Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.S_Pi_NA_instance_values = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Pi_mumax_pmax_instance_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_NA_values;
                 Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.S_Pi_Omega_instance_values = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Pi_mumax_pmax_instance_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_Omega_values;
                 Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.S_Pi_RBS_k0_sigma0_instance_values = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Pi_mumax_pmax_instance_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_k0_sigma0_values;
                 Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.S_Pi_RBS_inv_sigma0_instance_values = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Pi_mumax_pmax_instance_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_inv_sigma0_values;
               end % if not empty
            end % for each instance
            % WELLS
            for i = 1:instance_count
               if ~(isempty(Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Mu_mumax_pmax_instance_mean)) %some experiments can be empty if all wells of a construct are outliers
                  no_wells = length(Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Wells);
                  for j = 1:no_wells
                      Mu_vector = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.Mu_mumax_pmax;
                      Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,Results_Tensor_Lib6_ALL_reduced{p,q,r}.RBS_k0_sigma0_mean,...
                                           RBS_inv_sigma0,Results_Tensor_Lib6_ALL_reduced{p,q,r}.Omega_mean,Results_Tensor_Lib6_ALL_reduced{p,q,r}.Gene_cn_mean);
                      Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.S_Pi_NA_well_values = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.Pi_mumax_pmax.*Synthesis_predictions.Relative_sensitivities.S_Pi_NA_values;
                      Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.S_Pi_Omega_well_values = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.Pi_mumax_pmax.*Synthesis_predictions.Relative_sensitivities.S_Pi_Omega_values;
                      Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.S_Pi_RBS_k0_sigma0_well_values = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.Pi_mumax_pmax.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_k0_sigma0_values; 
                      Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.S_Pi_RBS_inv_sigma0_well_values = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Instances{i}.Wells{j}.Pi_mumax_pmax.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_inv_sigma0_values;
                  end % for each well
               end % if not empty
           end % for each instance


             % Get synthesis predictions for a range of values of growth rate:
             Mu_vector =  0.005:0.003:0.032;
             Results_Tensor_Lib6_ALL_reduced{p,q,r}.Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,Results_Tensor_Lib6_ALL_reduced{p,q,r}.RBS_k0_sigma0_mean,...
                                                           Results_Tensor_Lib6_ALL_reduced{p,q,r}.RBS_inv_sigma0_mean,Results_Tensor_Lib6_ALL_reduced{p,q,r}.Omega_mean,Results_Tensor_Lib6_ALL_reduced{p,q,r}.Gene_cn_mean);
             Results_Tensor_Lib6_ALL_reduced{p,q,r}.Mu_slices_values = Mu_vector;
            
            % MONTE CARLO predictions for a range of values of growth rate:
            % First generate distributions for the parameters
            % To avoid values less or equal to zero, uses Tim Benham (2025). Truncated multivariate normal (https://www.mathworks.com/matlabcentral/fileexchange/34402-truncated-multivariate-normal), MATLAB Central File Exchange.
             A=-1;
             B=-1e-6;
             Results_Tensor_Lib6_ALL_reduced{p,q,r}.Gene_cn_MC_samples = rmvnrnd(Results_Tensor_Lib6_ALL_reduced{p,q,r}.Gene_cn_mean, Results_Tensor_Lib6_ALL_reduced{p,q,r}.Gene_cn_std.^2, n_samples_MC,A,B);
             Results_Tensor_Lib6_ALL_reduced{p,q,r}.Omega_MC_samples = rmvnrnd(Results_Tensor_Lib6_ALL_reduced{p,q,r}.Omega_mean, Results_Tensor_Lib6_ALL_reduced{p,q,r}.Omega_std.^2, n_samples_MC,A,B);
             Results_Tensor_Lib6_ALL_reduced{p,q,r}.RBS_k0_sigma0_MC_samples = rmvnrnd(Results_Tensor_Lib6_ALL_reduced{p,q,r}.RBS_k0_sigma0_mean, Results_Tensor_Lib6_ALL_reduced{p,q,r}.RBS_k0_sigma0_std.^2, n_samples_MC,A,B);
             MC_indexes = 1:1:n_samples_MC; 
             for j=1:n_samples_MC
                 sample_gen_cn = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Gene_cn_MC_samples(MC_indexes(j));
                 sample_Omega = Results_Tensor_Lib6_ALL_reduced{p,q,r}.Omega_MC_samples(MC_indexes(j));
                 sample_RBS_k0_sigma0 = Results_Tensor_Lib6_ALL_reduced{p,q,r}.RBS_k0_sigma0_MC_samples(MC_indexes(j));
                 Results_Tensor_Lib6_ALL_reduced{p,q,r}.MC_samples(j).Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,...
                                                              sample_RBS_k0_sigma0,RBS_inv_sigma0,sample_Omega,sample_gen_cn);
             end % for MC samples

            % Next, for each value of growth rate, we can get the probability density function of the predicted synthesis rate and get statistics from it
           MC_Pi_pred_mu = NaN*ones(n_samples_MC,length(Mu_vector)); %each column is a value of growth rate, each row is a MC sample
           for j=1:n_samples_MC
               MC_Pi_pred_mu(j,:) = Results_Tensor_Lib6_ALL_reduced{p,q,r}.MC_samples(j).Synthesis_predictions.Pi_pred_values';
           end
         
           for m=1:length(Mu_vector)
            pd = fitdist(MC_Pi_pred_mu(:,m),'kernel','Kernel','epanechnikov');
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.MC_mu_slices(m).Pi_pdf = pd;
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.MC_mu_slices(m).Mu_slice = Mu_vector(m);
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.MC_mu_slices(m).Pi_pred_mean = mean(pd);
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.MC_mu_slices(m).Pi_pred_std = std(pd);
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.MC_mu_slices(m).Pi_pred_q50 = icdf(pd,0.5);
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.MC_mu_slices(m).Pi_pred_q25 = icdf(pd,0.25);
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.MC_mu_slices(m).Pi_pred_q75 = icdf(pd,0.75);
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.MC_mu_slices(m).Pi_pred_q2p5 = icdf(pd,0.025);
            Results_Tensor_Lib6_ALL_reduced{p,q,r}.MC_mu_slices(m).Pi_pred_q97p5 = icdf(pd,0.975);
           end

        end %rbss
    end %promoters
end %plasmids

% SAVE the results
% --- Save generated results in the same folder as this script (portable) ---
file_name  = "Results_Tensor_Lib6_ALL_reduced_" + Use_mean + ".mat";
file_tensor = fullfile(SCRIPT_DIR, file_name);
save(file_tensor, "Results_Tensor_Lib6_ALL_reduced", "-mat");


