%% Generate_Results_Lib30_L1O_full_model
% SynTwin workflow script: post-processing and compilation of Lib30 L1O full-model results.
%
% DESCRIPTION
%   Loads the per-construct optimization outputs produced by
%   Estimate_Lib30_L1O_full_model and compiles them into a unified Results tensor
%   (Global / Instances / Wells), ready for plotting and downstream analysis.
%
% INPUTS
%   None (loads results from ./Estimated_results/ and uses local configuration).
%
%   Configuration is set inside the script.
%   The user must set the desired <Use_mean> option. Options are: 
%       - 'Global' (global mean for each construct),
%       - 'Instances' (mean of each experiment for each construct),
%       - 'Wells' (use data of all individual culture wells)
%
% OUTPUTS (saved to disk)
%   Results_Tensor_Lib30_L1O_full_model_<Use_mean>.mat
%     Example: Results_Tensor_Lib30_L1O_full_model_Instances.mat
%
% DEPENDENCIES
%   - SynTwin initialization recommended:
%       ROOT = init_SynTwin(...);
%   - Requires that ./Estimated_results/ contains the expected optimization files.
%
% USAGE
%   Generate_Results_Lib30_L1O_full_model
%
% NOTES
%   - Use_mean controls the aggregation level: 'Global', 'Instances', or 'Wells'.
%   - The output file is typically loaded by Show_Results_* scripts.

% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin('experimental',true);
% Ensure this script folder is on the path (for local functions in this folder)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);

% --- Load data using portable absolute paths ---
load(SynTwin_path('Generate_HEM','HEM_Surrogate','HEM_Surrogate.mat'));          % loads data of the Host Equivalent Model (HEM) 
load(SynTwin_path('Experimental_Data','ExpData_Tensor_lib30_micro.mat'));        % loads data of Lib30: ExpData_Tensor_lib30_micro

% Options are: 
% - 'Global' (global mean for each construct),
% - 'Instances' (mean of each experiment for each construct),
% - 'No' (use data of all individual culture wells)

Use_mean = 'Instances';  
%Use_mean = 'Global';  
%Use_mean = 'Wells'; 

indices_plasmids_lib30 = [1,2];
num_plasmids_Lib30 = length(indices_plasmids_lib30);
indices_promoters_lib30 = [1,2,3];
num_promoters_Lib30 = length(indices_promoters_lib30);
indices_rbss_lib30 = [1,2,3,4,5];
num_rbss_Lib30 = length(indices_rbss_lib30);
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

% Loads the estimated parameters and generates the structure Estimated_parameters
Estimated_parameters.TU = {};
Matrix_tempo_All =[];
for p=1:length(indices_plasmids_lib30)
    for q=1:length(indices_promoters_lib30)
        for r=1:length(indices_rbss_lib30)
            Construct_2_LO = [indices_plasmids_lib30(p),indices_promoters_lib30(q),indices_rbss_lib30(r)];
            % --- Load estimation results from the local Estimated_results folder (portable) ---
            results_dir = fullfile(SCRIPT_DIR,'Estimated_results');
            if ~exist(results_dir,'dir')
                 error('Estimated_results folder not found: %s', results_dir);
            end
            file_name = "Results_BADS_Lib30_L1O_complete_" + ...
                            num2str(indices_plasmids_lib30(p)) + ...
                            num2str(indices_promoters_lib30(q)) + ...
                            num2str(indices_rbss_lib30(r)) + "_" + ...
                            Use_mean + ".mat";
            file_tensor = fullfile(results_dir, file_name);
            load(file_tensor, "Results_BADS_Lib30_L1O_complete_v2", "-mat");
            num_runs = length(Results_BADS_Lib30_L1O_complete_v2);
            tempo = [];
            for i=1:num_runs
               tempo=[tempo;Results_BADS_Lib30_L1O_complete_v2{i}.results(:,1:14)];
               Matrix_tempo_All=[Matrix_tempo_All;Results_BADS_Lib30_L1O_complete_v2{i}.results(:,1:14)];
            end
            Estimated_parameters.TU{p,q,r}.raw = tempo;
            Estimated_parameters.TU{p,q,r}.mean = mean(tempo,1);
            Estimated_parameters.TU{p,q,r}.std = std(tempo, 0, 1);
        end %rbss
    end %promoters
end %plasmids
Estimated_parameters.ALL_raw = Matrix_tempo_All;
Estimated_parameters.ALL_mean = mean(Matrix_tempo_All,1);
Estimated_parameters.ALL_std = std(Matrix_tempo_All,0,1);

% Adds the experimental data to the structure
% Results_Tensor_Lib30_L1O_complete
Results_Tensor_Lib30_L1O_complete = ExpData_Tensor_lib30_micro;

% Adds the estimated parameters and their statistics to the structure Results_Tensor_Lib30_L1O_complete
for p=1:length(indices_plasmids_lib30)
    for q=1:length(indices_promoters_lib30)
        for r=1:length(indices_rbss_lib30)
            Construct_2_LO = [indices_plasmids_lib30(p),indices_promoters_lib30(q),indices_rbss_lib30(r)];
            if strcmp(Use_mean,'Wells' )   
                 Results_Tensor_Lib30_L1O_complete{p,q,r}.Use_mean = 'Wells'; 
            elseif strcmp(Use_mean,'Instances')
                 Results_Tensor_Lib30_L1O_complete{p,q,r}.Use_mean = 'Instances'; 
            elseif strcmp(Use_mean,'Global')
                 Results_Tensor_Lib30_L1O_complete{p,q,r}.Use_mean ='Global'; 
            else
                  error('Invalid option for Use_mean. Choose either "Wells", "Instances", or "Global".');
            end
            %Local estimations are the estimated parameters for all TUs but
            %in the (p,q,r)th-LOOCV iteration
            Results_Tensor_Lib30_L1O_complete{p,q,r}.Parameters_local_raw = Estimated_parameters.TU{p,q,r}.raw;
            Results_Tensor_Lib30_L1O_complete{p,q,r}.Parameters_local_mean = Estimated_parameters.TU{p,q,r}.mean;
            Results_Tensor_Lib30_L1O_complete{p,q,r}.Parameters_local_std = Estimated_parameters.TU{p,q,r}.std;
           %Estimated parameters obtained using the resuls of all the LOOCV iterations:
           %Plasmids p=1 high copy (pGreen), p=2-> low copy (pSC101)
            Results_Tensor_Lib30_L1O_complete{p,q,r}.Gene_cn_raw = model_c.N_pSC101*(p-1 + mod(p,2)*Estimated_parameters.ALL_raw(:,num_promoters_Lib30+2*num_rbss_Lib30+1)); 
            Results_Tensor_Lib30_L1O_complete{p,q,r}.Gene_cn_mean = model_c.N_pSC101*(p-1 + mod(p,2)*Estimated_parameters.ALL_mean(:,num_promoters_Lib30+2*num_rbss_Lib30+1)); 
            Results_Tensor_Lib30_L1O_complete{p,q,r}.Gene_cn_std =  model_c.N_pSC101*mod(p,2)*Estimated_parameters.ALL_std(:,num_promoters_Lib30+2*num_rbss_Lib30+1); 
            Results_Tensor_Lib30_L1O_complete{p,q,r}.Omega_raw = Estimated_parameters.ALL_raw(:,q); 
            Results_Tensor_Lib30_L1O_complete{p,q,r}.Omega_mean =  Estimated_parameters.ALL_mean(:,q); 
            Results_Tensor_Lib30_L1O_complete{p,q,r}.Omega_std =  Estimated_parameters.ALL_std(:,q); 
            Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_k0_sigma0_raw = Estimated_parameters.ALL_raw(:,num_promoters_Lib30+r); 
            Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_k0_sigma0_mean = Estimated_parameters.ALL_mean(:,num_promoters_Lib30+r); 
            Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_k0_sigma0_std = Estimated_parameters.ALL_std(:,num_promoters_Lib30+r);
            Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_inv_sigma0_raw = Estimated_parameters.ALL_raw(:,num_promoters_Lib30+num_rbss_Lib30+r); 
            Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_inv_sigma0_mean = Estimated_parameters.ALL_mean(:,num_promoters_Lib30+num_rbss_Lib30+r); 
            Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_inv_sigma0_std = Estimated_parameters.ALL_std(:,num_promoters_Lib30+num_rbss_Lib30+r);

            %Get the sensitivities of the synthesis rate w.r.t. the estimated parameters at the experimental values of growth rate:
              % GLOBAL 
              Mu_vector = Results_Tensor_Lib30_L1O_complete{p,q,r}.Mu_mumax_pmax_global_mean;
              Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_k0_sigma0_mean,...
                                           Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_inv_sigma0_mean,Results_Tensor_Lib30_L1O_complete{p,q,r}.Omega_mean,Results_Tensor_Lib30_L1O_complete{p,q,r}.Gene_cn_mean);
              Results_Tensor_Lib30_L1O_complete{p,q,r}.S_Pi_NA_global_values = Results_Tensor_Lib30_L1O_complete{p,q,r}.Pi_mumax_pmax_global_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_NA_values;
              Results_Tensor_Lib30_L1O_complete{p,q,r}.S_Pi_Omega_global_values = Results_Tensor_Lib30_L1O_complete{p,q,r}.Pi_mumax_pmax_global_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_Omega_values;
              Results_Tensor_Lib30_L1O_complete{p,q,r}.S_Pi_RBS_k0_sigma0_global_values = Results_Tensor_Lib30_L1O_complete{p,q,r}.Pi_mumax_pmax_global_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_k0_sigma0_values;
              Results_Tensor_Lib30_L1O_complete{p,q,r}.S_Pi_RBS_inv_sigma0_global_values = Results_Tensor_Lib30_L1O_complete{p,q,r}.Pi_mumax_pmax_global_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_inv_sigma0_values;
             % INSTANCES
             instance_count = length(Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances);
             for i = 1:instance_count
               if ~(isempty(Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Mu_mumax_pmax_instance_mean)) %some experiments can be empty if all wells of a construct are outliers
                 Mu_vector = Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Mu_mumax_pmax_instance_mean;
                 Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_k0_sigma0_mean,...
                                           Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_inv_sigma0_mean,Results_Tensor_Lib30_L1O_complete{p,q,r}.Omega_mean,Results_Tensor_Lib30_L1O_complete{p,q,r}.Gene_cn_mean);
                 Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.S_Pi_NA_instance_values = Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Pi_mumax_pmax_instance_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_NA_values;
                 Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.S_Pi_Omega_instance_values = Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Pi_mumax_pmax_instance_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_Omega_values;
                 Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.S_Pi_RBS_k0_sigma0_instance_values = Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Pi_mumax_pmax_instance_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_k0_sigma0_values;
                 Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.S_Pi_RBS_inv_sigma0_instance_values = Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Pi_mumax_pmax_instance_mean.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_inv_sigma0_values;
               end % if not empty
            end % for each instance
            % WELLS
            for i = 1:instance_count
               if ~(isempty(Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Mu_mumax_pmax_instance_mean)) %some experiments can be empty if all wells of a construct are outliers
                  no_wells = length(Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Wells);
                  for j = 1:no_wells
                      Mu_vector = Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Wells{j}.Mu_mumax_pmax;
                      Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_k0_sigma0_mean,...
                                           Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_inv_sigma0_mean,Results_Tensor_Lib30_L1O_complete{p,q,r}.Omega_mean,Results_Tensor_Lib30_L1O_complete{p,q,r}.Gene_cn_mean);
                      Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Wells{j}.S_Pi_NA_well_values = Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Wells{j}.Pi_mumax_pmax.*Synthesis_predictions.Relative_sensitivities.S_Pi_NA_values;
                      Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Wells{j}.S_Pi_Omega_well_values = Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Wells{j}.Pi_mumax_pmax.*Synthesis_predictions.Relative_sensitivities.S_Pi_Omega_values;
                      Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Wells{j}.S_Pi_RBS_k0_sigma0_well_values = Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Wells{j}.Pi_mumax_pmax.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_k0_sigma0_values; 
                      Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Wells{j}.S_Pi_RBS_inv_sigma0_well_values = Results_Tensor_Lib30_L1O_complete{p,q,r}.Instances{i}.Wells{j}.Pi_mumax_pmax.*Synthesis_predictions.Relative_sensitivities.S_Pi_RBS_inv_sigma0_values;
                  end % for each well
               end % if not empty
           end % for each instance

            % Gets synthesis predictions for a range of values of growth rate:
            Mu_vector =  0.005:0.003:0.032;
            Results_Tensor_Lib30_L1O_complete{p,q,r}.Mu_slices_values = Mu_vector;
            Results_Tensor_Lib30_L1O_complete{p,q,r}.Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_k0_sigma0_mean,...
                                           Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_inv_sigma0_mean,Results_Tensor_Lib30_L1O_complete{p,q,r}.Omega_mean,Results_Tensor_Lib30_L1O_complete{p,q,r}.Gene_cn_mean);

            % MONTE CARLO predictions for a range of values of growth rate:
            % First generate distributions for the parameters
            % To avoid values less or equal to zero, uses Tim Benham (2025). Truncated multivariate normal (https://www.mathworks.com/matlabcentral/fileexchange/34402-truncated-multivariate-normal), MATLAB Central File Exchange.
             A=-1;
             B=-1e-6;
             Results_Tensor_Lib30_L1O_complete{p,q,r}.Gene_cn_MC_samples = rmvnrnd(Results_Tensor_Lib30_L1O_complete{p,q,r}.Gene_cn_mean, Results_Tensor_Lib30_L1O_complete{p,q,r}.Gene_cn_std.^2, n_samples_MC,A,B);
             Results_Tensor_Lib30_L1O_complete{p,q,r}.Omega_MC_samples = rmvnrnd(Results_Tensor_Lib30_L1O_complete{p,q,r}.Omega_mean, Results_Tensor_Lib30_L1O_complete{p,q,r}.Omega_std.^2, n_samples_MC,A,B);
             Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_k0_sigma0_MC_samples = rmvnrnd(Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_k0_sigma0_mean, Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_k0_sigma0_std.^2, n_samples_MC,A,B);
             Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_inv_sigma0_MC_samples = rmvnrnd(Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_inv_sigma0_mean, Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_inv_sigma0_std.^2, n_samples_MC,A,B);
             MC_indexes = 1:1:n_samples_MC; 
             for j=1:n_samples_MC
                 sample_gen_cn = Results_Tensor_Lib30_L1O_complete{p,q,r}.Gene_cn_MC_samples(MC_indexes(j));
                 sample_Omega = Results_Tensor_Lib30_L1O_complete{p,q,r}.Omega_MC_samples(MC_indexes(j));
                 sample_RBS_k0_sigma0 = Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_k0_sigma0_MC_samples(MC_indexes(j));
                 sample_RBS_inv_sigma0 = Results_Tensor_Lib30_L1O_complete{p,q,r}.RBS_inv_sigma0_MC_samples(MC_indexes(j));
                 Results_Tensor_Lib30_L1O_complete{p,q,r}.MC_samples(j).Synthesis_predictions = Get_synthesis_predictions(HEM,model_c,Mu_vector,...
                                                              sample_RBS_k0_sigma0,sample_RBS_inv_sigma0,sample_Omega,sample_gen_cn);
             end % for MC samples
            % Next, for each value of growth rate, we can get the probability density function of the predicted synthesis rate and get statistics from it
            MC_Pi_pred_mu = NaN*ones(n_samples_MC,length(Mu_vector)); %each column is a value of growth rate, each row is a MC sample
               for j=1:n_samples_MC
                   MC_Pi_pred_mu(j,:) = Results_Tensor_Lib30_L1O_complete{p,q,r}.MC_samples(j).Synthesis_predictions.Pi_pred_values';
               end
               for m=1:length(Mu_vector)
                    pd = fitdist(MC_Pi_pred_mu(:,m),'kernel','Kernel','epanechnikov');
                    Results_Tensor_Lib30_L1O_complete{p,q,r}.MC_mu_slices(m).Pi_pdf = pd;
                    Results_Tensor_Lib30_L1O_complete{p,q,r}.MC_mu_slices(m).Mu_slice = Mu_vector(m);
                    Results_Tensor_Lib30_L1O_complete{p,q,r}.MC_mu_slices(m).Pi_pred_mean = mean(pd);
                    Results_Tensor_Lib30_L1O_complete{p,q,r}.MC_mu_slices(m).Pi_pred_std = std(pd);
                    Results_Tensor_Lib30_L1O_complete{p,q,r}.MC_mu_slices(m).Pi_pred_q50 = icdf(pd,0.5);
                    Results_Tensor_Lib30_L1O_complete{p,q,r}.MC_mu_slices(m).Pi_pred_q25 = icdf(pd,0.25);
                    Results_Tensor_Lib30_L1O_complete{p,q,r}.MC_mu_slices(m).Pi_pred_q75 = icdf(pd,0.75);
                    Results_Tensor_Lib30_L1O_complete{p,q,r}.MC_mu_slices(m).Pi_pred_q2p5 = icdf(pd,0.025);
                    Results_Tensor_Lib30_L1O_complete{p,q,r}.MC_mu_slices(m).Pi_pred_q97p5 = icdf(pd,0.975);
               end
          
        end %rbss
    end %promoters
end %plasmids

% SAVE the results
% --- Save generated results in the same folder as this script (portable) ---
file_name  = "Results_Tensor_Lib30_L1O_full_model_" + Use_mean + ".mat";
file_tensor = fullfile(SCRIPT_DIR, file_name);
save(file_tensor, "Results_Tensor_Lib30_L1O_complete", "-mat");
