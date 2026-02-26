function Sum_prediction_error = J4_LogPI_Lib6_ALL_reduced(params,model_c,Use_mean,ExpData_Tensor_lib30_micro,ParamsData_lib24, HEM,num_run)
%% J4_LogPI_Lib6_ALL_reduced
% SynTwin objective function (cost) for reduced-model parameter estimation
% of the L6 sublibrary (ALL scheme).
%
% CONTEXT (L6 sublibrary)
%   L6 is a 6-TU sublibrary defined by:
%     - 2 plasmid origins (Ori indices: [1 2])
%     - 3 promoters      (Promoter indices: [1 2 3])
%     - 1 shared RBS     (RBS index: [3])
%
%   In this estimation, only the shared RBS parameter is fitted. All Ori- and
%   promoter-related parameters are inherited from previously estimated Lib24
%   results (provided via ParamsData_lib24), while known constants (e.g. the
%   pSC101 baseline copy number) are provided in model_c.
%
% DESCRIPTION
%   Computes the scalar objective used by the optimizer (BADS/MEIGO).
%   For each TU in the L6 sublibrary, SynTwin digital-twin predictions of the
%   synthesis rate Pi are generated using:
%     - the Host Equivalent Model (HEM),
%     - experimentally measured growth-rate trajectories Mu(t),
%     - inherited Lib24 parameter samples (Omega and Gene copy number),
%     - the RBS parameters (k0_sigma0 estimated; inv_sigma0 fixed).
%   Predicted and experimental Pi trajectories are compared and aggregated into
%   a single library-scale cost.
%
% COST FUNCTION
%   The mismatch is computed in base-10 logarithmic space using a magnitude-
%   dependent weighting:
%
%       w_log_abs_prediction_error =
%           abs( log10(Pi_exp) .* ( log10(Pi_pred) - log10(Pi_exp) ) );
%
%   The returned scalar cost Sum_prediction_error aggregates this weighted
%   log-error across all TUs in the L6 sublibrary and across time points.
%
% INPUTS
%   params:
%       Scalar parameter to be estimated for L6 (reduced model):
%         - RBS_k0_sigma0  (intrinsic initiation capacity of the shared RBS)
%       The sensitivity-related parameter inv_sigma0 is fixed in this script.
%
%   model_c:
%       Struct with a priori known/fixed quantities and constants used by the
%       digital twin (e.g., fixed pSC101 copy number baseline).
%
%   Use_mean:
%       Aggregation level for experimental data:
%         'Global', 'Instances', or 'Wells'.
%
%   ExpData_Tensor_lib30_micro:
%       Experimental data container from which the L6 TU measurements are
%       retrieved using the indices (Ori=[1 2], Promoters=[1 2 3], RBS=3).
%       (Note: the L6 sublibrary is embedded in the full Lib30 tensor format.)
%
%   ParamsData_lib24:
%       Cell array containing Lib24 inferred parameter samples used here as
%       inherited inputs (per Ori/Promoter context), e.g.:
%         - Omega_MC_samples(num_run)
%         - Gene_cn_MC_samples(num_run)
%
%   HEM:
%       Host Equivalent Model used to compute context-aware Pi predictions.
%
%   num_run:
%       Index selecting one Monte Carlo sample from ParamsData_lib24 to be used
%       consistently across all TUs in the current cost evaluation.
%
% OUTPUTS
%   Sum_prediction_error:
%       Scalar objective value (lower is better).
%
% NOTES
%   - Called repeatedly by the optimizer; must not produce plots or interactive
%     output (warnings are allowed for debugging/diagnostics).
%   - Reduced-model L6 estimation: estimates RBS_k0_sigma0 while keeping
%     inv_sigma0 fixed (as defined in the script).
%   - Deterministic given (params, model_c, Use_mean, ExpData_Tensor_lib30_micro,
%     ParamsData_lib24, HEM, num_run).

indices_plasmids = [1,2];
indices_promoters = [1,2,3];
indices_rbss = [3];
num_constructs = length(indices_plasmids)*length(indices_promoters);
Sum_prediction_error = 0; 
% NOTICE each combination (i,j,k) below will correspond to a library transcriptional unit (TU)
for i=1:length(indices_plasmids) %Plasmids i=1 high copy (pGreen), i=2-> low copy (pSC101)
    for j=1:length(indices_promoters)
        for k=1:length(indices_rbss)
                Omega = ParamsData_lib24{i,j}.Omega_MC_samples(num_run);
                RBS_k0_sigma0 = params;
                RBS_inv_sigma0 = 0.02;
                Gen_cn =  ParamsData_lib24{i,j}.Gene_cn_MC_samples(num_run);

                if strcmp(Use_mean,'Global') 
                % obtain the cost index using the mean data over the instances of the construct
                    Mu_vector = ExpData_Tensor_lib30_micro{indices_plasmids(i),indices_promoters(j),indices_rbss(k)}.Mu_mumax_pmax_global_mean;
                    Pi_experimental = ExpData_Tensor_lib30_micro{indices_plasmids(i),indices_promoters(j),indices_rbss(k)}.Pi_mumax_pmax_global_mean;
                    Pi_predictions = Get_synthesis_predictions_lite(HEM,model_c, Mu_vector,RBS_k0_sigma0,RBS_inv_sigma0,Omega,Gen_cn);
                    if any(isnan(Pi_predictions))
                        Sum_prediction_error = Sum_prediction_error + Inf; % Handle NaN case
                         if isinf(Sum_prediction_error)
                            warning('MyWarning:InfiniteResult', 'Pi prediction resulting in Inf in TU with indexes (%d,%d,%d) for RBS_k0_sigma0=%f, Omega=%f, and Gen_cn=%f.\n', i, j,k,RBS_k0_sigma0,Omega,Gen_cn);
                         end
                    else
                        Sum_prediction_error =  Sum_prediction_error +  sum(abs( log10(Pi_predictions).*log10(Pi_experimental) - log10(Pi_experimental).^2 ))/length(Mu_vector);
                    end % If prediction is NaN
                 elseif strcmp(Use_mean,'Instances' ) 
                     % obtain the cost index using the means of each instance (experiment with set of 10 wells) of the construct   
                      instance_count = length(ExpData_Tensor_lib30_micro{indices_plasmids(i),indices_promoters(j),indices_rbss(k)}.Instances);
                       Sum_prediction_error_instance=0;
                       valid_instances = 0;
                      for r = 1:instance_count
                         if ~(isempty(ExpData_Tensor_lib30_micro{indices_plasmids(i),indices_promoters(j),indices_rbss(k)}.Instances{r}.Pi_mumax_pmax_instance_mean))
                             Mu_vector = ExpData_Tensor_lib30_micro{indices_plasmids(i),indices_promoters(j),indices_rbss(k)}.Instances{r}.Mu_mumax_pmax_instance_mean;
                             Pi_experimental = ExpData_Tensor_lib30_micro{indices_plasmids(i),indices_promoters(j),indices_rbss(k)}.Instances{r}.Pi_mumax_pmax_instance_mean;
                             Pi_predictions = Get_synthesis_predictions_lite(HEM,model_c, Mu_vector,RBS_k0_sigma0,RBS_inv_sigma0,Omega,Gen_cn);
                             if any(isnan(Pi_predictions))
                                 Sum_prediction_error_instance = Sum_prediction_error_instance + Inf; % Handle NaN case
                                 if isinf(Sum_prediction_error_instance)
                                    warning('MyWarning:InfiniteResult', 'Pi prediction resulting in Inf in instance  %d of the TU with indexes (%d,%d,%d) for RBS_k0_sigma0=%f, Omega=%f, and Gen_cn=%f.\n', r, i, j,k,RBS_k0_sigma0,Omega,Gen_cn);
                                end
                             else
                                Sum_prediction_error_instance =  Sum_prediction_error_instance +  sum(abs( log10(Pi_predictions).*log10(Pi_experimental) - log10(Pi_experimental).^2 ))/length(Mu_vector);
                             end % If prediction is NaN
                             valid_instances = valid_instances + 1;
                         end %If instance has data
                      end % for each instance of the current TU
                      if valid_instances > 0
                        Sum_prediction_error_instance = Sum_prediction_error_instance/valid_instances;
                        Sum_prediction_error = Sum_prediction_error + Sum_prediction_error_instance;
                     end
                 else 
                 % obtain the cost index using the raw data of all instances (wells x instances) of the construct
                     instance_count = length(ExpData_Tensor_lib30_micro{indices_plasmids(i),indices_promoters(j),indices_rbss(k)}.Instances);
                     Sum_prediction_error_instance=0;
                     valid_instances=0;
                      for r = 1:instance_count
                         if ~(isempty(ExpData_Tensor_lib30_micro{indices_plasmids(i),indices_promoters(j),indices_rbss(k)}.Instances{r}.Pi_mumax_pmax_instance_mean))
                             wells_count = length(ExpData_Tensor_lib30_micro{indices_plasmids(i),indices_promoters(j),indices_rbss(k)}.Instances{r}.Wells);
                             Sum_prediction_error_wells=0;
                             valid_wells=0;
                             for w = 1:wells_count
                                 if ~isempty(ExpData_Tensor_lib30_micro{indices_plasmids(i),indices_promoters(j),indices_rbss(k)}.Instances{r}.Wells{w}.Pi_mumax_pmax)
                                     % Process the data for each well
                                     Mu_vector = ExpData_Tensor_lib30_micro{indices_plasmids(i),indices_promoters(j),indices_rbss(k)}.Instances{r}.Wells{w}.Mu_mumax_pmax;
                                     Pi_experimental = ExpData_Tensor_lib30_micro{indices_plasmids(i),indices_promoters(j),indices_rbss(k)}.Instances{r}.Wells{w}.Pi_mumax_pmax;
                                     Pi_predictions = Get_synthesis_predictions_lite(HEM,model_c, Mu_vector,RBS_k0_sigma0,RBS_inv_sigma0,Omega,Gen_cn);
                                     if any(isnan(Pi_predictions))
                                         Sum_prediction_error_wells= Sum_prediction_error_wells + Inf; % Handle NaN case
                                          if isinf(Sum_prediction_error_wells)
                                            warning('MyWarning:InfiniteResult', 'Pi prediction resulting in Inf in well %d of instance  %d of the TU with indexes (%d,%d,%d) for RBS_k0_sigma0=%f, Omega=%f, and Gen_cn=%f.\n', w,r, i, j,k,RBS_k0_sigma0,Omega,Gen_cn);
                                         end
                                     else
                                        Sum_prediction_error_wells =  Sum_prediction_error_wells +  sum(abs( log10(Pi_predictions).*log10(Pi_experimental) - log10(Pi_experimental).^2 ))/length(Mu_vector);
                                     end % If prediction is NaN
                                     valid_wells= valid_wells+1;
                                 end %If wells have data
                             end %for wells 
                              if valid_wells > 0
                                Sum_prediction_error_wells =  Sum_prediction_error_wells/valid_wells; %averaged Sum_prediction_error over number of wells of the current instance of current TU
                                Sum_prediction_error_instance = Sum_prediction_error_instance + Sum_prediction_error_wells;
                                valid_instances = valid_instances + 1;
                            end
                          end %If instance has data
                      end % for each instance of the current TU
                      if valid_instances > 0
                          Sum_prediction_error_instance = Sum_prediction_error_instance/valid_instances; %averaged Sum_prediction_error over number of instances of the current TU
                          Sum_prediction_error = Sum_prediction_error + Sum_prediction_error_instance;
                      end
                end % Which data is used to get the prediction error
       end %loop for RBS
    end %loop for promoter
end %loop for plasmids
Sum_prediction_error =  Sum_prediction_error/num_constructs;
end % of function

