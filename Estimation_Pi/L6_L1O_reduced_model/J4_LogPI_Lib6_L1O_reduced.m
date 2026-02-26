function Sum_prediction_error = J4_LogPI_Lib6_L1O_reduced(params,Construct_2_LO,model_c,Use_mean,ExpData_Tensor_lib30_micro,ParamsData_lib24, HEM,num_run)
%% J4_LogPI_Lib6_L1O_reduced
% SynTwin objective function (cost) for reduced-model parameter estimation
% of the L6 sublibrary under the Leave-One-Out (L1O) scheme.
%
% CONTEXT (L6 sublibrary)
%   L6 is a 6-TU sublibrary defined by:
%     - 2 plasmid origins (Ori indices: [1 2])
%     - 3 promoters       (Promoter indices: [1 2 3])
%     - 1 shared RBS      (RBS index: [3])
%
%   In this estimation:
%     - Only the shared RBS intrinsic initiation capacity (k0_sigma0) is estimated.
%     - Ori- and promoter-related parameters (Gene copy number and Omega)
%       are inherited from previously estimated Lib24 results.
%     - Known constants (e.g. pSC101 baseline copy number) are provided via model_c.
%
% DESCRIPTION
%   Computes the scalar objective used by the optimizer (BADS/MEIGO)
%   for one Leave-One-Out (L1O) configuration of L6.
%
%   For each TU in the L6 sublibrary except the construct defined in
%   Construct_2_LO, the function:
%     - retrieves experimentally measured growth-rate trajectories Mu(t),
%     - computes SynTwin digital-twin predictions of the synthesis rate Pi,
%     - compares predicted and experimental Pi values,
%     - aggregates the mismatch into a single library-scale cost.
%
%   The left-out construct is excluded from the cost evaluation and is
%   reserved for validation.
%
% COST FUNCTION
%   The mismatch is computed in base-10 logarithmic space using a magnitude-
%   dependent weighting:
%
%       w_log_abs_prediction_error =
%           abs( log10(Pi_exp) .* ( log10(Pi_pred) - log10(Pi_exp) ) );
%
%   The returned scalar cost Sum_prediction_error aggregates this weighted
%   log-error across:
%       - all training TUs (excluding the left-out one),
%       - all considered time points.
%
% INPUTS
%   params:
%       Scalar parameter to be estimated (reduced model):
%         - RBS_k0_sigma0  (intrinsic initiation capacity of the shared RBS)
%       The sensitivity-related parameter inv_sigma0 is fixed internally.
%
%   Construct_2_LO:
%       1×3 vector defining the TU to be excluded:
%         [plasmid_index, promoter_index, rbs_index]
%
%   model_c:
%       Struct with known/fixed quantities and digital-twin constants.
%
%   Use_mean:
%       Aggregation level for experimental data:
%         'Global', 'Instances', or 'Wells'.
%
%   ExpData_Tensor_lib30_micro:
%       Experimental data container from which L6 TU measurements are
%       retrieved using the L6 indices.
%       (Note: L6 is embedded within the Lib30 tensor structure.)
%
%   ParamsData_lib24:
%       Cell array containing inherited Lib24 parameter samples:
%         - Omega_MC_samples(num_run)
%         - Gene_cn_MC_samples(num_run)
%
%   HEM:
%       Host Equivalent Model used to compute context-aware Pi predictions.
%
%   num_run:
%       Monte Carlo index selecting one inherited parameter sample from
%       Lib24 (used consistently across all training TUs).
%
% OUTPUTS
%   Sum_prediction_error:
%       Scalar objective value (lower is better).
%
% NOTES
%   - Called repeatedly by the optimizer; must not produce plots or
%     interactive output.
%   - The left-out TU is excluded from cost computation.
%   - Reduced-model L6 estimation: estimates RBS_k0_sigma0 while keeping
%     inv_sigma0 fixed.
%   - Deterministic given (params, Construct_2_LO, model_c, Use_mean,
%     ExpData_Tensor_lib30_micro, ParamsData_lib24, HEM, num_run).

indices_plasmids = [1,2];
indices_promoters = [1,2,3];
indices_rbss = [3];
num_constructs = length(indices_plasmids)*length(indices_promoters)-1;
Sum_prediction_error = 0;
% NOTICE each combination (i,j,k) below will correspond to a library transcriptional unit (TU)
for i=1:length(indices_plasmids) %Plasmids i=1 high copy (pGreen), i=2-> low copy (pSC101)
    for j=1:length(indices_promoters)
        for k=1:length(indices_rbss)
            if  ~( (i==Construct_2_LO(1,1))&&(j==Construct_2_LO(1,2))&&(k==Construct_2_LO(1,3)) )% If construct to be considered    
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
                                        Sum_prediction_error_wells= Sum_prediction_error_wells +  sum(abs( log10(Pi_predictions).*log10(Pi_experimental) - log10(Pi_experimental).^2 ))/length(Mu_vector);
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
            end % If construct to be considered
       end %loop for RBS
    end %loop for promoter
end %loop for plasmids
Sum_prediction_error =  Sum_prediction_error/num_constructs;
end % of function

