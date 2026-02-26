function Sum_prediction_error = J4_LogPI_Lib24_L1O_reduced(params,Construct_2_LO,model_c,Use_mean,ExpData_Tensor_lib30_micro,HEM)
%% J4_LogPI_Lib24_L1O_reduced
% SynTwin objective function (cost) for reduced-model parameter estimation 
% (Lib24, L1O).
%
% DESCRIPTION
%   Computes the scalar objective J used by the optimizer (BADS/MEIGO).
%   Given a parameter vector and an input_data structure describing one 
%   L1O case (i.e. all transcriptional units in the library but one), the 
%   function simulates SynTwin digital-twin predictions of the synthesis 
%   rate Pi and compares them against experimental estimates.
%
% COST FUNCTION
%   The mismatch is computed in base-10 logarithmic space using a
%   magnitude-dependent weighting:
%
%       w_log_abs_prediction_error =
%           abs( log10(Pi_exp_matrix) .* ...
%                ( log10(predicted_Pi_matrix) - log10(Pi_exp_matrix) ) );
%
%   This formulation emphasizes discrepancies at higher synthesis-rate
%   regimes while preserving scale invariance in log space.
%
%   The final scalar objective J aggregates these weighted log-errors
%   across all TUs and time points.
%
% INPUTS
%   param:        Parameter vector to be estimated. 
%
%      reduced-model: for RBSs, estimates k0_sigma0 (intrinsic initiation
%      capacity) while keeping inv_sigma0 fixed (as defined in the script).
%
%   input_data:   Struct containing all information required to evaluate
%                 the library-scale cost function:
%
%       - Construct_2_LO:
%             Indices defining the left-out construct (if L1O) and the
%             corresponding training set.
%
%       - model_c:
%             Struct with a priori known biopart parameters
%             (e.g., fixed plasmid copy numbers such as pSC101).
%
%       - Use_mean:
%             Aggregation level for experimental data
%             ('Global', 'Instances', or 'Wells').
%
%       - ExpData_Tensor_lib24_micro:
%             Experimental time-series data and processed observables
%             for the full combinatorial library.
%
%       - HEM:
%             Host Equivalent Model (digital-twin component) used to
%             compute context-aware synthesis-rate predictions.
%
% OUTPUTS
%   J:            Scalar objective value (lower is better).
%
% NOTES
%   - This function is called repeatedly by the optimizer and must not
%     produce plots or interactive output.
%   - The reduced model estimates the intrinsic initiation capacity
%     k0_sigma0 while keeping inv_sigma0 fixed.
%   - Deterministic given (param, input_data).

indices_plasmids_lib24 = [1,2];
indices_promoters_lib24 = [1,2,3];
indices_rbss_lib24 = [1,2,4,5];
num_constructs = length(indices_plasmids_lib24)*length(indices_promoters_lib24)*length(indices_rbss_lib24)-1;
Sum_prediction_error = 0;
% NOTICE each combination (i,j,k) below will correspond to a library transcriptional unit (TU)
for i=1:length(indices_plasmids_lib24) %Plasmids i=1 high copy (pGreen), i=2-> low copy (pSC101)
    for j=1:length(indices_promoters_lib24)
        for k=1:length(indices_rbss_lib24)
            if  ~( (i==Construct_2_LO(1,1))&&(j==Construct_2_LO(1,2))&&(k==Construct_2_LO(1,3)) ) % If construct to be considered    
                Omega = params(j);
                RBS_k0_sigma0 = params(length(indices_promoters_lib24)+k);
                RBS_inv_sigma0 = 0.02;
                Gen_cn = model_c.N_pSC101*(i-1 + mod(i,2)*params(length(indices_promoters_lib24)+length(indices_rbss_lib24)+1) );
                if strcmp(Use_mean,'Global') 
                  
                % obtain the cost index using the mean data over the instances of the construct
                    Mu_vector = ExpData_Tensor_lib30_micro{indices_plasmids_lib24(i),indices_promoters_lib24(j),indices_rbss_lib24(k)}.Mu_mumax_pmax_global_mean;
                    Pi_experimental = ExpData_Tensor_lib30_micro{indices_plasmids_lib24(i),indices_promoters_lib24(j),indices_rbss_lib24(k)}.Pi_mumax_pmax_global_mean;
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
                      instance_count = length(ExpData_Tensor_lib30_micro{indices_plasmids_lib24(i),indices_promoters_lib24(j),indices_rbss_lib24(k)}.Instances);
                      Sum_prediction_error_instance=0;
                      valid_instances = 0;
                      for r = 1:instance_count
                         if ~(isempty(ExpData_Tensor_lib30_micro{indices_plasmids_lib24(i),indices_promoters_lib24(j),indices_rbss_lib24(k)}.Instances{r}.Pi_mumax_pmax_instance_mean))
                             Mu_vector = ExpData_Tensor_lib30_micro{indices_plasmids_lib24(i),indices_promoters_lib24(j),indices_rbss_lib24(k)}.Instances{r}.Mu_mumax_pmax_instance_mean;
                             Pi_experimental = ExpData_Tensor_lib30_micro{indices_plasmids_lib24(i),indices_promoters_lib24(j),indices_rbss_lib24(k)}.Instances{r}.Pi_mumax_pmax_instance_mean;
                             Pi_predictions = Get_synthesis_predictions_lite(HEM,model_c, Mu_vector,RBS_k0_sigma0,RBS_inv_sigma0,Omega,Gen_cn);
                             if any(isnan(Pi_predictions))
                                 Sum_prediction_error_instance = Sum_prediction_error_instance + Inf; % Handle NaN case
                                 if isinf(Sum_prediction_error_instance)
                                    warning('MyWarning:InfiniteResult', 'Pi prediction resulting in Inf in instance  %d of the TU with indexes (%d,%d,%d) for RBS_k0_sigma0=%f, Omega=%f, and Gen_cn=%f.\n', r, i, j,k,RBS_k0_sigma0,Omega,Gen_cn);
                                end
                             else
                                Sum_prediction_error_instance =  Sum_prediction_error_instance  +  sum(abs( log10(Pi_predictions).*log10(Pi_experimental) - log10(Pi_experimental).^2 ))/length(Mu_vector);
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
                     instance_count = length(ExpData_Tensor_lib30_micro{indices_plasmids_lib24(i),indices_promoters_lib24(j),indices_rbss_lib24(k)}.Instances);
                     Sum_prediction_error_instance=0;
                     valid_instances=0;
                      for r = 1:instance_count
                         if ~(isempty(ExpData_Tensor_lib30_micro{indices_plasmids_lib24(i),indices_promoters_lib24(j),indices_rbss_lib24(k)}.Instances{r}.Pi_mumax_pmax_instance_mean))
                             wells_count = length(ExpData_Tensor_lib30_micro{indices_plasmids_lib24(i),indices_promoters_lib24(j),indices_rbss_lib24(k)}.Instances{r}.Wells);
                             Sum_prediction_error_wells=0;
                             valid_wells=0;
                             for w = 1:wells_count
                                 if ~isempty(ExpData_Tensor_lib30_micro{indices_plasmids_lib24(i),indices_promoters_lib24(j),indices_rbss_lib24(k)}.Instances{r}.Wells{w}.Pi_mumax_pmax)
                                     % Process the data for each well
                                     Mu_vector = ExpData_Tensor_lib30_micro{indices_plasmids_lib24(i),indices_promoters_lib24(j),indices_rbss_lib24(k)}.Instances{r}.Wells{w}.Mu_mumax_pmax;
                                     Pi_experimental = ExpData_Tensor_lib30_micro{indices_plasmids_lib24(i),indices_promoters_lib24(j),indices_rbss_lib24(k)}.Instances{r}.Wells{w}.Pi_mumax_pmax;
                                     Pi_predictions = Get_synthesis_predictions_lite(HEM,model_c, Mu_vector,RBS_k0_sigma0,RBS_inv_sigma0,Omega,Gen_cn);
                                     if any(isnan(Pi_predictions))
                                         Sum_prediction_error_wells= Sum_prediction_error_wells + Inf;  % Handle NaN case
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
            end % If construct to be considered
       end %loop for RBS
    end %loop for promoter
end %loop for plasmids
Sum_prediction_error =  Sum_prediction_error/num_constructs;
end % of function

