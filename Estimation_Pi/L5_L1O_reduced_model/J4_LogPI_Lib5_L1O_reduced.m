function Sum_prediction_error = J4_LogPI_Lib5_L1O_reduced(params,Construct_2_LO,model_c,Use_mean,ExpData_Tensor_lib5_micro,Inherited_ParamsData, HEM,num_run)
%% J4_LogPI_Lib5_L1O_reduced
% SynTwin objective function (cost) for reduced-model parameter estimation
% of the Lib5 sublibrary under L1O (leave-one-out) cross-validation.
%
% CONTEXT (Lib5 sublibrary)
%   Lib5 is a 5-TU sublibrary defined by:
%     - 1 plasmid origin  (pGreen; Ori index fixed to 1)
%     - 1 promoter        (J23100; parameter to be estimated)
%     - 5 RBSs            (RBS indices: [1 2 3 4 5])
%
%   In this L1O scheme, one RBS-defined TU is left out at each iteration.
%   The promoter strength parameter Omega is estimated using the remaining
%   4 TUs, and the left-out TU is excluded from the cost computation.
%
%   As in the ALL version:
%     - Gene copy number (Gene_cn) is inherited.
%     - RBS k0_sigma0 parameters are inherited per RBS.
%       * RBS B0034 (index 3) is inherited from Lib6 inference.
%       * The remaining RBS parameters are inherited from Lib24.
%     - inv_sigma0 is fixed inside the script.
%
% DESCRIPTION
%   Computes the scalar objective used by the optimizer (BADS/MEIGO)
%   under an L1O cross-validation scheme.
%
%   For each TU in Lib5 (Ori=1, Promoter=1, RBS=k), except the left-out
%   construct defined by Construct_2_LO, SynTwin digital-twin predictions
%   of the synthesis rate Pi are generated using:
%     - the Host Equivalent Model (HEM),
%     - experimentally measured growth-rate trajectories Mu(t),
%     - inherited parameter samples (Gene_cn and RBS k0_sigma0),
%     - the promoter parameter Omega (estimated here),
%     - a fixed inv_sigma0 value.
%
%   Predicted and experimental Pi trajectories are compared and aggregated
%   into a library-scale cost using only the training TUs.
%
% COST FUNCTION
%   The mismatch is computed in base-10 logarithmic space using a
%   magnitude-dependent weighting:
%
%       w_log_abs_prediction_error =
%           abs( log10(Pi_pred) .* log10(Pi_exp) - log10(Pi_exp).^2 );
%
%   The returned scalar cost Sum_prediction_error aggregates this weighted
%   log-error across the 4 training TUs and across time points.
%
% INPUTS
%   params:
%       Scalar parameter to be estimated:
%         - Omega (promoter strength; J23100)
%
%   Construct_2_LO:
%       Scalar index (1–5) identifying the RBS-defined TU to be left out
%       in the current L1O iteration.
%
%   model_c:
%       Struct with a priori known/fixed quantities and constants used by
%       the digital twin.
%
%   Use_mean:
%       Aggregation level for experimental data:
%         'Global', 'Instances', or 'Wells'.
%
%   ExpData_Tensor_lib5_micro:
%       Experimental data container for Lib5. Data are retrieved at indices
%       {Ori=1, Promoter=1, RBS=k}.
%
%   Inherited_ParamsData:
%       Struct containing inherited parameter samples used as fixed inputs
%       for the current cost evaluation, e.g.:
%         - Inherited_ParamsData.Gene_cn_MC_samples(num_run)
%         - Inherited_ParamsData.RBS{k}.RBS_k0_sigma0_MC_samples(num_run,1)
%
%   HEM:
%       Host Equivalent Model used to compute context-aware Pi predictions.
%
%   num_run:
%       Index selecting one Monte Carlo sample from the inherited parameter
%       sets to be used consistently across all training TUs.
%
% OUTPUTS
%   Sum_prediction_error:
%       Scalar objective value (lower is better).
%
% NOTES
%   - Called repeatedly by the optimizer; must not produce plots or
%     interactive output (warnings allowed for diagnostics).
%   - Reduced-model Lib5 L1O estimation: estimates Omega while keeping
%     RBS parameters inherited and inv_sigma0 fixed.
%   - Deterministic given (params, Construct_2_LO, model_c, Use_mean,
%     ExpData_Tensor_lib5_micro, Inherited_ParamsData, HEM, num_run).

indices_rbss = [1,2,3,4,5];
num_constructs = length(indices_rbss);
Sum_prediction_error = 0;
% NOTICE each combination (1,1,k) below will correspond to a library transcriptional unit (TU)

for k=1:length(indices_rbss)
    if  ~( k==Construct_2_LO )% If construct to be considered    (varies as 1,2,4,5 , the considered RBSs)
        Omega = params;
        RBS_k0_sigma0 = Inherited_ParamsData.RBS{indices_rbss(k)}.RBS_k0_sigma0_MC_samples(num_run,1); 
        RBS_1_sigma0 = 0.02;
        Gen_cn =  Inherited_ParamsData.Gene_cn_MC_samples(num_run);

        if strcmp(Use_mean,'Global') 
        % obtain the cost index using the mean data over the instances of the construct
            Mu_vector = ExpData_Tensor_lib5_micro{1,1,indices_rbss(k)}.Mu_mumax_pmax_global_mean;
            Pi_experimental = ExpData_Tensor_lib5_micro{1,1,indices_rbss(k)}.Pi_mumax_pmax_global_mean;
            Pi_predictions = Get_synthesis_predictions_lite(HEM,model_c, Mu_vector,RBS_k0_sigma0,RBS_1_sigma0,Omega,Gen_cn);
            if any(isnan(Pi_predictions))
                Sum_prediction_error = Sum_prediction_error + Inf; % Handle NaN case
                 if isinf(Sum_prediction_error)
                    warning('MyWarning:InfiniteResult', 'Pi prediction resulting in Inf in TU with indexes (%d,%d,%d) for RBS_k0_sigma0=%f, Omega=%f, and Gen_cn=%f.\n', 1, 1,k,RBS_k0_sigma0,Omega,Gen_cn);
                 end
            else
                Sum_prediction_error =  Sum_prediction_error +  sum(abs( log10(Pi_predictions).*log10(Pi_experimental) - log10(Pi_experimental).^2 ))/length(Mu_vector);
            end % If prediction is NaN
         elseif strcmp(Use_mean,'Instances' ) 
             % obtain the cost index using the means of each instance (experiment with set of 10 wells) of the construct   
              instance_count = length(ExpData_Tensor_lib5_micro{1,1,indices_rbss(k)}.Instances);
               Sum_prediction_error_instance=0;
               valid_instances = 0;
              for r = 1:instance_count
                 if ~(isempty(ExpData_Tensor_lib5_micro{1,1,indices_rbss(k)}.Instances{r}.Pi_mumax_pmax_instance_mean))
                     Mu_vector = ExpData_Tensor_lib5_micro{1,1,indices_rbss(k)}.Instances{r}.Mu_mumax_pmax_instance_mean;
                     Pi_experimental = ExpData_Tensor_lib5_micro{1,1,indices_rbss(k)}.Instances{r}.Pi_mumax_pmax_instance_mean;
                     Pi_predictions = Get_synthesis_predictions_lite(HEM,model_c, Mu_vector,RBS_k0_sigma0,RBS_1_sigma0,Omega,Gen_cn);
                     if any(isnan(Pi_predictions))
                         Sum_prediction_error_instance = Sum_prediction_error_instance  + Inf; % Handle NaN case
                         if isinf(Sum_prediction_error_instance)
                            warning('MyWarning:InfiniteResult', 'Pi prediction resulting in Inf in instance  %d of the TU with indexes (%d,%d,%d) for RBS_k0_sigma0=%f, Omega=%f, and Gen_cn=%f.\n', r, 1,1,k,RBS_k0_sigma0,Omega,Gen_cn);
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
             instance_count = length(ExpData_Tensor_lib5_micro{1,1,indices_rbss(k)}.Instances);
             Sum_prediction_error_instance=0;
             valid_instances=0;
              for r = 1:instance_count
                 if ~(isempty(ExpData_Tensor_lib5_micro{1,1,indices_rbss(k)}.Instances{r}.Pi_mumax_pmax_instance_mean))
                     wells_count = length(ExpData_Tensor_lib5_micro{1,1,indices_rbss(k)}.Instances{r}.Wells);
                      Sum_prediction_error_wells=0;
                     valid_wells=0;
                     for w = 1:wells_count
                         if ~isempty(ExpData_Tensor_lib5_micro{1,1,indices_rbss(k)}.Instances{r}.Wells{w}.Pi_mumax_pmax)
                             % Process the data for each well
                             Mu_vector = ExpData_Tensor_lib5_micro{1,1,indices_rbss(k)}.Instances{r}.Wells{w}.Mu_mumax_pmax;
                             Pi_experimental = ExpData_Tensor_lib5_micro{1,1,indices_rbss(k)}.Instances{r}.Wells{w}.Pi_mumax_pmax;
                             Pi_predictions =  Get_synthesis_predictions_lite(HEM,model_c, Mu_vector,RBS_k0_sigma0,RBS_1_sigma0,Omega,Gen_cn);
                             if any(isnan(Pi_predictions))
                                 Sum_prediction_error_wells= Sum_prediction_error_wells+ Inf; % Handle NaN case
                                  if isinf(Sum_prediction_error_wells)
                                    warning('MyWarning:InfiniteResult', 'Pi prediction resulting in Inf in well %d of instance  %d of the TU with indexes (%d,%d,%d) for RBS_k0_sigma0=%f, Omega=%f, and Gen_cn=%f.\n', w,r, 1, 1,k,RBS_k0_sigma0,Omega,Gen_cn);
                                 end
                             else
                                Sum_prediction_error_wells =  Sum_prediction_error_wells+  sum(abs( log10(Pi_predictions).*log10(Pi_experimental) - log10(Pi_experimental).^2 ))/length(Mu_vector);
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
Sum_prediction_error =  Sum_prediction_error/num_constructs;
end % of function

