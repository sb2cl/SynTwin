function Sum_prediction_error = J5w_LogPI_Lib5_ALL_reduced(params,model_c,Use_mean,ExpData_Tensor_lib5_micro,Inherited_ParamsData,HEM,num_run,RBS_inv_sigma0,delta)
%% J5w_LogPI_Lib5_ALL_reduced
% Evaluate the all-construct Lib5 reduced-model objective.
%
% DESCRIPTION
%   Computes the weighted pseudo-Huber objective over the five Lib5
%   constructs. The fitted scalar is the J23100 promoter transcription
%   parameter, Omega.
%
%   Effective copy number and each RBS intrinsic initiation capacity are fixed
%   to the inherited parameter realization indexed by num_run.
%
% OUTPUT
%   Sum_prediction_error    Scalar objective value.
%
% The function is deterministic for fixed inputs and produces no figures or
% files.

indices_rbss = [1,2,3,4,5];
num_constructs = length(indices_rbss);
Sum_prediction_error = 0;
% NOTICE each combination (1,1,k) below will correspond to a library transcriptional unit (TU)

for k=1:length(indices_rbss)
        Omega = params;
        RBS_k0_sigma0 = Inherited_ParamsData.RBS{indices_rbss(k)}.RBS_k0_sigma0_MC_samples(num_run,1); 
        Gen_cn =  Inherited_ParamsData.Gene_cn_MC_samples(num_run);

        if strcmp(Use_mean,'Global') 
        % obtain the cost index using the mean data over the instances of the construct
            Mu_vector = ExpData_Tensor_lib5_micro{1,1,indices_rbss(k)}.Mu_mumax_pmax_global_mean;
            Pi_experimental = ExpData_Tensor_lib5_micro{1,1,indices_rbss(k)}.Pi_mumax_pmax_global_mean;
            Pi_predictions = Get_synthesis_predictions_lite(HEM,model_c, Mu_vector,RBS_k0_sigma0,RBS_inv_sigma0,Omega,Gen_cn);
            if any(isnan(Pi_predictions))
                Sum_prediction_error = Sum_prediction_error + Inf; % Handle NaN case
                 if isinf(Sum_prediction_error)
                    warning('MyWarning:InfiniteResult', 'Pi prediction resulting in Inf in TU with indexes (%d,%d,%d) for RBS_k0_sigma0=%f, Omega=%f, and Gen_cn=%f.\n', 1, 1,k,RBS_k0_sigma0,Omega,Gen_cn);
                 end
            else
                error = log10(Pi_predictions) - log10(Pi_experimental);
                rho = delta^2 * (sqrt(1 + (error./delta).^2) - 1);
                weight = sqrt(1 + log10(1 + Pi_experimental));
                Sum_prediction_error =  Sum_prediction_error +  sum(weight.*rho)/length(Mu_vector);
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
                     Pi_predictions = Get_synthesis_predictions_lite(HEM,model_c, Mu_vector,RBS_k0_sigma0,RBS_inv_sigma0,Omega,Gen_cn);
                     if any(isnan(Pi_predictions))
                         Sum_prediction_error_instance = Sum_prediction_error_instance  + Inf; % Handle NaN case
                         if isinf(Sum_prediction_error_instance)
                            warning('MyWarning:InfiniteResult', 'Pi prediction resulting in Inf in instance  %d of the TU with indexes (%d,%d,%d) for RBS_k0_sigma0=%f, Omega=%f, and Gen_cn=%f.\n', r, 1,1,k,RBS_k0_sigma0,Omega,Gen_cn);
                        end
                     else
                         error = log10(Pi_predictions) - log10(Pi_experimental);
                         rho = delta^2 * (sqrt(1 + (error./delta).^2) - 1);
                         weight = sqrt(1 + log10(1 + Pi_experimental));
                         Sum_prediction_error_instance =  Sum_prediction_error_instance +  sum(weight.*rho)/length(Mu_vector); 
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
                             Pi_predictions =  Get_synthesis_predictions_lite(HEM,model_c, Mu_vector,RBS_k0_sigma0,RBS_inv_sigma0,Omega,Gen_cn);
                             if any(isnan(Pi_predictions))
                                 Sum_prediction_error_wells= Sum_prediction_error_wells+ Inf; % Handle NaN case
                                  if isinf(Sum_prediction_error_wells)
                                    warning('MyWarning:InfiniteResult', 'Pi prediction resulting in Inf in well %d of instance  %d of the TU with indexes (%d,%d,%d) for RBS_k0_sigma0=%f, Omega=%f, and Gen_cn=%f.\n', w,r, 1, 1,k,RBS_k0_sigma0,Omega,Gen_cn);
                                 end
                             else
                                 error = log10(Pi_predictions) - log10(Pi_experimental);
                                rho = delta^2 * (sqrt(1 + (error./delta).^2) - 1);
                                weight = sqrt(1 + log10(1 + Pi_experimental));
                                Sum_prediction_error_wells =  Sum_prediction_error_wells +  sum(weight.*rho)/length(Mu_vector); 
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
Sum_prediction_error =  Sum_prediction_error/num_constructs;
end % of function

