function Sum_prediction_error = J5w_LogPI_Lib24_ALL_reduced(params,model_c,Use_mean,ExpData_Tensor_lib30_micro,HEM,RBS_inv_sigma0,delta)
%% J5w_LogPI_Lib24_ALL_reduced
% Evaluate the Lib24 reduced-model objective function.
%
% DESCRIPTION
%   Computes the scalar BADS objective for all 24 Lib24 constructs. The
%   parameter vector contains three promoter rates, four RBS intrinsic
%   initiation capacities, and the pGreen/pSC101 copy-number multiplier.
%
%   The non-ribosomal mRNA degradation rate is supplied as model_c.dm_c.
%   Therefore, the same objective function is used for every d_mA case.
%
% COST FUNCTION
%   Residuals are evaluated in base-10 logarithmic synthesis-rate space and
%   penalized with a weighted pseudo-Huber loss:
%
%       r = log10(Pi_pred) - log10(Pi_exp)
%       L = w*delta^2*(sqrt(1 + (r/delta)^2) - 1)
%       w = sqrt(1 + log10(1 + Pi_exp))
%
%   The objective is averaged over growth-rate points, wells or instances,
%   and the 24 constructs.
%
% INPUTS
%   params                  Estimated parameter vector.
%   model_c                 Circuit constants, including d_mA and N_pSC101.
%   Use_mean                'Global', 'Instances', or 'Wells'.
%   ExpData_Tensor_lib30_micro
%                           Experimental tensor used to select the Lib24 subset.
%   HEM                     Host Equivalent Model surrogate.
%   RBS_inv_sigma0          Fixed rho^0 value.
%   delta                   Pseudo-Huber transition parameter.
%
% OUTPUT
%   Sum_prediction_error    Scalar objective value.
%
% This function is deterministic and produces no figures or files.

indices_plasmids_lib24 = [1,2];
indices_promoters_lib24 = [1,2,3];
indices_rbss_lib24 = [1,2,4,5];
num_constructs = length(indices_plasmids_lib24)*length(indices_promoters_lib24)*length(indices_rbss_lib24);
Sum_prediction_error = 0;
% NOTICE each combination (i,j,k) below will correspond to a library transcriptional unit (TU)
for i=1:length(indices_plasmids_lib24) %Plasmids i=1 high copy (pGreen), i=2-> low copy (pSC101)
    for j=1:length(indices_promoters_lib24)
        for k=1:length(indices_rbss_lib24)
            Omega = params(j);
            RBS_k0_sigma0 = params(length(indices_promoters_lib24)+k);
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
                    error = log10(Pi_predictions) - log10(Pi_experimental);
                    rho = delta^2 * (sqrt(1 + (error./delta).^2) - 1);
                    weight = sqrt(1 + log10(1 + Pi_experimental));
                    Sum_prediction_error =  Sum_prediction_error +  sum(weight.*rho)/length(Mu_vector);
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
                             Sum_prediction_error_instance = Sum_prediction_error_instance + Inf;  % Handle NaN case
                             if isinf(Sum_prediction_error_instance)
                                warning('MyWarning:InfiniteResult', 'Pi prediction resulting in Inf in instance  %d of the TU with indexes (%d,%d,%d) for RBS_k0_sigma0=%f, Omega=%f, and Gen_cn=%f.\n', r, i, j,k,RBS_k0_sigma0,Omega,Gen_cn);
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
                                     Sum_prediction_error_wells= Sum_prediction_error_wells + Inf; % Handle NaN case
                                      if isinf(Sum_prediction_error_wells)
                                        warning('MyWarning:InfiniteResult', 'Pi prediction resulting in Inf in well %d of instance  %d of the TU with indexes (%d,%d,%d) for RBS_k0_sigma0=%f, Omega=%f, and Gen_cn=%f.\n', w,r, i, j,k,RBS_k0_sigma0,Omega,Gen_cn);
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
    end %loop for promoter
end %loop for plasmids
Sum_prediction_error =  Sum_prediction_error/num_constructs;
end % of function

