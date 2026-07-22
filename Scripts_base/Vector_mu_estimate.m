function [estimated_mu_matrix,mean_estimated_mu,std_estimated_mu] = Vector_mu_estimate(Input_time_series_matrix,sampling_time,wn)
%Estimates the specific growth rate for a matrix containing several column
%vectors of sampled data of absorbance or Particles
% Uses the function estimated_mu= Mu_estimate(Input_time_series,sampling_time)
%
% Input_time_series_matrix is the sampled data of absorbance or Particles
% for several cultures. Each column is a culture. Rows are time samples.
% Sampling is assumed to be in minutes
% estimated_mu_matrix gives the estimate of the specific growth rate in
% mins^{-1} for each culture 
% mean_estimated_mu and std_estimated_mu give the mean and std mu time profile of the set
% of cultures
num_cultures=size(Input_time_series_matrix,2);
time_samples= size(Input_time_series_matrix,1);
estimated_mu_matrix=NaN*ones(time_samples,num_cultures);
for i=1:num_cultures
   estimated_mu_matrix(:,i)=   Mu_estimate(Input_time_series_matrix(:,i),sampling_time,wn);
end
if nargout == 2
    mean_estimated_mu = mean(estimated_mu_matrix,2);
elseif nargout == 3
    mean_estimated_mu = mean(estimated_mu_matrix,2);
    std_estimated_mu = std(estimated_mu_matrix,0,2);
end
end