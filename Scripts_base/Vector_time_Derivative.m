function [derivative_time_series_matrix,mean_derivative_time_series,std_derivative_time_series] = Vector_time_Derivative(Input_time_series_matrix,sampling_time,wn)
%
% Estimates the time derivative for a matrix containing several column
% vectors of sampled data. 
% Uses the function derivative_time_series = Time_Derivative(Input_time_series,sampling_time)
%
% Input_time_series_matrix is the sampled data 
% for several cultures. Each column is a culture. Rows are time samples.
% Sampling is assumed to be in minutes
% derivative_time_series_matrix gives the estimate of the derivative in
% mins^{-1} for each culture
% mean_derivative_time_series and std_derivative_time_series give the mean
% and std 
%
num_cultures=size(Input_time_series_matrix,2);
time_samples= size(Input_time_series_matrix,1);
derivative_time_series_matrix=NaN*ones(time_samples,num_cultures);
for i=1:num_cultures
   derivative_time_series_matrix(:,i)=   Time_Derivative(Input_time_series_matrix(:,i),sampling_time,wn);
end
if nargout == 2
    mean_derivative_time_series = mean(derivative_time_series_matrix,2);
elseif nargout == 3
    mean_derivative_time_series = mean(derivative_time_series_matrix,2);
    std_derivative_time_series = std(mean_derivative_time_series,0,2);
end
end
