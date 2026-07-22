function filtered_time_matrix = Matrix_time_Filter(Input_time_series_matrix,time_series_samples,wn)
%
% Filters a matrix containing several column  vectors of sampled data. 
%
% Input_time_series_matrix is the sampled data 
% for several cultures. Each column is a culture. Rows are time samples.
%
% wn is the fraction of time-samples that will be used as window length by
% the gaussian smoother

%MODIFIED 29/05/2025 to include rlowess method to filter out outliers
num_cultures=size(Input_time_series_matrix,2);
time_samples= size(Input_time_series_matrix,1);
filtered_time_matrix=NaN*ones(time_samples,num_cultures);
%[b,a] = butter(2,wn); 
smoother_window = floor(length(time_series_samples)*wn);
for i=1:num_cultures
   Input_time_series_clean = smoothdata(Input_time_series_matrix(:,i),"rlowess",4);
   filtered_time_matrix(:,i) = smoothdata(Input_time_series_clean,'gaussian',smoother_window,'SamplePoints',time_series_samples);
 % OLD VERSION WITH NO oUTLIER DETECTION:  filtered_time_matrix(:,i) = smoothdata(Input_time_series_matrix(:,i),'gaussian',smoother_window,'SamplePoints',time_series_samples);
end

end


