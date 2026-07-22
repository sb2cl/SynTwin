function filtered_time_matrix = Matrix_Logy_filter(Input_time_series_matrix,interval_eval,interval_fix)
%
% To avoid problems caused by a sometimes wrong first measurements from the
% cytation, we use a simple procedure to fix this: 
%
% Input_time_series_matrix is the sampled data 
% for the replicas of the culture. Each column is a replica. Rows are time samples.
%
arguments
    Input_time_series_matrix double
    interval_eval double = [2,4]
    interval_fix double = [1,3] 
end
num_cultures=size(Input_time_series_matrix,2);
filtered_time_matrix=Input_time_series_matrix;
for i=1:num_cultures
   if filtered_time_matrix(3,i) < filtered_time_matrix(2,i) 
        filtered_time_matrix(2,i)  =filtered_time_matrix(3,i);
        filtered_time_matrix(3,i) = 0.65*filtered_time_matrix(2,i)+0.35*filtered_time_matrix(4,i);
   end 
   if filtered_time_matrix(2,i) < filtered_time_matrix(1,i) 
        filtered_time_matrix(1,i)  = filtered_time_matrix(2,i);
        filtered_time_matrix(2,i) = 0.65*filtered_time_matrix(1,i)+0.35*filtered_time_matrix(3,i);
   end
end

end

