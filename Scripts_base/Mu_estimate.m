function [estimated_mu] = Mu_estimate(Input_time_series,sampling_time,wn)
%Estimates the specific growth rate
% We  obtain the derivative as
%   dy(k)/dt = (-y(k+2) + 8*y(k+1)-8*y(k-1) + y(k-2))/12/dt  
% and then apply a Butterworth filter:
%   [b,a] = butter(2,4/33);
%   data_ff = filtfilt(b,a,data);
% Input_time_series is the sampled data of absorbance or Particles
% Sampling is assumed to be in minutes
% estimated_mu gives the estimate of the specific growth rate in mins^{-1}
%
% We first extend the vector on the left and right extremes:
% We first extend the vector on the left and right extremes:
poly =  polyfit([1,2,3]*sampling_time,Input_time_series(1:3),1);
xm1_value = polyval(poly,0);
xmm1_value = polyval(poly,-1*sampling_time);
poly =  polyfit([1,2,3]*sampling_time,Input_time_series(end-2:end),1);
xM1_value = polyval(poly,4*sampling_time);
xMM1_value = polyval(poly,5*sampling_time);

temp_data= [xmm1_value;xm1_value;Input_time_series;xM1_value;xMM1_value];
der_temp_data = (-temp_data(5:end,:)+8*temp_data(4:end-1,:)-8*temp_data(2:end-3,:)+temp_data(1:end-4,:))./(12*sampling_time); %in mins^{-1}
raw_estimated_mu = der_temp_data./Input_time_series;
raw_estimated_mu_clean = smoothdata(raw_estimated_mu,"rlowess",4);

time_series_samples = 0:sampling_time:sampling_time*(length(raw_estimated_mu_clean)-1);
smoother_window = floor(length(time_series_samples)*wn);
estimated_mu= smoothdata(raw_estimated_mu_clean,'gaussian',smoother_window,'SamplePoints',time_series_samples);

end
