function  varphi_circuit  = Find_varphi(HEM_mu_vector,HEM_varphi_vector, mu_common)
% ========================================================================
% SynTwin — Find_varphi
% ========================================================================
% PURPOSE
%   Interpolate the HEM surrogate mapping to obtain the resource flux varphi
%   corresponding to a given growth rate mu_common:
%       mu_common  ->  varphi_circuit
%
% Used by: Get_synthesis_predictions.m, Get_synthesis_predictions_lite.m
%
% INPUTS
%   HEM_mu_vector      Vector of mu values along the HEM surrogate curve
%                      at a fixed substrate fraction f_s (typically ordered
%                      from high to low as load increases).
%   HEM_varphi_vector  Vector of varphi values corresponding to HEM_mu_vector.
%   mu_common          Scalar growth-rate value for which varphi is requested.
%
% OUTPUTS
%   varphi_circuit     Scalar interpolated varphi value (NaN if mu is out of range).
%
% NOTES
%   - This function is used in "digital twin prediction mode", where mu is
%     taken from experimental measurements.
%   - Out-of-range mu_common returns NaN (guard for extrapolation).
%
% LICENSE / CITATION
%   Part of SynTwin. Please cite the SynTwin paper/software documentation.
% ========================================================================

% First check that mu_common lies in HEM_mu_vector
% Notice HEM_mu_vector is ordered from high to low values (from no loading
% to high loading)
if (mu_common>HEM_mu_vector(1))||(mu_common<HEM_mu_vector(end))
    varphi_circuit=NaN;
    return
end

% Get the closest left and right values of mu HEM_mu_vector
index_low =  find(HEM_mu_vector>=mu_common,1,'last');
index_high = find(HEM_mu_vector<=mu_common,1,'first');
if isempty(index_low) %We are in the first element
     index_low = index_high;
elseif isempty(index_high) %We are in the last element
    index_high = index_low;
end
% Obtain the corresponding values of mu in HEM_mu_vector:
varphi_circuit_low = HEM_varphi_vector(index_low);
varphi_circuit_high = HEM_varphi_vector(index_high);
% Get the interpolated value using Vq = interp1(X,V,Xq) 
    if (varphi_circuit_low==varphi_circuit_high)
        varphi_circuit = varphi_circuit_low;
    else
       varphi_circuit = interp1([HEM_mu_vector(index_low), HEM_mu_vector(index_high)],[varphi_circuit_low,varphi_circuit_high], mu_common);
    end
end
