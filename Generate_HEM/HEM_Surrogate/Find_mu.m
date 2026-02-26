function  mu_circuit  = Find_mu(HEM_mu_vector,HEM_varphi_vector, varphi_common)
% ========================================================================
% SynTwin — Find_mu
% ========================================================================
% PURPOSE
%   Interpolate the HEM surrogate mapping to obtain the growth rate mu
%   corresponding to a given resource flux varphi_common:
%       varphi_common  ->  mu_circuit
%
% INPUTS
%   HEM_mu_vector      Vector of mu values along the HEM surrogate curve
%                      at fixed substrate fraction f_s.
%   HEM_varphi_vector  Vector of varphi values corresponding to HEM_mu_vector
%                      (typically ordered from high to low as load increases).
%   varphi_common      Scalar resource-flux value.
%
% OUTPUTS
%   mu_circuit         Scalar interpolated growth rate.
%
% USAGE CONTEXT
%   Primarily used in forward simulation workflows (when mu is not imposed
%   experimentally). Included for completeness in this SynTwin release.
%
% LICENSE / CITATION
%   Part of SynTwin. Please cite the SynTwin paper/software documentation.
% ========================================================================

% Get the closest left and right values of varphi in HEM_varphi_vector
    % Notice HEM_varphi_vector is ordered from high to low values (from no loading
    % to high loading)
    index_high = find(HEM_varphi_vector<=varphi_common,1,'first');
    index_low =  find(HEM_varphi_vector>=varphi_common,1,'last');
    if isempty(index_low) %We are in the first element
         index_low = index_high;
    elseif isempty(index_high) %We are in the last element
        index_high = index_low;
    end
  % Obtain the corresponding values of mu in HEM_mu_vector:
    mu_circuit_low = HEM_mu_vector(index_low);
    mu_circuit_high = HEM_mu_vector(index_high);
% Get the interpolated value using Vq = interp1(X,V,Xq) 
    if (mu_circuit_low==mu_circuit_high)
        mu_circuit = mu_circuit_low;
    else
       mu_circuit = interp1([HEM_varphi_vector(index_low), HEM_varphi_vector(index_high)],[mu_circuit_low,mu_circuit_high],varphi_common);
    end
end