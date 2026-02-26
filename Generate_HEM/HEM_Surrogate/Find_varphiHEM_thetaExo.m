function  varphi_HEM  = Find_varphiHEM_thetaExo(HEM_varphi_vector,Theta_HEM, theta_exo)
% ========================================================================
% SynTwin — Find_varphiHEM_thetaExo
% ========================================================================
% PURPOSE
%   Interpolate the HEM surrogate mapping to obtain varphi corresponding to
%   an exogenous (circuit) load theta_exo at a given substrate fraction f_s:
%       theta_exo  ->  varphi_HEM
%
% Used by: Get_synthesis_predictions.m, Get_synthesis_predictions_lite.m
%
% INPUTS
%   HEM_varphi_vector  Vector of varphi values along the HEM surrogate curve
%                      at fixed f_s.
%   Theta_HEM          Vector of theta load values used to evaluate the HEM
%                      surrogate (typically ordered from low to high).
%   theta_exo          Scalar exogenous/circuit load value.
%
% OUTPUTS
%   varphi_HEM         Scalar interpolated varphi value (NaN if out of range).
%
% USAGE CONTEXT
%   Used to reconcile host and circuit constraints in the digital twin:
%   given a circuit load, find the compatible host resource flux.
%
% LICENSE / CITATION
%   Part of SynTwin. Please cite the SynTwin paper/software documentation.
% ========================================================================

% First check that theta_exo lies in Theta_HEM
% Notice Theta_HEM is ordered from low to high values (from no loading
% to high loading)
if (theta_exo<Theta_HEM(1))||(theta_exo>Theta_HEM(end))||isnan(theta_exo)
    varphi_HEM=NaN;
    return
end

% Get the closest left and right values of Theta_HEM
index_low =  find(Theta_HEM<=theta_exo,1,'last');
index_high = find(Theta_HEM>=theta_exo,1,'first');
if isempty(index_low) %We are in the first element
     index_low = index_high;
elseif isempty(index_high) %We are in the last element
    index_high = index_low;
end
% Obtain the corresponding values of mu in HEM_mu_vector:
varphi_HEM_low = HEM_varphi_vector(index_low);
varphi_HEM_high = HEM_varphi_vector(index_high);
% Get the interpolated value using Vq = interp1(X,V,Xq) 
    if (varphi_HEM_low==varphi_HEM_high)
        varphi_HEM = varphi_HEM_low;
    else
       varphi_HEM = interp1([Theta_HEM(index_low), Theta_HEM(index_high)],[varphi_HEM_low,varphi_HEM_high], theta_exo);
    end
end
