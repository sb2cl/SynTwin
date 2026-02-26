function  [varphi_common, theta_common, fs_common]  = Find_HEMLoad_varphi_fs_common(Varphi_exp,Varphi_HEM, Load_exp, f_s)
% ========================================================================
% SynTwin — Find_HEMLoad_varphi_fs_common
% ========================================================================
% PURPOSE
%   Reconcile (match) experimental/circuit-consistent and HEM-consistent
%   resource-flux curves across substrate fraction f_s, and return the
%   common operating point:
%       - varphi_common : reconciled resource flux
%       - theta_common  : reconciled circuit/host load
%       - fs_common     : substrate fraction at the intersection
%
% Used by: Get_synthesis_predictions.m, Get_synthesis_predictions_lite.m
%
% INPUTS
%   Varphi_exp    Vector varphi(f_s) obtained by mapping experimental mu
%                through the HEM surrogate (mu -> varphi) at each f_s.
%   Varphi_HEM    Vector varphi_HEM(f_s) obtained by mapping circuit load
%                theta through the HEM surrogate (theta -> varphi) at each f_s.
%   Load_exp      Vector theta(f_s) circuit load associated with Varphi_exp.
%   f_s           Vector of substrate fraction values (same grid as above).
%
% OUTPUTS
%   varphi_common Reconciled resource flux at the intersection point.
%   theta_common  Reconciled load at the same intersection point.
%   fs_common     Substrate fraction where reconciliation occurs.
%
% METHOD (IMPLEMENTATION NOTE)
%   The intersection is located by minimizing |Varphi_exp - Varphi_HEM| on
%   the valid (non-NaN) subset; a local linear approximation is used to
%   refine the f_s value around the closest point to zero difference.
%
% USAGE CONTEXT
%   Called internally by Get_synthesis_predictions(_lite) to obtain a
%   consistent host–circuit operating point for each experimental mu.
%
% LICENSE / CITATION
%   Part of SynTwin. Please cite the SynTwin paper/software documentation.
% ========================================================================


% Takes the experimental Varphi_exp (parametrized by f_s) with  the compatible circuit fluxes to a given Mu_exp, and the corresponding circuit loads Load_exp (parametrized by f_s).
% Takes the vector  Varphi_HEM (parametrized by f_s) with the compatible HEM fluxes to  Load_exp.
% Calculates the intersection point between Varphi_exp and Varphi_HEM for some value of the substrate  f_s. This gives the actual flux of resources varphi_common
% (a single reconciliated value) and the corresponding load theta_common at the single common value of f_s (fs_common) corresponding to that experimental Mu_exp.
 
% We first obtain the difference between 
difference_HEM_exp = Varphi_exp - Varphi_HEM;
% We obtain the subvector with no Nan elements
[~,indices] = find(not(isnan(difference_HEM_exp)));
vector_diff = difference_HEM_exp(indices);

% Get the indices of the closest values to zero (intersection point):
index_closest=find( min(abs(vector_diff)) == abs(vector_diff) );
f_s_restricted = f_s(indices);
if not(length(f_s_restricted)==1)
    if index_closest==1 %the intersection either is on the left of the f_s values or between the first and second values
        index = index_closest;
        index_second = index_closest + 1;
        slope = (vector_diff(index_second)- vector_diff(index))/( f_s_restricted(index_second) - f_s_restricted(index) );
    elseif index_closest==length(vector_diff) %the intersection either is on the right of the f_s values or between the last and one-to-last values
        index = index_closest -1;
        index_second = index_closest;
        slope = (vector_diff(index_second)- vector_diff(index))/( f_s_restricted(index_second) - f_s_restricted(index) );
    else
       index = index_closest;
       index_first = index_closest-1;
       index_second = index_closest+1;
       slope_before = (vector_diff(index)- vector_diff(index_first))/( f_s_restricted(index) - f_s_restricted(index_first) );
       slope_after = (vector_diff(index_second)- vector_diff(index))/( f_s_restricted(index_second) - f_s_restricted(index) );
       slope = 0.5*slope_before+0.5*slope_after;
    end
    % Now we find the intersection with zero using the linear approximation 
    % varphi_predicted(fs(x)) = f(index) + slope*(fs(x)-fs(index)) -> fs(x)_zero = fs(index) - f(index)/slope
    fs_zero = f_s_restricted(index) - vector_diff(index)/slope;
    if (fs_zero>1)&&(fs_zero-1< 0.01) %add a tolerance
        fs_zero = 1;
    elseif (fs_zero<0)&&(fs_zero> -0.01)
        fs_zero = 0;
    end
else
    fs_zero = f_s_restricted(index_closest);
end

% Now we see if the intersection falls between 0\leq fs_zero \leq 1
if not(isempty(fs_zero))&&(fs_zero>=0)&&(fs_zero<=1)
    fs_common = fs_zero;
    % Get the closest left and right values of mu HEM_mu_vector
    index_low =  find(f_s_restricted<=fs_zero,1,'last');
    index_high = find(f_s_restricted>=fs_zero,1,'first');
    if isempty(index_low) %We are in the first element
         index_low = index_high;
    elseif isempty(index_high) %We are in the last element
        index_high = index_low;
    end
    Varphi_exp_restricted = Varphi_exp(indices);
    Load_exp_restricted = Load_exp(indices);
    % Obtain the corresponding values of varphi:
    varphi_low = Varphi_exp_restricted(index_low);
    varphi_high = Varphi_exp_restricted(index_high);
    theta_low = Load_exp_restricted(index_low);
    theta_high = Load_exp_restricted(index_high);
    % Get the interpolated value using Vq = interp1(X,V,Xq) 
        if (varphi_low==varphi_high)
            varphi_common = varphi_low;
        else
           varphi_common = interp1([f_s_restricted(index_low), f_s_restricted(index_high)],[varphi_low,varphi_high], fs_common);
        end
          if (theta_low==theta_high)
            theta_common = theta_low;
        else
           theta_common = interp1([f_s_restricted(index_low), f_s_restricted(index_high)],[theta_low,theta_high], fs_common);
        end
  
else %there is no feasible solution
    fs_common = NaN;
    varphi_common = NaN;
    theta_common = NaN;
end

end
