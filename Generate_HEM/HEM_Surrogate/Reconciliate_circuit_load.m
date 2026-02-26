function [varphi_common, theta_circuit] = Reconciliate_circuit_load(HEM_varphi_vector,HEM_thetas, Load_theta)
% ========================================================================
% SynTwin — Reconciliate_circuit_load
% ========================================================================
% PURPOSE
%   Reconcile host and circuit constraints to obtain a consistent operating
%   point (varphi_common, theta_circuit) when the circuit load is provided
%   as a function of varphi (Load_theta). This routine iteratively searches
%   for a varphi value such that:
%       theta_circuit(varphi) matches the host theta(varphi) relation.
%
% INPUTS
%   HEM_varphi_vector  Vector of varphi values along the HEM surrogate curve
%                      at fixed substrate fraction f_s.
%   HEM_thetas         Vector of theta values corresponding to HEM_varphi_vector.
%   Load_theta         Vector theta(varphi) describing the exogenous circuit load
%                      sampled on the same varphi grid (or compatible indexing).
%
% OUTPUTS
%   varphi_common      Reconciled resource flux.
%   theta_circuit      Reconciled circuit load.
%
% USAGE CONTEXT
%   Mainly used in forward simulation settings (when mu is not imposed from
%   experiment). Kept in SynTwin for internal consistency and extensibility.
%
% LICENSE / CITATION
%   Part of SynTwin. Please cite the SynTwin paper/software documentation.
% ========================================================================

% Assume initial value of varphi:
% Given the shape of the Host theta-varphi function, we start from a low
% value of varphi_common
varphi_common = 0.05;

while true 
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
    % Obtain the corresponding values of load theta in Load_theta:
    theta_circuit_low = Load_theta(index_low);
    theta_circuit_high = Load_theta(index_high);
    % Get the interpolated value using Vq = interp1(X,V,Xq) 
    if (theta_circuit_low==theta_circuit_high)
        theta_circuit = theta_circuit_low;
    else
        theta_circuit = interp1([HEM_varphi_vector(index_low), HEM_varphi_vector(index_high)],[theta_circuit_low,theta_circuit_high],varphi_common);
    end
    %Get where (indices) the Host achieves that theta. Notice HEM_thetas
    %goes from zero to high
    index_low = find(HEM_thetas<=theta_circuit,1,'last');
    index_high =  find(HEM_thetas>=theta_circuit,1,'first');
    if isempty(index_low) %We are in the first element
         index_low = index_high;
    elseif isempty(index_high) %We are in the last element
        index_high = index_low;
    end
    % Get the corresponding values of varphi of the host
    varphi_host_low = HEM_varphi_vector(index_low);
    varphi_host_high = HEM_varphi_vector(index_high);
    if (varphi_host_low==varphi_host_high)
        varphi_common_new = varphi_host_low;
    else
        if  (theta_circuit_low==theta_circuit_high)
            varphi_common_new = (varphi_host_low+varphi_host_high)/2;
        else
            % Get the interpolated value
            varphi_common_new = interp1([theta_circuit_low,theta_circuit_high],[varphi_host_low,varphi_host_high],theta_circuit);
        end
    end
   % Obtain an error measure 
    if abs(varphi_common_new - varphi_common) < 1e-4
        varphi_common = varphi_common_new;
        break;
    else
        varphi_common = varphi_common_new; %Iterate with the new estimated varphi 
    end
end % while

end