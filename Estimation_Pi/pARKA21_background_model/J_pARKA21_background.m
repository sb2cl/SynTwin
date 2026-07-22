function J = J_pARKA21_background(params,Data_pARKA21,skip_data)
%% J_pARKA21_background
% Weighted equally across pARKA21 experimental instances.
%
% MODEL
%   x = 100*mu
%   Pi_bk = b0 + A*(1-exp(-x/k)) + q1*x + q2*x^2
%
% PARAMETER ORDER
%   params = [b0, A, k, q1, q2]
%
% COST
%   Mean squared prediction error is calculated within each experimental
%   instance and then averaged across valid instances. The same leading and
%   trailing points used by the visualization workflow are excluded through
%   skip_data.first and skip_data.last.

arguments
    params (1,:) double
    Data_pARKA21 struct
    skip_data struct = struct('first',4,'last',3)
end

if numel(params) ~= 5 || any(~isfinite(params)) || params(3) <= 0
    J = Inf;
    return;
end

if ~isfield(Data_pARKA21,'Stats_Mu') || ...
        ~isfield(Data_pARKA21.Stats_Mu,'List_instances')
    error('J_pARKA21_background:InvalidData', ...
        'Data_pARKA21 does not contain Stats_Mu.List_instances.');
end

num_instances = numel(Data_pARKA21.Stats_Mu.List_instances);
instance_cost = nan(num_instances,1);

for instance_idx = 1:num_instances
    mu = Data_pARKA21.Stats_Mu.List_instances{instance_idx}.Data_mumax_pmax_mean(:);
    pi_exp = Data_pARKA21.Stats_Pi_p.List_instances{instance_idx}.Data_mumax_pmax_mean(:);

    n = min(numel(mu),numel(pi_exp));
    first_idx = max(1,skip_data.first);
    last_idx = n - max(0,skip_data.last);
    if last_idx < first_idx
        continue;
    end

    mu = mu(first_idx:last_idx);
    pi_exp = pi_exp(first_idx:last_idx);
    valid = isfinite(mu) & isfinite(pi_exp);
    if ~any(valid)
        continue;
    end

    pi_pred = background_model(mu(valid),params);
    if any(~isfinite(pi_pred))
        J = Inf;
        return;
    end
    residual = pi_exp(valid) - pi_pred;
    instance_cost(instance_idx) = mean(residual.^2);
end

if all(~isfinite(instance_cost))
    J = Inf;
else
    J = mean(instance_cost,'omitnan');
end
end

function y = background_model(mu,params)
x = 100 .* max(0,mu);
b0 = params(1);
A = params(2);
k = params(3);
q1 = params(4);
q2 = params(5);
y = b0 + A .* (1-exp(-x./k)) + q1.*x + q2.*x.^2;
end
