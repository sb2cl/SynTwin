%% Generate_Results_pARKA21_background_model
% Compile pARKA21 background-model estimations and predictions.
%
% INPUTS
%   Experimental_Data/ExpData_pARKA21.mat
%   Estimated_results/Results_BADS_pARKA21_background_Wells.mat
%
% OUTPUT
%   Results_pARKA21_background_Wells.mat

clearvars;
close all;
dbstop if error
warning on

ROOT = init_SynTwin('experimental',true); %#ok<NASGU>
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);

data_file = SynTwin_path('Experimental_Data','ExpData_pARKA21.mat');
D = load(data_file,'Data_pARKA21');
if ~isfield(D,'Data_pARKA21')
    error('Generate_Results_pARKA21_background_model:MissingData', ...
        'Expected Data_pARKA21 in %s.',data_file);
end
Data_pARKA21 = D.Data_pARKA21;

estimation_file = fullfile(SCRIPT_DIR,'Estimated_results', ...
    'Results_BADS_pARKA21_background_Wells.mat');
E = load(estimation_file);
if ~isfield(E,'Results_BADS_pARKA21_background')
    error('Generate_Results_pARKA21_background_model:MissingEstimations', ...
        'Expected Results_BADS_pARKA21_background in %s.',estimation_file);
end

runs = E.Results_BADS_pARKA21_background;
num_runs = numel(runs);
run_matrix = nan(num_runs,6);
for run_idx = 1:num_runs
    if isempty(runs{run_idx}) || ~isfield(runs{run_idx},'results')
        continue;
    end
    row = runs{run_idx}.results;
    if numel(row) ~= 6
        error('Generate_Results_pARKA21_background_model:InvalidRun', ...
            'Run %d does not contain five parameters and one cost.',run_idx);
    end
    run_matrix(run_idx,:) = row(:).';
end
run_matrix = run_matrix(all(isfinite(run_matrix),2),:);
if isempty(run_matrix)
    error('Generate_Results_pARKA21_background_model:NoValidRuns', ...
        'No valid BADS runs were found.');
end

[J_sorted,order] = sort(run_matrix(:,6),'ascend');
params_sorted = run_matrix(order,1:5);
num_best_for_mean = min(5,size(params_sorted,1));

Results_pARKA21_background = struct();
Results_pARKA21_background.Use_mean = 'Wells';
Results_pARKA21_background.Data_pARKA21 = Data_pARKA21;
Results_pARKA21_background.ModelExpression = ...
    'x=100*mu; Pi_bk=b0+A*(1-exp(-x/k))+q1*x+q2*x^2';
Results_pARKA21_background.ParameterNames = {'b0','A','k','q1','q2'};
Results_pARKA21_background.Parameters_raw = run_matrix(:,1:5);
Results_pARKA21_background.J_raw = run_matrix(:,6);
Results_pARKA21_background.Parameters_best = params_sorted(1,:);
Results_pARKA21_background.J_best = J_sorted(1);
Results_pARKA21_background.Parameters_mean_bestN = ...
    mean(params_sorted(1:num_best_for_mean,:),1);
Results_pARKA21_background.NumBestForMean = num_best_for_mean;
Results_pARKA21_background.Parameters_mean = mean(run_matrix(:,1:5),1);
Results_pARKA21_background.Parameters_std = std(run_matrix(:,1:5),0,1);
Results_pARKA21_background.Bounds = struct( ...
    'lower',E.lb,'upper',E.ub, ...
    'plausible_lower',E.plb,'plausible_upper',E.pub);
Results_pARKA21_background.SkipData = E.skip_data;

params_for_prediction = Results_pARKA21_background.Parameters_mean_bestN;
mu_all = [];
num_instances = numel(Data_pARKA21.Stats_Mu.List_instances);
Results_pARKA21_background.Instances = cell(num_instances,1);

for instance_idx = 1:num_instances
    mu = Data_pARKA21.Stats_Mu.List_instances{instance_idx}.Data_mumax_pmax_mean(:);
    pi_exp = Data_pARKA21.Stats_Pi_p.List_instances{instance_idx}.Data_mumax_pmax_mean(:);
    n = min(numel(mu),numel(pi_exp));
    first_idx = max(1,E.skip_data.first);
    last_idx = n-max(0,E.skip_data.last);

    if last_idx >= first_idx
        mu = mu(first_idx:last_idx);
        pi_exp = pi_exp(first_idx:last_idx);
    else
        mu = [];
        pi_exp = [];
    end

    valid = isfinite(mu) & isfinite(pi_exp);
    mu = mu(valid);
    pi_exp = pi_exp(valid);
    pi_pred = background_model(mu,params_for_prediction);

    Results_pARKA21_background.Instances{instance_idx}.Mu = mu;
    Results_pARKA21_background.Instances{instance_idx}.Pi_exp = pi_exp;
    Results_pARKA21_background.Instances{instance_idx}.Pi_pred = pi_pred;
    Results_pARKA21_background.Instances{instance_idx}.Residual = pi_exp-pi_pred;
    Results_pARKA21_background.Instances{instance_idx}.RMSE = ...
        sqrt(mean((pi_exp-pi_pred).^2,'omitnan'));

    source = Data_pARKA21.Stats_Particles.List_instances{instance_idx};
    if isfield(source,'Experiment_IDnum')
        Results_pARKA21_background.Instances{instance_idx}.Experiment_IDnum = ...
            source.Experiment_IDnum;
    end
    if isfield(source,'Experiment_date')
        Results_pARKA21_background.Instances{instance_idx}.Experiment_date = ...
            source.Experiment_date;
    end
    mu_all = [mu_all;mu]; %#ok<AGROW>
end

mu_min = min(mu_all,[],'omitnan');
mu_max = max(mu_all,[],'omitnan');
Results_pARKA21_background.Mu_grid = linspace(mu_min,mu_max,400).';
Results_pARKA21_background.Pi_grid = background_model( ...
    Results_pARKA21_background.Mu_grid,params_for_prediction);

% Backward-compatible parameter vector used by downstream background correction.
Params_global_BKmodel_mean = params_for_prediction; %#ok<NASGU>

output_file = fullfile(SCRIPT_DIR,'Results_pARKA21_background_Wells.mat');
save(output_file,'Results_pARKA21_background', ...
    'Params_global_BKmodel_mean','-mat');

fprintf('[pARKA21] Compiled results saved to %s\n',output_file);

function y = background_model(mu,params)
x = 100 .* max(0,mu);
y = params(1) + params(2).*(1-exp(-x./params(3))) + ...
    params(4).*x + params(5).*x.^2;
end
