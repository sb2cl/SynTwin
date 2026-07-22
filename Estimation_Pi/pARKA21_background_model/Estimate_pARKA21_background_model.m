%% Estimate_pARKA21_background_model
% Estimate the pARKA21 background model using BADS multistart.
%
% MODEL
%   x = 100*mu
%   Pi_bk = b0 + A*(1-exp(-x/k)) + q1*x + q2*x^2
%
% PARAMETER ORDER
%   [b0, A, k, q1, q2]
%
% INPUT
%   Experimental_Data/ExpData_pARKA21.mat
%
% OUTPUT
%   Estimated_results/Results_BADS_pARKA21_background_Wells.mat
%
%   The raw Cytation Excel files are not required by this distributed script.

clearvars;
close all;
dbstop if error
warning on

ROOT = init_SynTwin('experimental',true,'bads',true); %#ok<NASGU>
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);

data_file = SynTwin_path('Experimental_Data','ExpData_pARKA21.mat');
S = load(data_file,'Data_pARKA21');
if ~isfield(S,'Data_pARKA21')
    error('Estimate_pARKA21_background_model:MissingData', ...
        'Expected Data_pARKA21 in %s.',data_file);
end
Data_pARKA21 = S.Data_pARKA21;

num_runs = 100;
skip_data = struct('first',4,'last',3);

% Bounds for [b0, A, k, q1, q2].
lb  = [0.0,  20.0, 0.02, -50.0, 0.0];
ub  = [4.0, 140.0, 2.00, -20.0,10.0];
plb = [0.5,  40.0, 0.50, -40.0, 2.0];
pub = [2.0,  90.0, 1.80, -28.0, 8.0];

options = bads('defaults');
options.Display = 'final';

Results_BADS_pARKA21_background = cell(num_runs,1);

parfor run_idx = 1:num_runs
    x0 = plb + (pub-plb).*rand(1,numel(lb));
    objective = @(parameters) J_pARKA21_background( ...
        parameters,Data_pARKA21,skip_data);
    [params,Jmin] = bads(objective,x0,lb,ub,plb,pub,[],options);

    run_result = struct();
    run_result.results = [params,Jmin];
    run_result.x0 = x0;
    Results_BADS_pARKA21_background{run_idx} = run_result;
end

results_dir = fullfile(SCRIPT_DIR,'Estimated_results');
if ~exist(results_dir,'dir')
    mkdir(results_dir);
end

output_file = fullfile(results_dir, ...
    'Results_BADS_pARKA21_background_Wells.mat');

ModelExpression = ...
    'x=100*mu; Pi_bk=b0+A*(1-exp(-x/k))+q1*x+q2*x^2';
ParameterNames = {'b0','A','k','q1','q2'};

save(output_file, ...
    'Results_BADS_pARKA21_background', ...
    'ModelExpression','ParameterNames', ...
    'lb','ub','plb','pub','skip_data','-mat');

fprintf('[pARKA21] Estimation results saved to %s\n',output_file);
