%% Analyze_Lib24_reduced_Jacobian
% Analyze practical identifiability of the reduced Lib24 model.
%
% PURPOSE
%   Builds absolute and relative experimental Jacobians of the synthesis-rate
%   output Pi with respect to the free Lib24 parameters under the reduced
%   translation model.
%
% DATA LAYERS
%   Global, Instances, or Wells.
%
% REPORTED QUANTITIES
%   Exact rank, SVD-based effective rank, condition number, least identifiable
%   singular-vector direction, approximate parameter standard deviations,
%   parameter correlations, and the least identifiable parameter combination.
%
% INPUT
%   A Lib24 result tensor containing global, instance-level, and well-level
%   sensitivities. Configure RESULTS_FILE, RESULTS_VAR_NAME, and DATA_LAYER in
%   the USER CONFIGURATION section.
%
% DEFAULT RESULT FILE
%   Estimation_Pi/L24_ALL_reduced_model/
%       Results_Tensor_Lib24_ALL_reduced_model_Wells.mat
%
% DEPENDENCIES
%   init_SynTwin, Load_Library_ALL_TU_INFO, and build_jacobian_reduced.
%
% USAGE
%   Analyze_Lib24_reduced_Jacobian
%
% NOTES
%   The effective-rank threshold is tol = 1e-2*s(1).
%   See README_Jacobian_analysis.md for the complete folder workflow.

clearvars;
close all;
dbstop if error
warning on

%% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin('experimental',true);

% Ensure this script folder is on the path (for local helper functions)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);

%% ---------------------- USER CONFIGURATION ----------------------
% Choose data layer: 'Global' | 'Instances' | 'Wells'
%DATA_LAYER = 'Global';
%DATA_LAYER = 'Instances';
DATA_LAYER = 'Wells';

% Lib24 results file containing sensitivities.
% Adjust the folder/file. The file must contain a cell array Results_Tensor_Lib24_(ALL,L1O)_reduced{ori,prom,rbs} with
% fields:
%   .S_Pi_NA_global_values, .Exp_Data.S_Pi_Omega_global_values, ...
%   .Instances{...}.S_Pi_*_instance_values
%   .Instances{...}.Wells{...}.S_Pi_*_well_values
RESULTS_FILE = SynTwin_path('Estimation_Pi/L24_ALL_reduced_model','Results_Tensor_Lib24_ALL_reduced_model_Wells.mat');

% Variable name inside RESULTS_FILE (edit if different)
RESULTS_VAR_NAME = 'Results_Tensor_Lib24_ALL_reduced';
%% ---------------------------------------------------------------

%% 1) Load library meta-information (internal TU naming)
% Creates Cell_TU_info_ALL{ori,prom,rbs} (and other helpers) used by SynTwin.
Load_Library_ALL_TU_INFO;

%% 2) Load Lib24 results containing sensitivities
if ~exist(RESULTS_FILE,'file')
    error(['Lib24 results file not found:\n  %s\n\n' ...
           'Edit RESULTS_FILE in this script to point to the .mat file ' ...
           'containing the Lib24 sensitivity tensor.'], RESULTS_FILE);
end
S_loaded = load(RESULTS_FILE);
if ~isfield(S_loaded, RESULTS_VAR_NAME)
    error(['Expected variable "%s" not found in:\n  %s\n\n' ...
           'Available variables: %s'], RESULTS_VAR_NAME, RESULTS_FILE, strjoin(fieldnames(S_loaded),', '));
end
Results_Lib24 = S_loaded.(RESULTS_VAR_NAME);
clear S_loaded

%% 3) Define Lib24 biopart indices (Lib24 tensor order)
indices_plasmids_Lib24  = [1,2];        % 1: pGreen, 2: pSC101
indices_promoters_Lib24 = [1,2,3];      % 1: J23106, 2: J23102, 3: J23101
indices_rbss_Lib24      = [1,2,3,4];  % 1: B0030, 2: B0032, 3: J61100, 5: J61101

num_plasmids_Lib24  = numel(indices_plasmids_Lib24);
num_promoters_Lib24 = numel(indices_promoters_Lib24);
num_rbss_Lib24      = numel(indices_rbss_Lib24);

%% 4) TU IDs used in the paper (Lib24 order)
TUs_paper_L24 = { ...
    'L24_13','L24_14','L24_15','L24_16', ...
    'L24_17','L24_18','L24_19','L24_20', ...
    'L24_21','L24_22','L24_23','L24_24', ...
    'L24_01','L24_02','L24_03','L24_04', ...
    'L24_05','L24_06','L24_07','L24_08', ...
    'L24_09','L24_10','L24_11','L24_12'};

TUs_number_L24 = { ...
    13,14,15,16, ...
    17,18,19,20, ...
    21,22,23,24, ...
    1, 2, 3, 4, ...
    5, 6, 7, 8, ...
    9,10,11,12};

IDs_L24 = TUs_paper_L24;

%% 5) Build ORIs/Promoters/RBSs arrays in the same order as TUs_paper
plasmid_names  = {'pGreen','pSC101'};
promoter_names = {'J23106','J23102','J23101'};
rbs_names      = {'B0030','B0032','J61100','J61101'};

nTUs_L24 = num_plasmids_Lib24 * num_promoters_Lib24 * num_rbss_Lib24;
if nTUs_L24 ~= numel(TUs_paper_L24)
    error('Inconsistent TU count: nTUs_L24 (%d) ~= numel(TUs_paper) (%d)', nTUs_L24, numel(TUs_paper_L24));
end

ORIs      = cell(nTUs_L24,1);
Promoters = cell(nTUs_L24,1);
RBSs      = cell(nTUs_L24,1);

Sensitivities_L24_raw = cell(nTUs_L24,1);

counter = 1;
for ii = 1:num_plasmids_Lib24
    for jj = 1:num_promoters_Lib24
        for kk = 1:num_rbss_Lib24

            % Paper identifiers
            S.ID_label = TUs_paper_L24{counter};
            S.ID_num   = TUs_number_L24{counter};

            % Bioparts
            ORIs{counter}      = plasmid_names{ii};
            Promoters{counter} = promoter_names{jj};
            RBSs{counter}      = rbs_names{kk};

            % Extract sensitivities from the Lib24 results tensor
            R = Results_Lib24{ii,jj,kk};
            S.ID_TU_Bioparts = R.TU_Bioparts; %#ok<NASGU> (kept for potential consistency checks)

            % --- Absolute Global  Sensitivities partial_Pi/partial_parameter (one mean trajectory per TU) ---
            S.dPi_dN_global_absolute      = R.S_Pi_NA_global_values;
            S.dPi_domega_global_absolute  = R.S_Pi_Omega_global_values;
            S.dPi_dthetaK_global_absolute = R.S_Pi_RBS_k0_sigma0_global_values;

            % --- Relative Global Sensitivities parameter/Pi*partial_Pi/partial_parameter (one mean trajectory per TU) ---
            S.dPi_dN_global_relative      = R.Gene_cn_mean./R.Pi_mumax_pmax_global_mean.*R.S_Pi_NA_global_values;
            S.dPi_domega_global_relative  = R.Omega_mean./R.Pi_mumax_pmax_global_mean.*R.S_Pi_Omega_global_values;
            S.dPi_dthetaK_global_relative = R.RBS_k0_sigma0_mean./R.Pi_mumax_pmax_global_mean.*R.S_Pi_RBS_k0_sigma0_global_values;

            % --- Instances (concatenate all instance-level means) ---
            S.dPi_dN_instances_absolute      = [];
            S.dPi_domega_instances_absolute  = [];
            S.dPi_dthetaK_instances_absolute = [];
            S.dPi_dN_instances_relative      = [];
            S.dPi_domega_instances_relative  = [];
            S.dPi_dthetaK_instances_relative = [];

            for p = 1:numel(R.Instances)
                inst = R.Instances{p};
                S.dPi_dN_instances_absolute      = [S.dPi_dN_instances_absolute;      inst.S_Pi_NA_instance_values]; %#ok<AGROW>
                S.dPi_domega_instances_absolute  = [S.dPi_domega_instances_absolute;  inst.S_Pi_Omega_instance_values]; %#ok<AGROW>
                S.dPi_dthetaK_instances_absolute = [S.dPi_dthetaK_instances_absolute; inst.S_Pi_RBS_k0_sigma0_instance_values]; %#ok<AGROW>

                S.dPi_dN_instances_relative      = [S.dPi_dN_instances_relative;      R.Gene_cn_mean./inst.Pi_mumax_pmax_instance_mean.*inst.S_Pi_NA_instance_values]; %#ok<AGROW>
                S.dPi_domega_instances_relative  = [S.dPi_domega_instances_relative;  R.Omega_mean./inst.Pi_mumax_pmax_instance_mean.*inst.S_Pi_Omega_instance_values]; %#ok<AGROW>
                S.dPi_dthetaK_instances_relative = [S.dPi_dthetaK_instances_relative; R.RBS_k0_sigma0_mean./inst.Pi_mumax_pmax_instance_mean.*inst.S_Pi_RBS_k0_sigma0_instance_values]; %#ok<AGROW>
            end

            % --- Wells (concatenate all wells across all instances) ---
            S.dPi_dN_wells_absolute      = [];
            S.dPi_domega_wells_absolute  = [];
            S.dPi_dthetaK_wells_absolute = [];
            S.dPi_dN_wells_relative      = [];
            S.dPi_domega_wells_relative  = [];
            S.dPi_dthetaK_wells_relative = [];

            for p = 1:numel(R.Instances)
                inst = R.Instances{p};
                for q = 1:numel(inst.Wells)
                    well = inst.Wells{q};
                    S.dPi_dN_wells_absolute      = [S.dPi_dN_wells_absolute;      well.S_Pi_NA_well_values]; %#ok<AGROW>
                    S.dPi_domega_wells_absolute  = [S.dPi_domega_wells_absolute;  well.S_Pi_Omega_well_values]; %#ok<AGROW>
                    S.dPi_dthetaK_wells_absolute = [S.dPi_dthetaK_wells_absolute; well.S_Pi_RBS_k0_sigma0_well_values]; %#ok<AGROW>

                    S.dPi_dN_wells_relative      = [S.dPi_dN_wells_relative;      R.Gene_cn_mean./well.Pi_mumax_pmax.*well.S_Pi_NA_well_values]; %#ok<AGROW>
                    S.dPi_domega_wells_relative  = [S.dPi_domega_wells_relative;  R.Omega_mean./well.Pi_mumax_pmax.*well.S_Pi_Omega_well_values]; %#ok<AGROW>
                    S.dPi_dthetaK_wells_relative = [S.dPi_dthetaK_wells_relative; R.RBS_k0_sigma0_mean./well.Pi_mumax_pmax.*well.S_Pi_RBS_k0_sigma0_well_values]; %#ok<AGROW>
                end
            end

            Sensitivities_L24_raw{counter} = S;
            counter = counter + 1;
        end
    end
end

%% 6) Convert Sensitivities_L24_raw into the flat structure expected by build_jacobian
Sensitivities_L24 = struct('ID',{},'dPi_dN_absolute',{},'dPi_domega_absolute',{},'dPi_dthetaK_absolute',{}...
                                  ,'dPi_dN_relative',{},'dPi_domega_relative',{},'dPi_dthetaK_relative',{});

switch lower(DATA_LAYER)
    case 'global'
        for k = 1:nTUs_L24
            d = Sensitivities_L24_raw{k};
            Sensitivities_L24(k).ID          = d.ID_label;
            Sensitivities_L24(k).dPi_dN_absolute      = d.dPi_dN_global_absolute;
            Sensitivities_L24(k).dPi_domega_absolute  = d.dPi_domega_global_absolute;
            Sensitivities_L24(k).dPi_dthetaK_absolute = d.dPi_dthetaK_global_absolute;
            Sensitivities_L24(k).dPi_dN_relative      = d.dPi_dN_global_relative;
            Sensitivities_L24(k).dPi_domega_relative  = d.dPi_domega_global_relative;
            Sensitivities_L24(k).dPi_dthetaK_relative = d.dPi_dthetaK_global_relative;
        end
    case 'instances'
        for k = 1:nTUs_L24
            d = Sensitivities_L24_raw{k};
            Sensitivities_L24(k).ID          = d.ID_label;
            Sensitivities_L24(k).dPi_dN_absolute      = d.dPi_dN_instances_absolute;
            Sensitivities_L24(k).dPi_domega_absolute  = d.dPi_domega_instances_absolute;
            Sensitivities_L24(k).dPi_dthetaK_absolute = d.dPi_dthetaK_instances_absolute;
            Sensitivities_L24(k).dPi_dN_relative      = d.dPi_dN_instances_relative;
            Sensitivities_L24(k).dPi_domega_relative  = d.dPi_domega_instances_relative;
            Sensitivities_L24(k).dPi_dthetaK_relative = d.dPi_dthetaK_instances_relative;
        end
    case 'wells'
        for k = 1:nTUs_L24
            d = Sensitivities_L24_raw{k};
            Sensitivities_L24(k).ID          = d.ID_label;
            Sensitivities_L24(k).dPi_dN_absolute      = d.dPi_dN_wells_absolute;
            Sensitivities_L24(k).dPi_domega_absolute  = d.dPi_domega_wells_absolute;
            Sensitivities_L24(k).dPi_dthetaK_absolute = d.dPi_dthetaK_wells_absolute;
            Sensitivities_L24(k).dPi_dN_relative      = d.dPi_dN_wells_relative;
            Sensitivities_L24(k).dPi_domega_relative  = d.dPi_domega_wells_relative;
            Sensitivities_L24(k).dPi_dthetaK_relative = d.dPi_dthetaK_wells_relative;
        end
    otherwise
        error('Invalid DATA_LAYER: %s. Use "Global", "Instances", or "Wells".', DATA_LAYER);
end

%% 7) Build and analyze Jacobian for Lib24
incORIs_L24      = {'pGreen'};   % sólo N_pGreen libre
incProm_L24 = {'J23106','J23102','J23101'};
incRBS_L24       = {'B0030','B0032','J61100','J61101'};

param_biopart_Names = [incORIs_L24 incProm_L24 incRBS_L24];
        

[J_L24_absolute, J_L24_relative, params_L24] = build_jacobian_reduced( ...
    IDs_L24, TUs_paper_L24, ORIs, Promoters, RBSs, ...
    'IncludeORIs',      incORIs_L24, ...
    'IncludePromoters', incProm_L24, ...
    'IncludeRBS',       incRBS_L24, ...
    'Sensitivities',    Sensitivities_L24);

print_rank_summary(J_L24_absolute, params_L24, sprintf('L24 absolute (%s)', DATA_LAYER),param_biopart_Names);

print_rank_summary(J_L24_relative, params_L24, sprintf('L24 relative (%s)', DATA_LAYER),param_biopart_Names);


%% ----------------------------- helpers -----------------------------
function print_rank_summary(J, params, label,param_biopart_Names)

    % Exact rank and effective rank based on singular values
    [~,S,V] = svd(J,'econ');
    s = diag(S);

    r_exact = rank(J);
    condJ = s(1)/s(end);
    v_bad = V(:,end);   % combinación peor identificada

    if isempty(s)
        tol = NaN;
        r_eff = 0;
    else
        tol = 1e-2 * s(1);
        r_eff = sum(s > tol);
    end

    F = J'*J;
    C = inv(F);
    stderr = sqrt(diag(C));
    Corr = C ./ sqrt(diag(C)*diag(C)');

    
    % ----- impresión -----

    fprintf('--- %s ---\n', label);
    fprintf('  #params        : %d\n', numel(params));
    fprintf('  rank(J)        : %d\n', r_exact);
    fprintf('  effective rank : %d (tol = %.2e)\n\n', r_eff, tol);

    fprintf('  conditioning   : %.3g\n', condJ);

    fprintf('  bad direction  : ');
    fprintf('%.3f ', abs(v_bad));
    fprintf('\n\n');

    fprintf('Parameter CV (relative std):\n')

    for i = 1:length(stderr)
        fprintf('  %-10s : %.3g\n', param_biopart_Names{i}, stderr(i));
    end

    fprintf('\nCorrelation matrix:\n');
    disp(Corr)

    % ---- worst identifiable direction ----

fprintf('\nLeast identifiable parameter combination (|V(:,end)| sorted):\n')

[~,idx] = sort(abs(v_bad),'descend');

for k = 1:length(idx)
    i = idx(k);
    fprintf('  %-10s : %+ .3f\n', param_biopart_Names{i}, v_bad(i));
end

end






