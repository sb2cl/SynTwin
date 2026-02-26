%% Analyze_Libs_Jacobian_ranks
% SynTwin script: practical-identifiability analysis via Jacobian rank.
%
% PURPOSE
%   Computes the (exact) rank and an SVD-based effective rank of the
%   experimental Jacobian of the synthesis-rate model output (Pi) with respect
%   to the free parameters of three libraries:
%       - Lib30
%       - Lib24
%       - Lib6
%
%   The Jacobian is built from the *absolute sensitivities* dPi/dtheta stored
%   in a Lib30 results tensor (Global / Instances / Wells). The analysis
%   provides a practical identifiability proxy: if the Jacobian is rank
%   deficient, some parameter combinations are locally unidentifiable given the
%   available data layer.
%
% DATA LAYERS
%   Select which experimental aggregation layer to use when constructing the
%   Jacobian rows:
%       - 'Global'    : one trajectory per TU (global mean)
%       - 'Instances' : concatenated per-experiment means
%       - 'Wells'     : concatenated per-well trajectories (highest detail)
%
% INPUTS
%   None (configuration is set inside the script).
%   Edit the section "USER CONFIGURATION" below.
%
% OUTPUTS
%   Prints summary statistics to the MATLAB command window:
%       - number of free parameters
%       - exact rank(J)
%       - effective rank(J) based on singular values
%
% DEPENDENCIES
%   - SynTwin must be initialized before running this script:
%       ROOT = init_SynTwin('experimental',true);
%   - This script expects the helper functions used to build the Jacobian to be
%     accessible on the MATLAB path (typically via SynTwin initialization and/or
%     the local folder):
%       * Load_Library_ALL_TU_INFO
%       * get_free_params_for_library
%       * build_jacobian
%
%   - Results file providing Lib30 sensitivities.
%     By default this script expects a variable named:
%       Results_Tensor_Lib30_L1O_reduced
%     stored in a .mat file.
%
% USAGE
%   Analyze_Libs_Jacobian_ranks
%
% NOTES
%   - The effective rank is computed using tol = 1e-2 * s(1), where s(1) is the
%     largest singular value.
%   - If you distribute SynTwin publicly, ensure that the referenced results
%     file is included (or adjust RESULTS_FILE below accordingly).

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
DATA_LAYER = 'Wells';

% Lib30 results file containing sensitivities.
% Adjust the folder/file if your distribution uses a different naming.
% The file must contain a cell array Results_Tensor_Lib30_L1O_reduced{ori,prom,rbs} with
% fields:
%   .S_Pi_NA_global_values, .Exp_Data.S_Pi_Omega_global_values, ...
%   .Instances{...}.S_Pi_*_instance_values
%   .Instances{...}.Wells{...}.S_Pi_*_well_values
RESULTS_FILE = SynTwin_path('Estimation_Pi/L30_L1O_reduced_model','Results_Tensor_Lib30_L1O_reduced_Wells.mat');

% Variable name inside RESULTS_FILE (edit if different)
RESULTS_VAR_NAME = 'Results_Tensor_Lib30_L1O_reduced';
%% ---------------------------------------------------------------

%% 1) Load library meta-information (internal TU naming)
% Creates Cell_TU_info_ALL{ori,prom,rbs} (and other helpers) used by SynTwin.
Load_Library_ALL_TU_INFO;

%% 2) Load Lib30 results containing sensitivities
if ~exist(RESULTS_FILE,'file')
    error(['Lib30 results file not found:\n  %s\n\n' ...
           'Edit RESULTS_FILE in this script to point to the .mat file ' ...
           'containing the Lib30 sensitivity tensor.'], RESULTS_FILE);
end
S_loaded = load(RESULTS_FILE);
if ~isfield(S_loaded, RESULTS_VAR_NAME)
    error(['Expected variable "%s" not found in:\n  %s\n\n' ...
           'Available variables: %s'], RESULTS_VAR_NAME, RESULTS_FILE, strjoin(fieldnames(S_loaded),', '));
end
Results_Lib30_L1O_wells = S_loaded.(RESULTS_VAR_NAME);
clear S_loaded

%% 3) Define Lib30 biopart indices (Lib30 tensor order)
indices_plasmids_Lib30  = [1,2];        % 1: pGreen, 2: pSC101
indices_promoters_Lib30 = [1,2,3];      % 1: J23106, 2: J23102, 3: J23101
indices_rbss_Lib30      = [1,2,3,4,5];  % 1: B0030, 2: B0032, 3: B0034, 4: J61100, 5: J61101

num_plasmids_Lib30  = numel(indices_plasmids_Lib30);
num_promoters_Lib30 = numel(indices_promoters_Lib30);
num_rbss_Lib30      = numel(indices_rbss_Lib30);

%% 4) TU IDs used in the paper (Lib30 order = Lib24 + extensions)
TUs_paper = { ...
    'L24_13','L24_14','L30_28','L24_15','L24_16', ...
    'L24_17','L24_18','L30_29','L24_19','L24_20', ...
    'L24_21','L24_22','L30_30','L24_23','L24_24', ...
    'L24_01','L24_02','L30_25','L24_03','L24_04', ...
    'L24_05','L24_06','L30_26','L24_07','L24_08', ...
    'L24_09','L24_10','L30_27','L24_11','L24_12'};

TUs_number = { ...
    13,14,28,15,16, ...
    17,18,29,19,20, ...
    21,22,30,23,24, ...
    1, 2, 25,3, 4, ...
    5, 6, 26,7, 8, ...
    9,10,27,11,12};

IDs_L30 = TUs_paper;

%% 5) Build ORIs/Promoters/RBSs arrays in the same order as TUs_paper
plasmid_names  = {'pGreen','pSC101'};
promoter_names = {'J23106','J23102','J23101'};
rbs_names      = {'B0030','B0032','B0034','J61100','J61101'};

nTUs_L30 = num_plasmids_Lib30 * num_promoters_Lib30 * num_rbss_Lib30;
if nTUs_L30 ~= numel(TUs_paper)
    error('Inconsistent TU count: nTUs_L30 (%d) ~= numel(TUs_paper) (%d)', nTUs_L30, numel(TUs_paper));
end

ORIs      = cell(nTUs_L30,1);
Promoters = cell(nTUs_L30,1);
RBSs      = cell(nTUs_L30,1);

Sensitivities_L30_raw = cell(nTUs_L30,1);

counter = 1;
for ii = 1:num_plasmids_Lib30
    for jj = 1:num_promoters_Lib30
        for kk = 1:num_rbss_Lib30

            % Paper identifiers
            S.ID_label = TUs_paper{counter};
            S.ID_num   = TUs_number{counter};

            % Bioparts
            ORIs{counter}      = plasmid_names{ii};
            Promoters{counter} = promoter_names{jj};
            RBSs{counter}      = rbs_names{kk};

            % Extract sensitivities from Lib30 results tensor
            R = Results_Lib30_L1O_wells{ii,jj,kk};
            S.ID_TU_Bioparts = R.TU_Bioparts; %#ok<NASGU> (kept for potential consistency checks)

            % --- Global (one mean trajectory per TU) ---
            S.dPi_dN_global      = R.S_Pi_NA_global_values;
            S.dPi_domega_global  = R.S_Pi_Omega_global_values;
            S.dPi_dthetaK_global = R.S_Pi_RBS_k0_sigma0_global_values;
            S.dPi_dthetaS_global = R.S_Pi_RBS_inv_sigma0_global_values;

            % --- Instances (concatenate all instance-level means) ---
            S.dPi_dN_instances      = [];
            S.dPi_domega_instances  = [];
            S.dPi_dthetaK_instances = [];
            S.dPi_dthetaS_instances = [];

            for p = 1:numel(R.Instances)
                inst = R.Instances{p};
                S.dPi_dN_instances      = [S.dPi_dN_instances;      inst.S_Pi_NA_instance_values]; %#ok<AGROW>
                S.dPi_domega_instances  = [S.dPi_domega_instances;  inst.S_Pi_Omega_instance_values]; %#ok<AGROW>
                S.dPi_dthetaK_instances = [S.dPi_dthetaK_instances; inst.S_Pi_RBS_k0_sigma0_instance_values]; %#ok<AGROW>
                S.dPi_dthetaS_instances = [S.dPi_dthetaS_instances; inst.S_Pi_RBS_inv_sigma0_instance_values]; %#ok<AGROW>
            end

            % --- Wells (concatenate all wells across all instances) ---
            S.dPi_dN_wells      = [];
            S.dPi_domega_wells  = [];
            S.dPi_dthetaK_wells = [];
            S.dPi_dthetaS_wells = [];

            for p = 1:numel(R.Instances)
                inst = R.Instances{p};
                for q = 1:numel(inst.Wells)
                    well = inst.Wells{q};
                    S.dPi_dN_wells      = [S.dPi_dN_wells;      well.S_Pi_NA_well_values]; %#ok<AGROW>
                    S.dPi_domega_wells  = [S.dPi_domega_wells;  well.S_Pi_Omega_well_values]; %#ok<AGROW>
                    S.dPi_dthetaK_wells = [S.dPi_dthetaK_wells; well.S_Pi_RBS_k0_sigma0_well_values]; %#ok<AGROW>
                    S.dPi_dthetaS_wells = [S.dPi_dthetaS_wells; well.S_Pi_RBS_inv_sigma0_well_values]; %#ok<AGROW>
                end
            end

            Sensitivities_L30_raw{counter} = S;
            counter = counter + 1;
        end
    end
end

%% 6) Convert Sensitivities_L30_raw into the flat structure expected by build_jacobian
Sensitivities_L30 = struct('ID',{},'dPi_domega',{},'dPi_dthetaK',{},'dPi_dthetaS',{},'dPi_dN',{});

switch lower(DATA_LAYER)
    case 'global'
        for k = 1:nTUs_L30
            d = Sensitivities_L30_raw{k};
            Sensitivities_L30(k).ID          = d.ID_label;
            Sensitivities_L30(k).dPi_domega  = d.dPi_domega_global;
            Sensitivities_L30(k).dPi_dthetaK = d.dPi_dthetaK_global;
            Sensitivities_L30(k).dPi_dthetaS = d.dPi_dthetaS_global;
            Sensitivities_L30(k).dPi_dN      = d.dPi_dN_global;
        end
    case 'instances'
        for k = 1:nTUs_L30
            d = Sensitivities_L30_raw{k};
            Sensitivities_L30(k).ID          = d.ID_label;
            Sensitivities_L30(k).dPi_domega  = d.dPi_domega_instances;
            Sensitivities_L30(k).dPi_dthetaK = d.dPi_dthetaK_instances;
            Sensitivities_L30(k).dPi_dthetaS = d.dPi_dthetaS_instances;
            Sensitivities_L30(k).dPi_dN      = d.dPi_dN_instances;
        end
    case 'wells'
        for k = 1:nTUs_L30
            d = Sensitivities_L30_raw{k};
            Sensitivities_L30(k).ID          = d.ID_label;
            Sensitivities_L30(k).dPi_domega  = d.dPi_domega_wells;
            Sensitivities_L30(k).dPi_dthetaK = d.dPi_dthetaK_wells;
            Sensitivities_L30(k).dPi_dthetaS = d.dPi_dthetaS_wells;
            Sensitivities_L30(k).dPi_dN      = d.dPi_dN_wells;
        end
    otherwise
        error('Invalid DATA_LAYER: %s. Use "Global", "Instances", or "Wells".', DATA_LAYER);
end

%% 7) Build and analyze Jacobian for Lib30
[incProm_L30, incRBS_L30, incORIs_L30] = get_free_params_for_library('L30');

[J_L30_exp, params_L30] = build_jacobian( ...
    IDs_L30, TUs_paper, ORIs, Promoters, RBSs, ...
    'IncludePromoters', incProm_L30, ...
    'IncludeRBS',       incRBS_L30, ...
    'IncludeORIs',      incORIs_L30, ...
    'Sensitivities',    Sensitivities_L30);

print_rank_summary(J_L30_exp, params_L30, sprintf('L30 (%s)', DATA_LAYER));

%% 8) Build and analyze Jacobian for Lib24 (subset of Lib30)
TUs_paper_L24 = { ...
    'L24_13','L24_14','L24_15','L24_16', ...
    'L24_17','L24_18','L24_19','L24_20', ...
    'L24_21','L24_22','L24_23','L24_24', ...
    'L24_01','L24_02','L24_03','L24_04', ...
    'L24_05','L24_06','L24_07','L24_08', ...
    'L24_09','L24_10','L24_11','L24_12'};
IDs_L24 = TUs_paper_L24;

[incProm_L24, incRBS_L24, incORIs_L24] = get_free_params_for_library('L24');

[J_L24_exp, params_L24] = build_jacobian( ...
    IDs_L24, TUs_paper_L24, ORIs, Promoters, RBSs, ...
    'IncludePromoters', incProm_L24, ...
    'IncludeRBS',       incRBS_L24, ...
    'IncludeORIs',      incORIs_L24, ...
    'Sensitivities',    Sensitivities_L30);

print_rank_summary(J_L24_exp, params_L24, sprintf('L24 (%s)', DATA_LAYER));

%% 9) Build and analyze Jacobian for Lib6 (subset of Lib30)
TUs_paper_L6 = {'L30_28','L30_29','L30_30','L30_25','L30_26','L30_27'};
IDs_L6 = TUs_paper_L6;

[incProm_L6, incRBS_L6, incORIs_L6] = get_free_params_for_library('L6');

[J_L6_exp, params_L6] = build_jacobian( ...
    IDs_L6, TUs_paper_L6, ORIs, Promoters, RBSs, ...
    'IncludePromoters', incProm_L6, ...
    'IncludeRBS',       incRBS_L6, ...
    'IncludeORIs',      incORIs_L6, ...
    'Sensitivities',    Sensitivities_L30);

print_rank_summary(J_L6_exp, params_L6, sprintf('L6 (%s)', DATA_LAYER));

%% ----------------------------- helpers -----------------------------
function print_rank_summary(J, params, label)
    % Exact rank and effective rank based on singular values
    [~,S,~] = svd(J,'econ');
    s = diag(S);
    r_exact = rank(J);
    if isempty(s)
        tol = NaN;
        r_eff = 0;
    else
        tol = 1e-2 * s(1);
        r_eff = sum(s > tol);
    end

    fprintf('--- %s ---\n', label);
    fprintf('  #params        : %d\n', numel(params));
    fprintf('  rank(J)        : %d\n', r_exact);
    fprintf('  effective rank : %d (tol = %.2e)\n\n', r_eff, tol);
end
