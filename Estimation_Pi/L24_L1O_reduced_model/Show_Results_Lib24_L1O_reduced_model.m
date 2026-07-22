%% Show_Results_Lib24_L1O_reduced_model
% Visualize the distributed Lib24 leave-one-out reduced-model results.
%
% Generates publication-oriented plots of:
%   - effective pGreen copy number;
%   - promoter transcription rates;
%   - RBS intrinsic initiation capacities;
%   - experimental and predicted synthesis-rate trajectories.
%
% The script calls Plot_Results_Lib_reduced_model using its backward-compatible
% seven-argument interface followed by name-value options.
%
% USAGE
%   Show_Results_Lib24_L1O_reduced_model
%
% See README.md for the complete workflow.

% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin(); %#ok<NASGU>
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);

% Model constants retained for compatibility with the original workflow.
model_c.lp_c = 240;
model_c.le_c = model_c.lp_c^0.097/0.0703;
model_c.dm_c  =  0.2;
model_c.Em_c =  model_c.lp_c/model_c.le_c*(1- (model_c.lp_c/(model_c.lp_c+model_c.le_c))^(model_c.lp_c/model_c.le_c));
model_c.WEm_c =  1 + 1/model_c.Em_c;
model_c.N_pSC101 = 5;

Lower_Bounds = [0.025,0.05,0.05, ...
          0.05, 2.5e-3,  2.5e-3,2.5e-3, ...
          2.5];

Upper_Bounds = [0.25,0.45,0.45, ...
          2.0, 20e-3,  20e-3,5e-2, ...
          8.25];

% Lib24 bioparts indices in the canonical color/order convention:
% RBS canonical order is [B0030, B0032, B0034, J61100, J61101],
% so Lib24 uses [1,2,4,5].
indices_plasmids  = [1,2];
indices_promoters = [1,2,3];
indices_rbss      = [1,2,4,5];

% Select data level.
instances = false;
wells = true;

OUTPUT_DIR = fullfile(SCRIPT_DIR,'Figures');
if ~exist(OUTPUT_DIR,'dir')
    mkdir(OUTPUT_DIR);
end

if instances
    results_file = fullfile(SCRIPT_DIR,'Results_Tensor_Lib24_L1O_reduced_Instances.mat');
    S = load(results_file);
    if ~isfield(S,'Results_Tensor_Lib24_L1O_reduced')
        error('Show_Results_Lib24_L1O_reduced_model:MissingVariable', ...
            'Expected variable Results_Tensor_Lib24_L1O_reduced in %s.', results_file);
    end
    Results_Tensor_Lib24_L1O_reduced = S.Results_Tensor_Lib24_L1O_reduced;
    title_text = '';
    Plot_Results_Lib_reduced_model(Results_Tensor_Lib24_L1O_reduced, ...
        Lower_Bounds,Upper_Bounds,indices_plasmids,indices_promoters,indices_rbss,title_text, ...
        'OutputDir',OUTPUT_DIR, ...
        'FigurePrefix','Lib24_L1O_reduced_Instances', ...
        'SaveFigures',true, ...
        'BaseFontSize',26, ...
        'TileFontSize',18, ...
        'AxisLabelFontSize',24);
end

if wells
    results_file = fullfile(SCRIPT_DIR,'Results_Tensor_Lib24_L1O_reduced_Wells.mat');
    S = load(results_file);
    if ~isfield(S,'Results_Tensor_Lib24_L1O_reduced')
        error('Show_Results_Lib24_L1O_reduced_model:MissingVariable', ...
            'Expected variable Results_Tensor_Lib24_L1O_reduced in %s.', results_file);
    end
    Results_Tensor_Lib24_L1O_reduced = S.Results_Tensor_Lib24_L1O_reduced;
    title_text = '';
    Plot_Results_Lib_reduced_model(Results_Tensor_Lib24_L1O_reduced, ...
        Lower_Bounds,Upper_Bounds,indices_plasmids,indices_promoters,indices_rbss,title_text, ...
        'OutputDir',OUTPUT_DIR, ...
        'FigurePrefix','Lib24_L1O_reduced_Wells', ...
        'SaveFigures',true, ...
        'BaseFontSize',26, ...
        'TileFontSize',18, ...
        'AxisLabelFontSize',24,...
        'FigCopySizeCm',[30 20],...
        'FigPromoterSizeCm',[65 30],...
        'FigRBSSizeCm',[110 35],...
        'FigTilesSizeCm',[150 75]);
end
