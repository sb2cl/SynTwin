%% Show_Results_Lib6_ALL_reduced_model
% Visualize the distributed Lib6 all-construct reduced-model results.
%
% Loads Results_Tensor_Lib6_ALL_reduced_Wells.mat and calls
% Plot_Results_Lib6_reduced_model to generate the B0034 parameter histogram
% and the six experimental-versus-predicted synthesis-rate panels.
%
% USAGE
%   Show_Results_Lib6_ALL_reduced_model
%
% See README.md for the complete workflow.

% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin();  % Scripts_base + HEM_Surrogate are added by default
% Absolute path to the folder where this script lives (model folder)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);

% Getting the Results structure corresponding to the estimated parameters
% of Lib6 with the reduced model 
% CHANGE BETWEEN Global, Instances and Wells appropriately
results_file = fullfile(SCRIPT_DIR,'Results_Tensor_Lib6_ALL_reduced_Wells.mat');
S = load(results_file);
if ~isfield(S,'Results_Tensor_Lib6_ALL_reduced')
    error('Show_Results_Lib6_ALL_reduced_model:MissingVariable', ...
        'Expected variable Results_Tensor_Lib6_ALL_reduced in %s.', results_file);
end
Results_Tensor_Lib6_ALL_reduced = S.Results_Tensor_Lib6_ALL_reduced;

OUTPUT_DIR = fullfile(SCRIPT_DIR,'Figures');
if ~exist(OUTPUT_DIR,'dir')
    mkdir(OUTPUT_DIR);
end


% Consider an exogenous circuit. We use a Transcriptional Unit (TU) expressing GFP:
model_c.lp_c = 240; %Length of GFP protein (aa)     
model_c.le_c = model_c.lp_c^0.097/0.0703; 	%Ribosome occupancy length (aa) %lea=la^0.097/0.0703
model_c.dm_c  =  0.2; %Mean degradation rate of non-ribosomal mRNA (1/min)
model_c.Em_c =  model_c.lp_c/model_c.le_c*(1- (model_c.lp_c/(model_c.lp_c+model_c.le_c))^(model_c.lp_c/model_c.le_c)) ;  
model_c.WEm_c =  1 + 1/model_c.Em_c;  
model_c.N_pSC101 = 5; %known

 indices_plasmids = [1,2];
 indices_promoters = [1,2,3];
 indices_rbss = [3];
 
% param bounds:
 Lower_Bounds = 0.05; %  lower expected bounds for RBS_k0_sigma0, B0034
 Upper_Bounds   =  0.4;  %upper expected bounds for RBS_k0_sigma0,  B0034

    
title_text = '';
Plot_Results_Lib6_reduced_model(Results_Tensor_Lib6_ALL_reduced, ...
    Lower_Bounds,Upper_Bounds,indices_plasmids,indices_promoters,indices_rbss,title_text, ...
    'OutputDir',OUTPUT_DIR, ...
    'FigurePrefix','Lib6_ALL_reduced_Wells', ...
    'SaveFigures',true, ...
    'BaseFontSize',26, ...
    'TileFontSize',18, ...
    'AxisLabelFontSize',24, ...
    'FigRBSSizeCm',[30 20], ...
    'FigTilesSizeCm',[76 48]);

