%% Show_Results_Lib30_ALL_full_model
% Visualize the distributed Lib30 complete-model Wells results.
%
% Loads Results_Tensor_Lib30_ALL_full_model_Wells.mat and calls
% Plot_Results_Lib_complete_model to generate parameter, correlation, and
% experimental-versus-predicted synthesis-rate figures.
%
% Usage:
%   Show_Results_Lib30_ALL_full_model
%
% See README.md for the complete workflow.

% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin();  % Scripts_base + HEM_Surrogate are added by default
% Absolute path to the folder where this script lives (model folder)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);

OUTPUT_DIR = fullfile(SCRIPT_DIR,'Figures');
if ~exist(OUTPUT_DIR,'dir')
    mkdir(OUTPUT_DIR);
end

% Getting the Results structure corresponding to the estimated parameters
% of Lib30 with the complete model using
% Generate_Results_Lib30_ALL_full_model.m and stored in this same folder
% CHANGE BETWEEN Global, Instances and Wells appropriately
results_file = fullfile(SCRIPT_DIR,'Results_Tensor_Lib30_ALL_full_model_Wells.mat');
S = load(results_file);
if ~isfield(S,'Results_Tensor_Lib30_ALL_complete')
    error('Show_Results_Lib30_ALL_full_model:MissingVariable', ...
        'Expected variable Results_Tensor_Lib30_ALL_complete in %s.', results_file);
end
Results_Tensor_Lib30_ALL_complete = S.Results_Tensor_Lib30_ALL_complete;


% Consider an exogenous circuit. We use a Transcriptional Unit (TU) expressing GFP:
model_c.lp_c = 240; %Length of GFP protein (aa)     
model_c.le_c = model_c.lp_c^0.097/0.0703; 	%Ribosome occupancy length (aa) %lea=la^0.097/0.0703
model_c.dm_c  =  0.2; %Mean degradation rate of non-ribosomal mRNA (1/min)
model_c.Em_c =  model_c.lp_c/model_c.le_c*(1- (model_c.lp_c/(model_c.lp_c+model_c.le_c))^(model_c.lp_c/model_c.le_c)) ;  
model_c.WEm_c =  1 + 1/model_c.Em_c;  
model_c.N_pSC101 = 5; %known

 indices_plasmids = [1,2];
 indices_promoters = [1,2,3];
 indices_rbss = [1,2,3,4,5];
 
% param bounds:
Lower_Bounds = [0.025,0.05,0.05... %lower expected bounds for omega, promoter strengths,
      0.05, 2.5e-3, 0.05, 2.5e-3,2.5e-3 ...  %lower expected bounds for RBS_k0_sigma0, 
      0.002, 0.002, 0.002, 0.002, 0.002 ...   % lower expected boundsfor RBS_inv_sigma0
      2.5 ];  %lower expected bound for high copy number multiplying term 

Upper_Bounds = [0.15,0.35,0.35... %upper expected bounds for omega, promoter strengths,
      2.0, 10e-3, 0.35, 8e-3,2e-2 ...  %upper expected bounds for RBS_k0_sigma0, 
      0.05, 0.05, 0.05, 0.05, 0.05 ...   % upper expected bounds for RBS_inv_sigma0
      6.0 ];  %upper expected bound for high copy number multiplying term 
 
title_text = '';

Plot_Results_Lib_complete_model(Results_Tensor_Lib30_ALL_complete,Lower_Bounds,Upper_Bounds,indices_plasmids,indices_promoters,indices_rbss,title_text, ...
    'OutputDir',OUTPUT_DIR, ...
    'FigurePrefix','Lib30_ALL_full_model_Wells', ...
    'SaveFigures',true, ...
    'BaseFontSize',26, ...
    'TileFontSize',18, ...
    'AxisLabelFontSize',24, ...
    'TitleFontSize',22, ...
    'FigCopySizeCm',[30 20], ...
    'FigPromoterSizeCm',[65 30], ...
    'FigRBS2SizeCm',[150 60], ...
    'FigRBS3SizeCm',[180 75], ...
    'FigCorrSizeCm',[85 35], ...
    'FigTilesSizeCm',[150 75]);



