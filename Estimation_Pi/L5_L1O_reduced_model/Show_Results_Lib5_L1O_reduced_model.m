%% Show_Results_Lib5_L1O_reduced_model
% Visualize the distributed Lib5 leave-one-out reduced-model results.
%
% Loads Results_Tensor_Lib5_L1O_reduced_Wells.mat and calls the shared
% Scripts_base function Plot_Results_Lib5_reduced_model.
%
% USAGE
%   Show_Results_Lib5_L1O_reduced_model
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
% of Lib6 with the reduced model 
% CHANGE BETWEEN Global, Instances and Wells appropriately
results_file = fullfile(SCRIPT_DIR,'Results_Tensor_Lib5_L1O_reduced_Wells.mat');
S = load(results_file);
if ~isfield(S,'Results_Tensor_Lib5_L1O_reduced')
    error('Show_Results_Lib5_L1O_reduced_model:MissingVariable', ...
        'Expected Results_Tensor_Lib5_L1O_reduced in %s.', results_file);
end
Results_Tensor_Lib5_L1O_reduced = S.Results_Tensor_Lib5_L1O_reduced;

% Consider an exogenous circuit. We use a Transcriptional Unit (TU) expressing GFP:
model_c.lp_c = 240; %Length of GFP protein (aa)     
model_c.le_c = model_c.lp_c^0.097/0.0703; 	%Ribosome occupancy length (aa) %lea=la^0.097/0.0703
model_c.dm_c  =  0.2; %Mean degradation rate of non-ribosomal mRNA (1/min)
model_c.Em_c =  model_c.lp_c/model_c.le_c*(1- (model_c.lp_c/(model_c.lp_c+model_c.le_c))^(model_c.lp_c/model_c.le_c)) ;  
model_c.WEm_c =  1 + 1/model_c.Em_c;  
model_c.N_pSC101 = 5; %known

indices_plasmids_lib5 = 1;
num_plasmids_Lib5 = length(indices_plasmids_lib5);
indices_promoters_lib5 = 4;
num_promoters_Lib5 = length(indices_promoters_lib5);
indices_rbss_lib5 = [1,2,3,4,5];
num_rbss_Lib5 = length(indices_rbss_lib5);

 Lower_Bounds = 0.02;   %lower expected bounds   Promoter Omega,     J23100          
 Upper_Bounds = 1.5;  %upper expected bounds  Promoter Omega,     J23100    

 title_text = '';

 Plot_Results_Lib5_reduced_model(Results_Tensor_Lib5_L1O_reduced, ...
        Lower_Bounds,Upper_Bounds,indices_plasmids_lib5,indices_promoters_lib5,indices_rbss_lib5,title_text, ...
        'OutputDir',OUTPUT_DIR, ...
        'FigurePrefix','Lib5_L1O_reduced_Wells', ...
        'SaveFigures',true, ...
        'BaseFontSize',26, ...
        'TileFontSize',18, ...
        'AxisLabelFontSize',24,...
        'FigPromoterSizeCm',[45 30],...
        'FigTilesSizeCm',[75 25]);

