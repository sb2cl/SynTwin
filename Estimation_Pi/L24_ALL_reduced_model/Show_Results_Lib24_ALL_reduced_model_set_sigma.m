%% Show_Results_Lib24_ALL_reduced_model_set_sigma
% Inspect one Lib24 result tensor generated at a selected fixed rho^0.
%
% PURPOSE
%   Loads one rho^0-specific result tensor, calls
%   Plot_Results_Lib_reduced_model, and visualizes optimizer convergence
%   diagnostics for promoter, RBS, and plasmid-copy-number estimates.
%
% CONFIGURATION
%   - Use_mean: 'Global', 'Instances', or 'Wells'
%   - RBS_inv_sigma0: select one value from set_RBS_inv_sigma0
%
% INPUT
%   Results_Tensor_Lib24_ALL_reduced_model_<Use_mean>_<rho_tag>.mat
%
% USAGE
%   Show_Results_Lib24_ALL_reduced_model_set_sigma

% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin();  % Scripts_base + HEM_Surrogate are added by default
% Absolute path to the folder where this script lives (model folder)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);

% Getting the Results structure corresponding to the estimated parameters
% of Lib24 with the reduced model 
% CHANGE BETWEEN Global, Instances and Wells appropriately

%Use_mean = 'Instances';  
%Use_mean = 'Global';  
Use_mean = 'Wells'; 

set_RBS_inv_sigma0 = [0.005,0.01,0.02,0.03,0.05];

% CHANGE  appropriately
RBS_inv_sigma0 = set_RBS_inv_sigma0(5);

file_name = "Results_Tensor_Lib24_ALL_reduced_model_"  + Use_mean + "_" + extractAfter(num2str(RBS_inv_sigma0), '.') + ".mat";
load(fullfile(SCRIPT_DIR,file_name));

% Consider an exogenous circuit. We use a Transcriptional Unit (TU) expressing GFP:
model_c.lp_c = 240; %Length of GFP protein (aa)     
model_c.le_c = model_c.lp_c^0.097/0.0703; 	%Ribosome occupancy length (aa) %lea=la^0.097/0.0703
model_c.dm_c  =  0.2; %Mean degradation rate of non-ribosomal mRNA (1/min)
model_c.Em_c =  model_c.lp_c/model_c.le_c*(1- (model_c.lp_c/(model_c.lp_c+model_c.le_c))^(model_c.lp_c/model_c.le_c)) ;  
model_c.WEm_c =  1 + 1/model_c.Em_c;  
model_c.N_pSC101 = 5; %known

 indices_plasmids = [1,2];
 indices_promoters = [1,2,3];
 indices_rbss = [1,2,4,5];


 Lower_Bounds = [0.025,0.05,0.05... %lower expected bounds for omega, promoter strengths,
          0.05, 2.5e-3, 0.05, 2.5e-3,2.5e-3 ...  %lower expected bounds for RBS_k0_sigma0, 
          2.5 ];  %lower expected bound for high copy number multiplying term 

    Upper_Bounds = [0.15,0.35,0.35... %upper expected bounds for omega, promoter strengths,
          2.0, 10e-3, 0.35, 8e-3,2e-2 ...  %upper expected bounds for RBS_k0_sigma0, 
          8.25 ];  %upper expected bound for high copy number multiplying term 
 


title_text = 'Lib24. ALL. RM.';

Plot_Results_Lib_reduced_model(Results_Tensor_Lib24_ALL_reduced, ...
    Lower_Bounds,Upper_Bounds,indices_plasmids,indices_promoters,indices_rbss,title_text, ...
    'OutputDir',fullfile(SCRIPT_DIR,'Figures'), ...
    'FigurePrefix',sprintf('Lib24_ALL_reduced_rho_%s',extractAfter(num2str(RBS_inv_sigma0),'.')), ...
    'SaveFigures',true)

figure('Name','Estimate vs initial guess Promoters','Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
for i = 1:3
    nexttile
    scatter(Estimated_parameters.X0_raw(:,i),Estimated_parameters.ALL_raw(:,i),35,Estimated_parameters.J_raw,'filled');
    hold on
    mn = min([Estimated_parameters.X0_raw(:,i); Estimated_parameters.ALL_raw(:,i)]);
    mx = max([Estimated_parameters.X0_raw(:,i); Estimated_parameters.ALL_raw(:,i)]);
    plot([mn mx],[mn mx],'k--');
    xlabel('$x_0$','Interpreter','latex');
    ylabel('$\hat{\theta}$','Interpreter','latex');
    grid on;
end

figure('Name','Estimate vs Cost Promoter','Color','w');
tiledlayout(1,3,'Padding','compact','TileSpacing','compact');
for i = 1:3
    nexttile
    scatter(Estimated_parameters.ALL_raw(:,i),Estimated_parameters.J_raw,35,Estimated_parameters.J_raw,'filled');
    xlabel('$\hat{\theta}$','Interpreter','latex');
    ylabel('$J$','Interpreter','latex');
    grid on;
end

figure('Name','Estimate vs initial guess RBSs','Color','w');
tiledlayout(1,4,'Padding','compact','TileSpacing','compact');
for i = 1:4
    nexttile
    scatter(Estimated_parameters.X0_raw(:,3+i),Estimated_parameters.ALL_raw(:,3+i),35,Estimated_parameters.J_raw,'filled');
    hold on
    mn = min([Estimated_parameters.X0_raw(:,3+i); Estimated_parameters.ALL_raw(:,3+i)]);
    mx = max([Estimated_parameters.X0_raw(:,3+i); Estimated_parameters.ALL_raw(:,3+i)]);
    plot([mn mx],[mn mx],'k--');
    xlabel('$x_0$','Interpreter','latex');
    ylabel('$\hat{\theta}$','Interpreter','latex');
    grid on;
end

figure('Name','Estimate vs Cost RBSs','Color','w');
tiledlayout(1,4,'Padding','compact','TileSpacing','compact');
for i = 1:4
    nexttile
    scatter(Estimated_parameters.ALL_raw(:,3+i),Estimated_parameters.J_raw,35,Estimated_parameters.J_raw,'filled');
    xlabel('$\hat{\theta}$','Interpreter','latex');
    ylabel('$J$','Interpreter','latex');
    grid on;
end

 

figure('Name','Estimate vs initial guess Ori','Color','w');
tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
    nexttile
    scatter(5*Estimated_parameters.X0_raw(:,8),5*Estimated_parameters.ALL_raw(:,8),35,Estimated_parameters.J_raw,'filled');
    hold on
    mn = min(5*[Estimated_parameters.X0_raw(:,8); Estimated_parameters.ALL_raw(:,8)]);
    mx = max(5*[Estimated_parameters.X0_raw(:,8); Estimated_parameters.ALL_raw(:,8)]);
    plot([mn mx],[mn mx],'k--');
    xlabel('$x_0$','Interpreter','latex');
    ylabel('$\hat{\theta}$','Interpreter','latex');
    grid on;
    

    figure('Name','Estimate vs cost Ori','Color','w');
    tiledlayout(1,1,'Padding','compact','TileSpacing','compact');
    nexttile
    scatter(5*Estimated_parameters.ALL_raw(:,8),Estimated_parameters.J_raw,35,Estimated_parameters.J_raw,'filled');
    xlabel('$\hat{\theta}$','Interpreter','latex');
    ylabel('$J$','Interpreter','latex');
    grid on;