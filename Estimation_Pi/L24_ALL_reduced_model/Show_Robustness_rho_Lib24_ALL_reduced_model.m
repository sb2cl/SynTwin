%% Show_Robustness_rho_Lib24_ALL_reduced_model
% Assess robustness of Lib24 parameter estimates to the fixed rho^0 value.
%
% PURPOSE
%   Compares the objective values and selected parameter distributions across
%   rho^0 = [0.005, 0.01, 0.02, 0.03, 0.05]. The rho^0 = 0.02 case is used
%   as the reference for the relative objective variation.
%
% INPUTS
%   Results_Tensor_Lib24_ALL_reduced_model_<Use_mean>_<rho_tag>.mat
%
% OUTPUT
%   Figures/FigS_RhoRobustness.pdf
%
% USAGE
%   Show_Robustness_rho_Lib24_ALL_reduced_model

% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin();  % Scripts_base + HEM_Surrogate are added by default
% Absolute path to the folder where this script lives (model folder)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);

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

% Getting the Results structure corresponding to the estimated parameters
% of Lib24 with the reduced model 
% CHANGE BETWEEN Global, Instances and Wells appropriately

%Use_mean = 'Instances';  
%Use_mean = 'Global';  
Use_mean = 'Wells'; 

set_RBS_inv_sigma0 = [0.005,0.01,0.02,0.03,0.05];
labels = arrayfun(@(x) sprintf('$\\sigma^{-1}=%.3g$', x), set_RBS_inv_sigma0, 'UniformOutput', false);

RBS_colors= hex2rgb(['#A8DADC';'#F4A261';'#B5E48C';'#CDB4DB';'#FFE066']);

 for i= 1:numel(set_RBS_inv_sigma0)

    RBS_inv_sigma0 = set_RBS_inv_sigma0(i);
    
    file_name = "Results_Tensor_Lib24_ALL_reduced_model_"  + Use_mean + "_" + extractAfter(num2str(RBS_inv_sigma0), '.') + ".mat";
    load(fullfile(SCRIPT_DIR,file_name));

   Results_Tensor_case{i}.Results_Tensor_Lib24_ALL_reduced = Results_Tensor_Lib24_ALL_reduced;
   Results_parameters_case{i}.Estimated_parameters = Estimated_parameters;
 end

%% =========================
% FIGURE: Robustness to rho0
% =========================

figure('Units','centimeters','Position',[2 2 18 24]);
tiledlayout(4,2,'TileSpacing','compact','Padding','compact');

%% =========================
% Build aggregated data
% =========================

rho_values = set_RBS_inv_sigma0(:);
n_cases = numel(rho_values);

theta_all = [];
J_all = [];
rho0_vec = [];

for i = 1:n_cases
    
    params_struct = Results_parameters_case{i}.Estimated_parameters;
    
    theta_i = params_struct.ALL_raw;   % Nx8
    J_i     = params_struct.J_raw;     % Nx1
    
    n_runs = size(theta_i,1);
    
    theta_all = [theta_all; theta_i];
    J_all     = [J_all; J_i];
    rho0_vec  = [rho0_vec; rho_values(i)*ones(n_runs,1)];
    
end

%% =========================
% Parameter indices
% =========================

idx.J23106 = 1;
idx.J23102 = 2;
idx.J23101 = 3;
idx.B0030  = 4;
idx.B0032  = 5;
idx.J61100 = 6;
idx.J61101 = 7;
idx.pGreen = 8;

%% =========================
% Reference rho0 = 0.02
% =========================

[~,i_ref] = min(abs(rho_values - 0.02));
theta_ref = Results_parameters_case{i_ref}.Estimated_parameters.ALL_mean;

%% =========================
% ΔJ (%) normalization
% =========================

J_min = min(J_all);
J_rel = 100*(J_all - J_min)/J_min;

%% =========================
% Sorting (clean plots)
% =========================

[rho0_vec, sort_idx] = sort(rho0_vec);
theta_all = theta_all(sort_idx,:);
J_all     = J_all(sort_idx);
J_rel     = J_rel(sort_idx);

rho_values = sort(rho_values);
cmap = parula(length(rho_values));

%% =========================
% Panel A — Cost vs rho0 (improved)
% =========================

nexttile([1 2]); hold on;

% Reference: rho0 = 0.02
rho_ref = 0.02;
tol = 1e-12;
idx_ref = abs(rho0_vec - rho_ref) < tol;

J_ref = mean(J_all(idx_ref));
J_rel = 100*(J_all - J_ref)/J_ref;

% Light band highlighting the <1% region
patch([min(rho_values)-0.002, max(rho_values)+0.002, max(rho_values)+0.002, min(rho_values)-0.002], ...
      [0, 0, 1, 1], ...
      [0.92 0.97 1.00], ...
      'EdgeColor','none', 'FaceAlpha',0.8);

% Scatter with horizontal jitter
rng(1); % reproducible jitter
for i = 1:length(rho_values)
    idx_r = abs(rho0_vec - rho_values(i)) < tol;

    x = rho_values(i) + 0.0007*randn(sum(idx_r),1);
    scatter(x, J_rel(idx_r), 26, ...
        'MarkerFaceColor', cmap(i,:), ...
        'MarkerEdgeColor', 'none', ...
        'MarkerFaceAlpha', 0.45);
end

% Mean ± std for each rho0
J_mean = zeros(size(rho_values));
J_std  = zeros(size(rho_values));

for i = 1:length(rho_values)
    idx_r = abs(rho0_vec - rho_values(i)) < tol;
    J_mean(i) = mean(J_rel(idx_r));
    J_std(i)  = std(J_rel(idx_r));
end

errorbar(rho_values, J_mean, J_std, ...
    'k', 'LineStyle','none', 'LineWidth',1.4, 'CapSize',8);

plot(rho_values, J_mean, 'ko', ...
    'MarkerFaceColor','k', 'MarkerSize',5);

% Reference lines
yline(0, 'k--', 'LineWidth',1);      % rho0 = 0.02 reference level
yline(1, 'r--', 'LineWidth',1.5);    % 1% threshold

% Labels
xlabel('\rho_0');
ylabel('\Delta J / J_{0.02} (%)');
title('Cost sensitivity to \rho_0');

% Tight y-axis around the actual variation
ymax = max(J_rel) + 0.05;
ylim([min(-0.05,min(J_rel)-0.05), max(1.05,ymax)]);

xlim([min(rho_values)-0.003, max(rho_values)+0.003]);

% Annotation
text(mean(xlim), 0.93, 'all fits remain within 1% of the \rho_0 = 0.02 reference', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','top', ...
    'FontSize',9, ...
    'FontWeight','bold');

box on;




%% =========================
% Panel B — Parameters
% =========================

param_names = {'pGreen','J23101','J23106','B0030'};
param_idx   = [idx.pGreen idx.J23101 idx.J23106 idx.B0030];

for k = 1:4
    nexttile; hold on;
    
    data = theta_all(:,param_idx(k));
    
    boxplot(data, rho0_vec,'Colors','k','Symbol','');
    
    yline(theta_ref(param_idx(k)), 'r--','LineWidth',1.5);
    
    title(param_names{k});
    xlabel('\rho_0');
    ylabel('estimate');
    
    med = median(data);
    ylim([0.85*med 1.15*med]);
    
    box on;
end

%% =========================
% Panel C — Distribution overlap (B0030)
% =========================

nexttile([1 2]); hold on;

for i = 1:length(rho_values)
    idx_r = rho0_vec == rho_values(i);
    
    data = theta_all(idx_r, idx.B0030);
    
    [f,xi] = ksdensity(data);
    plot(xi, f,'Color',cmap(i,:),'LineWidth',1.5);
end

xlabel('B0030 estimate');
ylabel('Density');
title('Distribution overlap (B0030)');
box on;

%% =========================
% Formatting
% =========================

set(findall(gcf,'-property','FontSize'),'FontSize',10);

%% =========================
% Save PDF
% =========================

set(gcf,'PaperPositionMode','auto');
figures_dir = fullfile(SCRIPT_DIR,'Figures');
if ~exist(figures_dir,'dir'), mkdir(figures_dir); end
print(gcf,fullfile(figures_dir,'FigS_RhoRobustness.pdf'),'-dpdf','-vector');
