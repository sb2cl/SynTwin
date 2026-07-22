%% Show_Results_pARKA21_background_model
% Plot the initial pARKA21 experimental trajectories, fitted model, parameter
% distributions, and residuals.
%
% INPUT
%   Results_pARKA21_background_Wells.mat
%
% OUTPUT
%   Figures/*.png
%   Figures/*.svg

clearvars;
close all;
dbstop if error
warning on

ROOT = init_SynTwin(); %#ok<NASGU>
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);

results_file = fullfile(SCRIPT_DIR,'Results_pARKA21_background_Wells.mat');
S = load(results_file,'Results_pARKA21_background');
if ~isfield(S,'Results_pARKA21_background')
    error('Show_Results_pARKA21_background_model:MissingResults', ...
        'Expected Results_pARKA21_background in %s.',results_file);
end
R = S.Results_pARKA21_background;

figures_dir = fullfile(SCRIPT_DIR,'Figures');
if ~exist(figures_dir,'dir')
    mkdir(figures_dir);
end

colors = turbo(max(3,numel(R.Instances)));
colors = colors(round(linspace(1,size(colors,1),numel(R.Instances))),:);

% Initial experimental trajectories.
fig1 = figure('Color','w','Units','centimeters','Position',[1 1 48 24]);
ax1 = axes(fig1); hold(ax1,'on');
for i = 1:numel(R.Instances)
    plot(ax1,R.Instances{i}.Mu,R.Instances{i}.Pi_exp, ...
        'LineWidth',1.5,'Color',colors(i,:));
end
grid(ax1,'on'); box(ax1,'on');
xlabel(ax1,'$\mu\;(\mathrm{min}^{-1})$','Interpreter','latex');
ylabel(ax1,'$\Pi_{bk}$','Interpreter','latex');
title(ax1,'pARKA21 experimental background trajectories');
export_pair(fig1,figures_dir,'pARKA21_background_experimental');

% Fitted model and experimental trajectories.
fig2 = figure('Color','w','Units','centimeters','Position',[1 1 48 24]);
ax2 = axes(fig2); hold(ax2,'on');
for i = 1:numel(R.Instances)
    plot(ax2,R.Instances{i}.Mu,R.Instances{i}.Pi_exp, ...
        'LineWidth',1.2,'Color',colors(i,:));
end
plot(ax2,R.Mu_grid,R.Pi_grid,'k-','LineWidth',3.2);
grid(ax2,'on'); box(ax2,'on');
xlabel(ax2,'$\mu\;(\mathrm{min}^{-1})$','Interpreter','latex');
ylabel(ax2,'$\Pi_{bk}$','Interpreter','latex');
title(ax2,'pARKA21 background-model fit');
export_pair(fig2,figures_dir,'pARKA21_background_fit');

% Parameter distributions.
fig3 = figure('Color','w','Units','centimeters','Position',[1 1 55 24]);
tl = tiledlayout(fig3,1,5,'Padding','compact','TileSpacing','compact');
for p = 1:5
    ax = nexttile(tl);
    histogram(ax,R.Parameters_raw(:,p),'BinMethod','auto');
    xline(ax,R.Parameters_mean_bestN(p),'k-','LineWidth',2);
    title(ax,R.ParameterNames{p},'Interpreter','none');
    grid(ax,'on'); box(ax,'on');
end
export_pair(fig3,figures_dir,'pARKA21_background_parameters');

% Residuals.
fig4 = figure('Color','w','Units','centimeters','Position',[1 1 48 24]);
ax4 = axes(fig4); hold(ax4,'on');
for i = 1:numel(R.Instances)
    plot(ax4,R.Instances{i}.Mu,R.Instances{i}.Residual, ...
        'LineWidth',1.3,'Color',colors(i,:));
end
yline(ax4,0,'k-');
grid(ax4,'on'); box(ax4,'on');
xlabel(ax4,'$\mu\;(\mathrm{min}^{-1})$','Interpreter','latex');
ylabel(ax4,'$\Pi_{bk}-\widehat{\Pi}_{bk}$','Interpreter','latex');
title(ax4,'pARKA21 background-model residuals');
export_pair(fig4,figures_dir,'pARKA21_background_residuals');

function export_pair(fig,outdir,name)
drawnow;
exportgraphics(fig,fullfile(outdir,[name '.png']), ...
    'Resolution',300,'BackgroundColor','white');
exportgraphics(fig,fullfile(outdir,[name '.svg']), ...
    'ContentType','vector','BackgroundColor','white');
end
