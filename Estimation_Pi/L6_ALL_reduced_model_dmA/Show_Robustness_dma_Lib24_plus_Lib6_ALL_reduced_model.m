%% Show_Robustness_dma_Lib24_plus_Lib6_ALL_reduced_model
% Combine the Lib24 and Lib6 d_mA robustness analyses.
%
% Shared promoter, plasmid, and four RBS parameters are read from the sibling
% L24_ALL_reduced_model_dmA workflow. B0034 is read from this Lib6 workflow.
%
% Both sets of Results tensors must be regenerated before running this script.
%
% Outputs:
%   Figures/
%   Tables/

clear; clc;

color_blue = hex2rgb('#00008B');
color_grey = hex2rgb('#4b4b4b'); %#ok<NASGU>
color_grey_boira = hex2rgb('#F2F2F2'); %#ok<NASGU>
color_grey_neutre = hex2rgb('#CCCCCC'); %#ok<NASGU>
RBS_colors_hex = {'#A8DADC','#F4A261','#B5E48C','#CDB4DB','#FFE066'}; % B0030, B0032, J61100, J61101, B0034
Promoter_colors_hex = {'#264653','#E76F51','#2A9D8F'}; % J23106, J23102, J23101

ROOT = init_SynTwin(); %#ok<NASGU>
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);
FIGURES_DIR = fullfile(SCRIPT_DIR,'Figures');
TABLES_DIR = fullfile(SCRIPT_DIR,'Tables');
if ~exist(FIGURES_DIR,'dir'), mkdir(FIGURES_DIR); end
if ~exist(TABLES_DIR,'dir'), mkdir(TABLES_DIR); end

% =========================
% Model / library setup
% =========================
model_c.lp_c = 240; % GFP protein length (aa)
model_c.le_c = model_c.lp_c^0.097/0.0703;
model_c.Em_c = model_c.lp_c/model_c.le_c * ...
    (1 - (model_c.lp_c/(model_c.lp_c+model_c.le_c))^(model_c.lp_c/model_c.le_c));
model_c.WEm_c = 1 + 1/model_c.Em_c;
model_c.N_pSC101 = 5; % known

Use_mean = 'Wells';
set_dma = [0.15, 0.17, 0.20, 0.22, 0.25];
dm_ref = 0.20;

dma_values = set_dma(:);
ndma = numel(dma_values);

% =========================
% Parameter indices / names / groups
% Combined view: Lib24 for shared parameters + Lib6 for B0034
% =========================
idx.J23106 = 1;
idx.J23102 = 2;
idx.J23101 = 3;
idx.B0030  = 4;
idx.B0032  = 5;
idx.J61100 = 6;
idx.J61101 = 7;
idx.pGreen = 8;
idx.B0034  = 9;

param_names = {'J23106','J23102','J23101','B0030','B0032','J61100','J61101','pGreen','B0034'};
param_family = {'Promoter','Promoter','Promoter','RBS','RBS','RBS','RBS','Plasmid','RBS'};
param_source = {'Lib24','Lib24','Lib24','Lib24','Lib24','Lib24','Lib24','Lib24','Lib6'};
group_map = {1:3, [4:7,9], 8};
group_titles = {'Promoters','RBSs','pGreen'};
n_params = numel(param_names);

% =========================
% Colors by parameter (consistent with prior figures)
% =========================
param_color = cell(1,n_params);
for j = 1:3
    param_color{j} = hex2rgb(Promoter_colors_hex{j});
end
for j = 1:4
    param_color{3+j} = hex2rgb(RBS_colors_hex{j});
end
param_color{idx.pGreen} = color_blue;
param_color{idx.B0034} = hex2rgb(RBS_colors_hex{5});

% =========================
% Results directories
% =========================
DIR_LIB6 = SCRIPT_DIR;
DIR_LIB24 = resolve_existing_dir({fullfile(fileparts(SCRIPT_DIR),'L24_ALL_reduced_model_dmA')});

% =========================
% Load Lib24 and Lib6 results independently
% =========================
lib24 = load_library_results_lib24(DIR_LIB24, dma_values, Use_mean, idx.pGreen);
lib6  = load_library_results_lib6(DIR_LIB6,  dma_values, Use_mean);

% =========================
% Summaries by dmA for each library
% =========================
lib24 = summarize_library(lib24, dm_ref, dma_values);
lib6  = summarize_library(lib6,  dm_ref, dma_values);

[~, i_ref24] = min(abs(dma_values - dm_ref));
[~, i_ref6]  = min(abs(dma_values - dm_ref));

% =========================
% Combined parameter set for Figure 2/3
% Lib24 for shared params; Lib6 for B0034
% =========================
theta_by_case_combined = cell(ndma,1);
theta_ref_combined = zeros(1, n_params);

for i = 1:ndma
    theta24_i = lib24.theta_by_case{i};
    theta6_i  = lib6.theta_by_case{i};

    theta_by_case_combined{i} = cell(1,n_params);
    for j = 1:n_params
        if strcmp(param_source{j}, 'Lib24')
            theta_by_case_combined{i}{j} = theta24_i(:,j);
        else
            theta_by_case_combined{i}{j} = theta6_i(:,1);
        end
    end
end

for j = 1:n_params
    if strcmp(param_source{j}, 'Lib24')
        theta_ref_combined(j) = lib24.Theta_median(i_ref24,j);
    else
        theta_ref_combined(j) = lib6.Theta_median(i_ref6,1);
    end
end

% =========================
% Balanced collapsed sampling per source library
% =========================
rng(1, 'twister');

target_runs24 = min(lib24.runs_per_case);
target_runs6  = min(lib6.runs_per_case);

balanced24 = balanced_theta(lib24.theta_by_case, target_runs24);
balanced6  = balanced_theta(lib6.theta_by_case,  target_runs6);

collapsed_balanced = cell(1, n_params);
for j = 1:n_params
    vals_abs = [];
    if strcmp(param_source{j}, 'Lib24')
        for i = 1:ndma
            vals_abs = [vals_abs; balanced24{i}(:,j)]; %#ok<AGROW>
        end
    else
        for i = 1:ndma
            vals_abs = [vals_abs; balanced6{i}(:,1)]; %#ok<AGROW>
        end
    end
    collapsed_balanced{j} = vals_abs;
end

% =========================
% Export summaries
% =========================
write_library_tables(lib24, TABLES_DIR, 'Lib24', {'J23106','J23102','J23101','B0030','B0032','J61100','J61101','pGreen'}, ...
    {'Promoter','Promoter','Promoter','RBS','RBS','RBS','RBS','Plasmid'}, dma_values);
write_library_tables(lib6,  TABLES_DIR, 'Lib6',  {'B0034'}, {'RBS'}, dma_values);

collapsed_summary = table(param_names(:), param_family(:), param_source(:), theta_ref_combined(:), ...
    cellfun(@median, collapsed_balanced(:)), ...
    cellfun(@(x) prctile(x,25), collapsed_balanced(:)), ...
    cellfun(@(x) prctile(x,75), collapsed_balanced(:)), ...
    'VariableNames', {'Parameter','Family','SourceLibrary','RefMedian_dmA_0p2', ...
    'CollapsedMedian_balanced','CollapsedQ25_balanced','CollapsedQ75_balanced'});
writetable(collapsed_summary, fullfile(TABLES_DIR, 'Table_dmA_collapsed_grouped_balanced_combined.csv'));

% =========================
% Shared aesthetics
% =========================
baseFont = 13;
titleFont = 14;
labelFont = 13;
axisLineWidth = 1.0;
gridAlpha = 0.16;
figColor = 'w';
cmap_dma = parula(ndma);

% =========================
% FIGURE 1: Cost function values + sensitivities (Lib24 and Lib6 separate)
% =========================
fig1 = figure('Units','centimeters','Position',[2 2 30 21], 'Color', figColor);
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

plot_cost_panel(lib24, dma_values, baseFont, titleFont, labelFont, axisLineWidth, gridAlpha, 'Lib24');
plot_sens_panel(lib24, {'J23106','J23102','J23101','B0030','B0032','J61100','J61101','pGreen'}, ...
    baseFont, titleFont, labelFont, axisLineWidth, gridAlpha, 'Lib24');
plot_cost_panel(lib6, dma_values, baseFont, titleFont, labelFont, axisLineWidth, gridAlpha, 'Lib6');
plot_sens_panel(lib6, {'B0034'}, baseFont, titleFont, labelFont, axisLineWidth, gridAlpha, 'Lib6');

export_figure(fig1, fullfile(FIGURES_DIR, 'Fig_dmA_CostAndSensitivity_Lib24_plus_Lib6'));

% =========================
% FIGURE 2: Individual parameter distributions vs dmA
% Combined parameter view: Lib24 shared params + Lib6 B0034
% =========================
fig2 = figure('Units','centimeters','Position',[2 2 33 22], 'Color', figColor);
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');

for j = 1:n_params
    nexttile; hold on;
    ymin_local = inf;
    ymax_local = -inf;

    for i = 1:ndma
        data = theta_by_case_combined{i}{j};
        data = data(isfinite(data));
        simple_violin(data, i, cmap_dma(i,:), 0.33, true, 0.16);
        if ~isempty(data)
            ymin_local = min(ymin_local, min(data));
            ymax_local = max(ymax_local, max(data));
        end
    end

    plot([1, ndma], [theta_ref_combined(j), theta_ref_combined(j)], '--', ...
        'Color', lighten_color(param_color{j}, 0.15), 'LineWidth', 1.2);

    title(sprintf('%s (%s)', param_names{j}, param_source{j}), 'Interpreter', 'none', 'FontSize', titleFont);
    xticks(1:ndma);
    xticklabels(compose('%.2f', dma_values));
    xlabel('$d_{mA}$','Interpreter','latex','FontSize',labelFont);
    ylabel('estimate','FontSize',labelFont);
    xlim([0.45, ndma+0.55]);
    if isfinite(ymin_local) && ymin_local > 0
        pad = 0.08*(ymax_local - ymin_local + eps);
        ylim([max(0, ymin_local - pad), ymax_local + pad]);
    end
    set(gca,'FontSize',baseFont,'LineWidth',axisLineWidth,'Layer','top');
    grid on; box on;
    ax = gca; ax.GridAlpha = gridAlpha; ax.MinorGridAlpha = gridAlpha*0.9;
end

export_figure(fig2, fullfile(FIGURES_DIR, 'Fig_dmA_IndividualViolins_ByParameter_Lib24_plus_Lib6'));

% =========================
% FIGURE 3: Grouped collapsed violins (absolute values, log scale)
% =========================
fig3 = figure('Units','centimeters','Position',[2 2 29 10.5], 'Color', figColor);
tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

for g = 1:3
    idxs = group_map{g};
    nexttile; hold on;

    local_vals = [];
    for k = 1:numel(idxs)
        j = idxs(k);

        data_collapsed = collapsed_balanced{j};
        data_collapsed = data_collapsed(data_collapsed > 0 & isfinite(data_collapsed));

        if strcmp(param_source{j}, 'Lib24')
            data_ref = lib24.theta_by_case{i_ref24}(:,j);
        else
            data_ref = lib6.theta_by_case{i_ref6}(:,1);
        end
        data_ref = data_ref(data_ref > 0 & isfinite(data_ref));

        if strcmp(param_names{j}, 'pGreen')
            col_dark = color_blue;
            col_light = lighten_color(color_blue, 0.62);
        else
            col_dark = param_color{j};
            col_light = lighten_color(param_color{j}, 0.58);
        end

        simple_violin(data_collapsed, k, col_light, 0.34, true, 0.08);
        simple_violin(data_ref,       k, col_dark,  0.18, true, 0.18);

        plot(k, median(data_ref), 'kd', ...
            'MarkerFaceColor', 'k', 'MarkerSize', 4.8, 'LineWidth', 1.0);

        local_vals = [local_vals; data_collapsed(:); data_ref(:)]; %#ok<AGROW>
    end

    if g == 2
        xline(4.5, '-', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.1);
    end

    set(gca, 'YScale', 'log');
    if ~isempty(local_vals)
        ymin = min(local_vals);
        ymax = max(local_vals);
        if ymin > 0 && ymax > ymin
            ylim([10^(log10(ymin)-0.15), 10^(log10(ymax)+0.15)]);
        end
    end

    xticks(1:numel(idxs));
    xticklabels(param_names(idxs));
    if numel(idxs) > 1
        xtickangle(25);
    end

    ylabel('estimate','FontSize',labelFont);
    title(group_titles{g}, 'Interpreter','none', 'FontSize', titleFont);
    xlim([0.45, numel(idxs)+0.55]);
    set(gca,'FontSize',baseFont,'LineWidth',axisLineWidth,'Layer','top');
    grid on; box on;
    ax = gca; ax.GridAlpha = gridAlpha; ax.MinorGridAlpha = gridAlpha*0.9;
end

export_figure(fig3, fullfile(FIGURES_DIR, 'Fig_dmA_CollapsedViolins_Grouped_Balanced_AbsoluteLog_OverlayRef_Lib24_plus_Lib6'));

fprintf('\n=== dmA robustness summary (combined Lib24 + Lib6) ===\n');
fprintf('Reference dmA = %.3f\n', dm_ref);
fprintf('Balanced collapsed analysis uses %d runs per dmA for Lib24 and %d runs per dmA for Lib6.\n', target_runs24, target_runs6);
disp(collapsed_summary);

% =========================
% Local helpers
% =========================
function lib = load_library_results_lib24(results_dir, dma_values, Use_mean, pgreen_idx)
    ndma = numel(dma_values);
    lib.tag = 'Lib24';
    lib.results_dir = results_dir;
    lib.theta_by_case = cell(ndma,1);
    lib.J_by_case = cell(ndma,1);
    lib.runs_per_case = zeros(ndma,1);

    for i = 1:ndma
        dm_val = dma_values(i);
        file_name = "Results_Tensor_Lib24_ALL_reduced_model_" + ...
                    Use_mean + "_dma_" + extractAfter(num2str(dm_val), '.') + ".mat";
        S = load(fullfile(results_dir, file_name));
        theta_i = S.Estimated_parameters.ALL_raw;
        theta_i(:,pgreen_idx) = 5 * theta_i(:,pgreen_idx);
        lib.theta_by_case{i} = theta_i;
        lib.J_by_case{i} = S.Estimated_parameters.J_raw(:);
        lib.runs_per_case(i) = size(theta_i,1);
    end
end

function lib = load_library_results_lib6(results_dir, dma_values, Use_mean)
    ndma = numel(dma_values);
    lib.tag = 'Lib6';
    lib.results_dir = results_dir;
    lib.theta_by_case = cell(ndma,1);
    lib.J_by_case = cell(ndma,1);
    lib.runs_per_case = zeros(ndma,1);

    for i = 1:ndma
        dm_val = dma_values(i);
        file_name = "Results_Tensor_Lib6_ALL_reduced_model_" + ...
                    Use_mean + "_dma_" + extractAfter(num2str(dm_val), '.') + ".mat";
        S = load(fullfile(results_dir, file_name));
        theta_i = S.Estimated_parameters.ALL_raw;
        if isvector(theta_i)
            theta_i = theta_i(:);
        end
        if size(theta_i,2) ~= 1 && size(theta_i,1) == 1
            theta_i = theta_i.';
        end
        if size(theta_i,2) ~= 1
            error('Lib6 loader expected a single estimated parameter (B0034), but got %d columns.', size(theta_i,2));
        end
        lib.theta_by_case{i} = theta_i;
        lib.J_by_case{i} = S.Estimated_parameters.J_raw(:);
        lib.runs_per_case(i) = size(theta_i,1);
    end
end

function lib = summarize_library(lib, dm_ref, dma_values)
    ndma = numel(lib.theta_by_case);
    n_params_here = size(lib.theta_by_case{1},2);
    lib.Theta_median = zeros(ndma, n_params_here);
    lib.Theta_q25 = zeros(ndma, n_params_here);
    lib.Theta_q75 = zeros(ndma, n_params_here);
    lib.J_median = zeros(ndma,1);
    lib.J_q25 = zeros(ndma,1);
    lib.J_q75 = zeros(ndma,1);

    for i = 1:ndma
        lib.Theta_median(i,:) = median(lib.theta_by_case{i},1);
        lib.Theta_q25(i,:) = prctile(lib.theta_by_case{i},25,1);
        lib.Theta_q75(i,:) = prctile(lib.theta_by_case{i},75,1);
        lib.J_median(i) = median(lib.J_by_case{i});
        lib.J_q25(i) = prctile(lib.J_by_case{i},25);
        lib.J_q75(i) = prctile(lib.J_by_case{i},75);
    end

    [~, lib.i_ref] = min(abs(dma_values - dm_ref));
    lib.theta_ref = lib.Theta_median(lib.i_ref,:);
    lib.J_ref = lib.J_median(lib.i_ref);

    lib.J_rel_median = 100*(lib.J_median - lib.J_ref)/lib.J_ref;
    lib.J_rel_q25 = zeros(ndma,1);
    lib.J_rel_q75 = zeros(ndma,1);
    for i = 1:ndma
        J_rel_i = 100*(lib.J_by_case{i} - lib.J_ref)/lib.J_ref;
        lib.J_rel_q25(i) = prctile(J_rel_i,25);
        lib.J_rel_q75(i) = prctile(J_rel_i,75);
    end

    log_dm = log(dma_values(:));
    lib.Sens_loglog = nan(1, n_params_here);
    lib.Sens_loglog_R2 = nan(1, n_params_here);
    for j = 1:n_params_here
        y = log(lib.Theta_median(:,j));
        p = polyfit(log_dm, y, 1);
        yhat = polyval(p, log_dm);
        lib.Sens_loglog(j) = p(1);
        ss_res = sum((y - yhat).^2);
        ss_tot = sum((y - mean(y)).^2);
        lib.Sens_loglog_R2(j) = 1 - ss_res/ss_tot;
    end
end

function theta_balanced = balanced_theta(theta_by_case, target_runs)
    ndma = numel(theta_by_case);
    theta_balanced = cell(ndma,1);
    for i = 1:ndma
        n_i = size(theta_by_case{i},1);
        if n_i > target_runs
            idx_sel = randperm(n_i, target_runs);
        else
            idx_sel = 1:n_i;
        end
        theta_balanced{i} = theta_by_case{i}(idx_sel,:);
    end
end

function write_library_tables(lib, script_dir, prefix, param_names_here, param_family_here, dma_values)
    cost_summary = table(dma_values(:), lib.runs_per_case, lib.J_median, lib.J_q25, lib.J_q75, lib.J_rel_median, lib.J_rel_q25, lib.J_rel_q75, ...
        'VariableNames', {'dmA','Nruns','J_median','J_q25','J_q75','DeltaJ_median_percent','DeltaJ_q25_percent','DeltaJ_q75_percent'});

    sensitivity_summary = table(param_names_here(:), param_family_here(:), lib.Sens_loglog(:), lib.Sens_loglog_R2(:), ...
        'VariableNames', {'Parameter','Family','LogLogSensitivity','R2'});

    writetable(cost_summary, fullfile(script_dir, ['Table_dmA_cost_summary_' prefix '.csv']));
    writetable(sensitivity_summary, fullfile(script_dir, ['Table_dmA_loglog_sensitivity_' prefix '.csv']));
end

function plot_cost_panel(lib, dma_values, baseFont, titleFont, labelFont, axisLineWidth, gridAlpha, titleTag)
    nexttile; hold on;
    cmap = parula(numel(dma_values));
    for i = 1:numel(dma_values)
        plot([dma_values(i), dma_values(i)], [lib.J_rel_q25(i), lib.J_rel_q75(i)], 'k-', 'LineWidth', 1.8);
        plot(dma_values(i), lib.J_rel_median(i), 'o', ...
            'MarkerFaceColor', cmap(i,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 7.0);
    end
    plot(dma_values, lib.J_rel_median, 'k-', 'LineWidth', 1.0);
    yline(0,'k--','LineWidth',1.0);
    xlabel('$d_{mA}$','Interpreter','latex','FontSize',labelFont);
    ylabel('Cost change (%)','FontSize',labelFont);
    title(['Cost function - ' titleTag],'Interpreter','none','FontSize',titleFont);
    set(gca,'FontSize',baseFont,'LineWidth',axisLineWidth,'Layer','top');
    grid on; box on;
    ax = gca; ax.GridAlpha = gridAlpha; ax.MinorGridAlpha = gridAlpha*0.9;
end

function plot_sens_panel(lib, tick_labels, baseFont, titleFont, labelFont, axisLineWidth, gridAlpha, titleTag)
    nexttile; hold on;
    bar(lib.Sens_loglog, 'FaceColor', [0.78 0.84 0.92], 'EdgeColor', 'k', 'LineWidth', 0.8);
    yline(0,'k-','LineWidth',0.8);
    xticks(1:numel(lib.Sens_loglog));
    xticklabels(tick_labels(1:numel(lib.Sens_loglog)));
    xtickangle(45);
    ylabel('$S_\theta^{d_{mA}}$','Interpreter','latex','FontSize',labelFont);
    title(['Log-log sensitivities - ' titleTag],'Interpreter','none','FontSize',titleFont);
    set(gca,'FontSize',baseFont,'LineWidth',axisLineWidth,'Layer','top');
    grid on; box on;
    ax = gca; ax.GridAlpha = gridAlpha; ax.MinorGridAlpha = gridAlpha*0.9;
end

function simple_violin(data, xpos, faceColor, halfWidth, showMedian, pointAlpha)
    data = data(:);
    data = data(isfinite(data));
    if isempty(data)
        return;
    end

    n = numel(data);
    if n == 1 || max(data) == min(data)
        plot(xpos, data, 'o', 'Color', faceColor, 'MarkerFaceColor', faceColor, 'MarkerSize', 5);
        return;
    end

    try
        [f, yi] = ksdensity(data, 'NumPoints', 200);
    catch
        yi = linspace(min(data), max(data), 200);
        bw = 1.06 * std(data) * n^(-1/5);
        if ~isfinite(bw) || bw <= 0
            bw = (max(data)-min(data))/20 + eps;
        end
        f = zeros(size(yi));
        for ii = 1:n
            f = f + exp(-0.5*((yi-data(ii))/bw).^2);
        end
        f = f/(n*bw*sqrt(2*pi));
    end

    if max(f) > 0
        f = f / max(f) * halfWidth;
    end

    patch([xpos + f, fliplr(xpos - f)], [yi, fliplr(yi)], faceColor, ...
        'FaceAlpha', 0.45, 'EdgeColor', faceColor, 'LineWidth', 1.0);

    jitter = (rand(n,1)-0.5) * halfWidth * 0.50;
    scatter(xpos + jitter, data, 8, 'MarkerFaceColor', faceColor, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', pointAlpha);

    q = prctile(data, [25 50 75]);
    plot([xpos-halfWidth*0.22, xpos+halfWidth*0.22], [q(1) q(1)], 'k-', 'LineWidth', 1.0);
    plot([xpos-halfWidth*0.22, xpos+halfWidth*0.22], [q(3) q(3)], 'k-', 'LineWidth', 1.0);
    plot([xpos, xpos], [q(1), q(3)], 'k-', 'LineWidth', 1.0);

    if showMedian
        plot([xpos-halfWidth*0.28, xpos+halfWidth*0.28], [q(2) q(2)], 'k-', 'LineWidth', 1.5);
    end
end

function export_figure(fig_handle, file_base)
    set(fig_handle, 'PaperPositionMode', 'auto');
    drawnow;
    exportgraphics(fig_handle, [file_base '.png'], 'Resolution', 300);
    exportgraphics(fig_handle, [file_base '.svg'], 'ContentType', 'vector');
end

function c = hex2rgb(hex)
    hex = char(hex);
    hex = strrep(hex, '#', '');
    if numel(hex) ~= 6
        error('hex2rgb:InvalidInput', 'Expected a 6-character hex code.');
    end
    c = [hex2dec(hex(1:2)), hex2dec(hex(3:4)), hex2dec(hex(5:6))] / 255;
end

function c2 = lighten_color(c1, frac_to_white)
    frac_to_white = max(0, min(1, frac_to_white));
    c2 = c1 + (1 - c1) * frac_to_white;
end

function outdir = resolve_existing_dir(candidates)
    outdir = '';
    for i = 1:numel(candidates)
        if isfolder(candidates{i})
            outdir = candidates{i};
            return;
        end
    end
    error('Could not resolve results directory from candidate paths.');
end
