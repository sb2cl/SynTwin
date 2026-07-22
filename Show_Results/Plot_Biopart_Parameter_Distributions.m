%% Plot_Biopart_Parameter_Distributions.m
% Plot the estimated distributions of plasmid, promoter, and RBS parameters.
%
% Associated publication:
%   "Host-aware Identification of Intrinsic Gene Expression Biopart
%    Parameters using Combinatorial Libraries"
%
% Main-text association:
%   - Figure 3a: effective plasmid copy-number distribution.
%   - Figure 3b: promoter transcription-rate distributions.
%   - Figure 3c: RBS intrinsic initiation-capacity distributions.
%
% The script does not rerun parameter estimation. It loads the complete
% stored estimation-result tensors and extracts the raw parameter samples
% used to construct the violin plots. Each violin is overlaid with the
% median, interquartile range, and 5th--95th percentile interval.
%
% Inputs:
%   Estimation_Pi/L24_L1O_reduced_model/
%       Results_Tensor_Lib24_L1O_reduced_Wells.mat
%   Estimation_Pi/L6_L1O_reduced_model/
%       Results_Tensor_Lib6_L1O_reduced_Wells.mat
%   Estimation_Pi/L5_ALL_reduced_model/
%       Results_Tensor_Lib5_ALL_reduced_Wells.mat
%
% Outputs in ./Figures:
%   Figure_3a_plasmid_copy_number.png/.svg/.pdf
%   Figure_3b_promoter_strengths.png/.svg/.pdf
%   Figure_3c_RBS_IIC.png/.svg/.pdf
%
% Requirements:
%   - MATLAB with violinplot and exportgraphics support.
%   - Statistics and Machine Learning Toolbox for prctile.
%
% Usage:
%   Run Plot_Biopart_Parameter_Distributions from any working directory
%   while the script remains inside the SynTwin repository.
%
% Data provenance:
%   The MAT files contain complete outputs from the corresponding
%   parameter-estimation workflows, including construct metadata, raw and
%   summary parameter estimates, sensitivities, synthesis predictions, and
%   Monte Carlo samples. They are not reduced plotting tables.
%
% See README_Show_Results.md for the tensor roles and panel mapping.

clearvars;
close all;
dbstop if error;

% Use a consistent sans-serif font for all non-mathematical text.
% Mathematical symbols rendered with the LaTeX interpreter retain the
% appropriate mathematical typography.
set(groot, ...
    'defaultAxesFontName','Arial', ...
    'defaultTextFontName','Arial', ...
    'defaultLegendFontName','Arial');

% Ensure this script folder is on the path (for local functions in this folder)
ROOT = init_SynTwin();
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);
FIGURES_ROOT = fileparts(SCRIPT_DIR);


%% ---- Colors ----
Ori_colors_hex      = '#F2F2F2';
RBS_colors_hex      = ['#A8DADC';'#F4A261';'#B5E48C';'#CDB4DB';'#FFE066'];
Promoter_colors_hex = ['#264653';'#E76F51';'#2A9D8F';'#6A0572'];

RBS_colors      = hex2rgb_local(RBS_colors_hex);
Promoter_colors = hex2rgb_local(Promoter_colors_hex);
Ori_colors      = hex2rgb_local(Ori_colors_hex);

%% ---- Output settings ----
save_figures = true;
fig_output_dir = fullfile(SCRIPT_DIR, 'Figures');
if save_figures && ~exist(fig_output_dir,'dir')
    mkdir(fig_output_dir);
end

% Publication aesthetics matched to ETR_kappa_surface_fixed_v8.m.
% These figures are intentionally created large and then scaled down in the manuscript.
font_axis  = 1.5*28;
font_label = 1.25*34;
font_title = 1.25*28;
line_width_axis = 1.25*1.2;
line_width_edge = 1.25*1.2;
line_width_p95  = 2*1.5;
line_width_iqr  = 2*4.2;
marker_median   = 1.25*12;
export_resolution = 600;
export_pause = 3.2;

% Violin geometry.
% widthFactor scales the native violin width; keep slightly narrow for compact panels.
widthFactor = 0.62;
face_alpha  = 0.78;

% Canvas sizes in cm. Width is adapted to the number of violins.
% For comparison: ETR script uses FigProjectionSizeCm = [64 38].
height_violin_cm = 34;
width_per_violin_cm = 10.5;
base_width_cm = 12.0;
min_width_cm = 24.0;
max_width_cm = 1.25*64.0;

%% ---- Load final estimation results ----
required_files = { ...
    'Results_Tensor_Lib24_L1O_reduced_Wells.mat', ...
    'Results_Tensor_Lib6_L1O_reduced_Wells.mat', ...
    'Results_Tensor_Lib5_ALL_reduced_Wells.mat'};


% Getting the data of Lib24 L1O reduced model: 
S24 = load(SynTwin_path('Estimation_Pi/L24_L1O_reduced_model','Results_Tensor_Lib24_L1O_reduced_Wells.mat'), ...
    'Results_Tensor_Lib24_L1O_reduced');
% Getting the estimated IIC (K/sigma) of RBS B0034 from Lib6: 
S6 = load(SynTwin_path('Estimation_Pi/L6_L1O_reduced_model','Results_Tensor_Lib6_L1O_reduced_Wells.mat'), ...
    'Results_Tensor_Lib6_L1O_reduced');
% Getting the estimated omega of promoter J23100 from Lib5: 
S5 = load(SynTwin_path('Estimation_Pi/L5_ALL_reduced_model','Results_Tensor_Lib5_ALL_reduced_Wells.mat'), ...
    'Results_Tensor_Lib5_ALL_reduced');

assert(isfield(S24,'Results_Tensor_Lib24_L1O_reduced'), ...
    'Variable Results_Tensor_Lib24_L1O_reduced is missing from %s.', required_files{1});
assert(isfield(S6,'Results_Tensor_Lib6_L1O_reduced'), ...
    'Variable Results_Tensor_Lib6_L1O_reduced is missing from %s.', required_files{2});
assert(isfield(S5,'Results_Tensor_Lib5_ALL_reduced'), ...
    'Variable Results_Tensor_Lib5_ALL_reduced is missing from %s.', required_files{3});

Results_Tensor_Lib24_L1O_reduced = S24.Results_Tensor_Lib24_L1O_reduced;
Results_Tensor_Lib6_L1O_reduced = S6.Results_Tensor_Lib6_L1O_reduced;
Results_Tensor_Lib5_ALL_reduced = S5.Results_Tensor_Lib5_ALL_reduced;

%% ---- Extract final distributions ----
RBS_names      = ["B0030","B0032","B0034","J61100","J61101"];
Promoter_names = ["J23106","J23102","J23101","J23100"];
Ori_names      = ["pGreen"];

Lib24_rbs_indices      = [1,2,NaN,3,4];
Lib24_promoter_indices = [1,2,3];

Ori_values = {};
Ori_values{1}.name_Ori = Ori_names(1);
Ori_values{1}.Gene_cn_raw  = Results_Tensor_Lib24_L1O_reduced{1,1,1}.Gene_cn_raw;
Ori_values{1}.Gene_cn_mean = Results_Tensor_Lib24_L1O_reduced{1,1,1}.Gene_cn_mean;
Ori_values{1}.Gene_cn_std  = Results_Tensor_Lib24_L1O_reduced{1,1,1}.Gene_cn_std;

RBS_values = {};
for i = 1:5
    RBS_values{i}.name_RBS = RBS_names(i);
    if i == 3
        src = Results_Tensor_Lib6_L1O_reduced{1,1};
    else
        src = Results_Tensor_Lib24_L1O_reduced{1,1,Lib24_rbs_indices(i)};
    end
    RBS_values{i}.RBS_k0_sigma0_raw  = src.RBS_k0_sigma0_raw;
    RBS_values{i}.RBS_k0_sigma0_mean = src.RBS_k0_sigma0_mean;
    RBS_values{i}.RBS_k0_sigma0_std  = src.RBS_k0_sigma0_std;
end

Promoter_values = {};
for i = 1:4
    Promoter_values{i}.name_Promoter = Promoter_names(i);
    if i == 4
        src = Results_Tensor_Lib5_ALL_reduced{1,1};
    else
        src = Results_Tensor_Lib24_L1O_reduced{1,i,1};
    end
    Promoter_values{i}.Omega_raw  = src.Omega_raw;
    Promoter_values{i}.Omega_mean = src.Omega_mean;
    Promoter_values{i}.Omega_std  = src.Omega_std;
end

%% ---- RBS violin plot: estimated kappa ----
rbs_data = cell(1,numel(RBS_names));
for i = 1:numel(RBS_names)
    rbs_data{i} = RBS_values{i}.RBS_k0_sigma0_raw(:);
end

plot_violin_publication( ...
    rbs_data, RBS_names, RBS_colors, ...
    'RBS IIC', '$\kappa^0$', ...
    fullfile(fig_output_dir,'Figure_3c_RBS_IIC'), ...
    [0.005 0.01 0.1 0.2 0.5 1], ...
    'log', ...
    save_figures, export_resolution, export_pause, ...
    font_axis, font_label, font_title, line_width_axis, line_width_edge, ...
    line_width_p95, line_width_iqr, marker_median, widthFactor, face_alpha, ...
    height_violin_cm, 1.2*width_per_violin_cm, 1.2*base_width_cm, min_width_cm, max_width_cm);

%% ---- Promoter violin plot: estimated omega ----
promoter_data = cell(1,numel(Promoter_names));
for i = 1:numel(Promoter_names)
    promoter_data{i} = Promoter_values{i}.Omega_raw(:);
end

plot_violin_publication( ...
    promoter_data, Promoter_names, Promoter_colors, ...
    'Promoter strengths', '$\omega$', ...
    fullfile(fig_output_dir,'Figure_3b_promoter_strengths'), ...
    [0.05 0.1 0.15 0.2], ...
    'log', ...
    save_figures, export_resolution, export_pause, ...
    font_axis, font_label, 1.2*font_title, line_width_axis, line_width_edge, ...
    line_width_p95, line_width_iqr, marker_median, widthFactor, face_alpha, ...
    height_violin_cm, width_per_violin_cm, base_width_cm, min_width_cm, max_width_cm);

%% ---- Plasmid copy-number violin plot ----
ori_data = cell(1,numel(Ori_names));
for i = 1:numel(Ori_names)
    ori_data{i} = Ori_values{i}.Gene_cn_raw(:);
end

plot_violin_publication( ...
    ori_data, Ori_names, Ori_colors, ...
    'Plasmid copy number', '$N$', ...
    fullfile(fig_output_dir,'Figure_3a_plasmid_copy_number'), ...
    [], ...
    'linear_auto3', ...
    save_figures, export_resolution, export_pause, ...
    font_axis, font_label, font_title, line_width_axis, line_width_edge, ...
    line_width_p95, line_width_iqr, marker_median, 0.5*widthFactor, face_alpha, ...
    height_violin_cm, 0.5*width_per_violin_cm, 0.5*base_width_cm, 0.5*min_width_cm, max_width_cm);

%% ========================================================================
%% Local functions
%% ========================================================================

function fig = plot_violin_publication(dataCell, groupNames, groupColors, ...
    titleText, yLabelText, basename, preferredTicks, yScaleMode, ...
    save_figures, export_resolution, export_pause, ...
    font_axis, font_label, font_title, line_width_axis, line_width_edge, ...
    line_width_p95, line_width_iqr, marker_median, widthFactor, face_alpha, ...
    height_cm, width_per_violin_cm, base_width_cm, min_width_cm, max_width_cm)

    nGroups = numel(groupNames);

    % Width adapted to the number of violins.
    % Single-violin panels still need enough width for large title/tick text
    % after manuscript scaling.
    if nGroups == 1
        % Single-violin panel: deliberately wider so the title, tick labels,
        % and violin do not feel boxed in after manuscript scaling.
        fig_width_cm = 22.5;
    else
        fig_width_cm = base_width_cm + width_per_violin_cm * nGroups;
        fig_width_cm = max(min_width_cm, min(max_width_cm, fig_width_cm));
    end
    fig_size_cm = [fig_width_cm height_cm];

    % Build numeric matrix, one column per group.
    [Ymat, cleanData] = build_y_matrix(dataCell);

    % Use numeric x positions rather than categorical x-data.
    % This avoids categorical-axis limitations in xlim while keeping
    % custom group labels through XTickLabel.
    xPos = 1:nGroups;

    fig = figure('Color','w','Name',char(regexprep(titleText,'[\$\{\}\\]','')), ...
        'Units','centimeters','Position',[1 1 fig_size_cm]);

    ax = axes(fig);
    hold(ax,'on');

    % Native MATLAB violinplot. This version returns matlab.graphics.chart.primitive.ViolinPlot.
    % Use numeric positions so that xlim and percentile overlays are robust.
    vp = violinplot(xPos, Ymat);

    % Hide raw data points if generated.
    sc = findobj(ax,'Type','Scatter');
    if ~isempty(sc)
        set(sc,'Visible','off');
    end

    % Style internal patch objects.
    for i = 1:numel(vp)
        p = findobj(vp(i),'Type','Patch');
        if ~isempty(p)
            for k = 1:numel(p)
                p(k).FaceColor = groupColors(i,:);
                p(k).FaceAlpha = face_alpha;
                p(k).EdgeColor = darken_color(groupColors(i,:),0.35);
                p(k).LineWidth = line_width_edge;
            end
        end

        % De-emphasize any internal lines created by violinplot.
        l = findobj(vp(i),'Type','Line');
        for k = 1:numel(l)
            l(k).LineWidth = 1.0;
            l(k).Color = 0.15*[1 1 1];
        end
    end

    shrink_violin_patches(ax, widthFactor);

    % Add median + percentile overlay:
    % thin 5-95 interval, thick IQR, black median dot.
    for i = 1:nGroups
        yi = cleanData{i};
        if isempty(yi), continue; end

        q5  = prctile(yi, 5);
        q25 = prctile(yi,25);
        q50 = prctile(yi,50);
        q75 = prctile(yi,75);
        q95 = prctile(yi,95);

        line(ax, [i i], [q5 q95], ...
            'Color', 0.18*[1 1 1], 'LineWidth', line_width_p95, ...
            'Clipping','on');

        line(ax, [i i], [q25 q75], ...
            'Color', 0.05*[1 1 1], 'LineWidth', line_width_iqr, ...
            'Clipping','on');

        plot(ax, i, q50, 'o', ...
            'MarkerSize', marker_median, ...
            'MarkerFaceColor', 0.05*[1 1 1], ...
            'MarkerEdgeColor', 'w', ...
            'LineWidth', 1.1, ...
            'Clipping','on');
    end

    % Y-axis scale and ticks.
    apply_y_axis_ticks(ax, Ymat, preferredTicks, yScaleMode);

    % X-axis limits and labels. Numeric x positions avoid categorical-axis errors.
    if nGroups == 1
        % Give the single-violin panel extra horizontal air; otherwise the
        % violin and large labels feel too close to the axes box.
        xlim(ax,[0.25 1.75]);
    else
        xlim(ax,[0.45 nGroups+0.55]);
    end
    ax.XTick = xPos;
    ax.XTickLabel = cellstr(groupNames);

    % Publication formatting, applied last because violinplot can reset axes.
    ax.Box = 'on';
    ax.LineWidth = line_width_axis;
    ax.TickDir = 'out';
    ax.FontUnits = 'points';
    ax.FontName = 'Arial';
    ax.FontSize = font_axis;
    ax.TickLabelInterpreter = 'latex';
    ax.XLabel.Interpreter = 'latex';
    ax.YLabel.Interpreter = 'latex';
    ax.Title.Interpreter = 'latex';

    title(ax,titleText, ...
        'Interpreter','latex', ...
        'FontName','Arial', ...
        'FontSize',font_title, ...
        'FontWeight','bold');
    ylabel(ax,yLabelText, ...
        'Interpreter','latex', ...
        'FontName','Arial', ...
        'FontSize',font_label);

    % No x-label: group labels are self-explanatory and this saves vertical space.
    xlabel(ax,'');

    % Grid: light and unobtrusive.
    grid(ax,'on');
    % Light grid. For multi-violin panels, keep a subtle x-grid to help
    % align each violin with its label without making the plot busy.
    if nGroups == 1
        ax.XGrid = 'off';
    else
        ax.XGrid = 'on';
    end
    ax.YGrid = 'on';
    ax.GridAlpha = 0.11;
    ax.MinorGridAlpha = 0.05;
    ax.YMinorGrid = 'on';

    % Improve layout and margins for large tick labels and LaTeX text.
    % Leave extra top margin for large titles and extra horizontal room in
    % single-violin panels.
    ax.Units = 'normalized';
    if nGroups == 1
        % More external margin and less filled plotting area for the plasmid
        % copy-number panel: visually less boxed-in, with room for title/ticks.
        ax.Position = [0.22 0.21 0.62 0.60];
    else
        ax.Position = [0.14 0.18 0.82 0.68];
    end

    % Ensure all text objects inherit final sizes if created after axes formatting.
    set(findall(fig,'Type','Text'), 'FontUnits','points');

    drawnow;

    if save_figures
        save_figure_publication(fig, basename, export_resolution, export_pause);
    end
end

function [Ymat, cleanData] = build_y_matrix(dataCell)
%BUILD_Y_MATRIX Convert one vector per group into a numeric matrix with NaN padding.
% Values <= 0 are removed because plots are displayed on a log-scale.
    nGroups = numel(dataCell);
    cleanData = cell(1,nGroups);
    lens = zeros(nGroups,1);

    for i = 1:nGroups
        yi = dataCell{i}(:);
        yi = yi(isfinite(yi) & yi > 0);
        cleanData{i} = yi;
        lens(i) = numel(yi);
    end

    maxN = max(lens);
    if isempty(maxN) || maxN == 0
        error('No positive finite values available for violin plot.');
    end

    Ymat = NaN(maxN,nGroups);
    for i = 1:nGroups
        yi = cleanData{i};
        Ymat(1:numel(yi),i) = yi;
    end
end

function shrink_violin_patches(ax, widthFactor)
%SHRINK_VIOLIN_PATCHES Reduce violin widths by scaling Patch XData around each center.
    if nargin < 2 || isempty(widthFactor)
        widthFactor = 0.7;
    end
    widthFactor = max(0.05, widthFactor);

    patches = findobj(ax,'Type','Patch');
    for k = 1:numel(patches)
        x = patches(k).XData;
        if isempty(x), continue; end
        xc = mean([min(x) max(x)]);
        patches(k).XData = xc + widthFactor*(x - xc);
    end
end

function apply_y_axis_ticks(ax, Ymat, preferredTicks, yScaleMode)
%APPLY_Y_AXIS_TICKS Apply publication-ready y-axis limits and labels.
% yScaleMode:
%   'log'          -> decimal ticks on log scale.
%   'linear_auto3' -> linear scale with three informative ticks:
%                     5th percentile, mean, 95th percentile.

    if nargin < 4 || isempty(yScaleMode)
        yScaleMode = 'log';
    end

    y = Ymat(:);
    y = y(isfinite(y) & y > 0);
    if isempty(y), return; end

    switch lower(char(yScaleMode))
        case 'linear_auto3'
            set(ax,'YScale','linear');

            ymin = min(y);
            ymax = max(y);
            yr = ymax - ymin;
            if yr <= 0
                yr = max(abs(ymax),1);
            end
            % Generous vertical padding for the single plasmid-copy violin so
            % the distribution and 5-95 overlay do not feel cramped.
            ylim(ax,[ymin - 0.25*yr, ymax + 0.30*yr]);

            % Three informative ticks: lower tail, central value, upper tail.
            % The lower/upper ticks avoid the absolute extremes, which are often
            % visually too close to the violin tips.
            tLow  = prctile(y,5);
            tMid  = mean(y);
            tHigh = prctile(y,95);

            ticks = unique([tLow tMid tHigh],'stable');
            if numel(ticks) < 3
                ticks = linspace(ymin,ymax,3);
            end

            ax.YTick = ticks;
            ax.YTickLabel = arrayfun(@natural_tick_label, ticks, 'UniformOutput', false);

        otherwise
            set(ax,'YScale','log');

            ymin = min(y);
            ymax = max(y);

            % Padding in log space.
            padFactor = 1.35;
            ylim(ax,[ymin/padFactor ymax*padFactor]);

            if nargin < 3 || isempty(preferredTicks)
              %  baseTicks = [0.001 0.002 0.005 0.01 0.02 0.05 0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000];
                baseTicks = [0.002 0.005 0.01  0.05 0.1 0.2 0.5 1 2 5 10 20 50 100 200 500 1000];
            else
                baseTicks = preferredTicks(:).';
            end

            yl = ylim(ax);
            ticks = baseTicks(baseTicks >= yl(1) & baseTicks <= yl(2));

            if numel(ticks) < 2
                lo = floor(log10(yl(1)));
                hi = ceil(log10(yl(2)));
                ticks = logspace(lo,hi,max(3,hi-lo+1));
            end

            ax.YTick = ticks;
            ax.YTickLabel = arrayfun(@natural_tick_label, ticks, 'UniformOutput', false);
    end
end

function s = natural_tick_label(x)
%NATURAL_TICK_LABEL Compact decimal labels without scientific notation.
    if x >= 100
        s = sprintf('%.0f',x);
    elseif x >= 10
        s = sprintf('%.0f',x);
    elseif x >= 1
        s = sprintf('%.2g',x);
    elseif x >= 0.1
        s = sprintf('%.2f',x);
    elseif x >= 0.01
        s = sprintf('%.3f',x);
    else
        s = sprintf('%.4f',x);
    end

    while contains(s,'.') && endsWith(s,'0')
        s = extractBefore(s,strlength(s));
    end
    if endsWith(s,'.')
        s = extractBefore(s,strlength(s));
    end
    s = char(s);
end

function c = darken_color(c, amount)
%DARKEN_COLOR Mix color c with black. amount=0 keeps c; amount=1 gives black.
    amount = max(0,min(1,amount));
    c = (1-amount).*c;
end

function save_figure_publication(fig,outdir_or_basename,basename_or_resolution,resolution_or_pause,maybe_pause)
%SAVE_FIGURE_PUBLICATION Export PNG, SVG, and PDF robustly.
% Supports:
%   save_figure_publication(fig, basename, resolution, pause_s)
%   save_figure_publication(fig, outdir, basename, resolution, pause_s)

    if nargin == 4
        basename = outdir_or_basename;
        resolution = basename_or_resolution;
        export_pause = resolution_or_pause;
    elseif nargin == 5
        outdir = outdir_or_basename;
        basename = fullfile(outdir, basename_or_resolution);
        resolution = resolution_or_pause;
        export_pause = maybe_pause;
    else
        error('save_figure_publication:InvalidInput','Unexpected number of inputs.');
    end

    [folder,~,~] = fileparts(basename);
    if ~isempty(folder) && ~exist(folder,'dir')
        mkdir(folder);
    end

    drawnow;
    pause(export_pause);

    pngfile = [basename '.png'];
    svgfile = [basename '.svg'];
    pdffile = [basename '.pdf'];

    try
        exportgraphics(fig,pngfile, ...
            'Resolution',resolution, ...
            'BackgroundColor','white');
    catch
        print(fig,pngfile,'-dpng',sprintf('-r%d',resolution));
    end

    drawnow;
    pause(export_pause);

    try
        set(fig,'Renderer','painters');
    catch
    end

    try
        exportgraphics(fig,svgfile, ...
            'ContentType','vector', ...
            'BackgroundColor','white');
    catch
        print(fig,svgfile,'-dsvg','-painters');
    end

    drawnow;
    pause(export_pause);

    try
        exportgraphics(fig,pdffile, ...
            'ContentType','vector', ...
            'BackgroundColor','white');
    catch
        print(fig,pdffile,'-dpdf','-painters');
    end

    fprintf('Saved figures:\n  %s\n  %s\n  %s\n', ...
        pngfile, svgfile, pdffile);
end

function rgb = hex2rgb_local(hex)
%HEX2RGB_LOCAL Convert #RRGGBB string(s) to RGB rows in [0,1].
    if isstring(hex)
        hex = char(hex);
    end
    if iscell(hex)
        hex = char(string(hex));
    end

    if size(hex,1) > 1
        rgb = zeros(size(hex,1),3);
        for i = 1:size(hex,1)
            rgb(i,:) = hex2rgb_local(strtrim(hex(i,:)));
        end
        return;
    end

    hex = strtrim(hex);
    if startsWith(hex,'#')
        hex = hex(2:end);
    end
    if numel(hex) ~= 6
        error('hex2rgb_local:InvalidHex','Hex color must be 6 characters: RRGGBB.');
    end

    rgb = [hex2dec(hex(1:2)) hex2dec(hex(3:4)) hex2dec(hex(5:6))] ./ 255;
end
