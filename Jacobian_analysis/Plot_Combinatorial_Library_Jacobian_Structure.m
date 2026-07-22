function Plot_Combinatorial_Library_Jacobian_Structure()
%PLOT_COMBINATORIAL_LIBRARY_JACOBIAN_STRUCTURE
% Visualize the combinatorial parameter--construct incidence structure of
% the complete L24 + L6 + L5 library.
%
% Columns correspond to 35 expression constructs. Rows correspond to global
% promoter, RBS, and plasmid parameters. Filled cells indicate structural
% dependence, not numerical sensitivity magnitude.
%
% Outputs:
%   Combinatorial_Library_Jacobian_Structure.png
%   Combinatorial_Library_Jacobian_Structure.svg
%
% See README_Jacobian_analysis.md for context.

%% =============================================================
% 1. Full library: L24 + L6 + L5
% =============================================================

IDs = { ...
    'L24_01','L24_02','L24_03','L24_04', ...
    'L24_05','L24_06','L24_07','L24_08', ...
    'L24_09','L24_10','L24_11','L24_12', ...
    'L24_13','L24_14','L24_15','L24_16', ...
    'L24_17','L24_18','L24_19','L24_20', ...
    'L24_21','L24_22','L24_23','L24_24', ...
    'L6_1','L6_2','L6_3','L6_4','L6_5','L6_6', ...
    'L5_1','L5_2','L5_4','L5_5','L5_3'};

ORIs = { ...
    'pSC101','pSC101','pSC101','pSC101', ...
    'pSC101','pSC101','pSC101','pSC101', ...
    'pSC101','pSC101','pSC101','pSC101', ...
    'pGreen','pGreen','pGreen','pGreen', ...
    'pGreen','pGreen','pGreen','pGreen', ...
    'pGreen','pGreen','pGreen','pGreen', ...
    'pSC101','pSC101','pSC101','pGreen','pGreen','pGreen', ...
    'pGreen','pGreen','pGreen','pGreen','pGreen'};

Promoters = { ...
    'J23106','J23106','J23106','J23106', ...
    'J23102','J23102','J23102','J23102', ...
    'J23101','J23101','J23101','J23101', ...
    'J23106','J23106','J23106','J23106', ...
    'J23102','J23102','J23102','J23102', ...
    'J23101','J23101','J23101','J23101', ...
    'J23106','J23102','J23101','J23106','J23102','J23101', ...
    'J23100','J23100','J23100','J23100','J23100'};

RBSs = { ...
    'B0030','B0032','J61100','J61101', ...
    'B0030','B0032','J61100','J61101', ...
    'B0030','B0032','J61100','J61101', ...
    'B0030','B0032','J61100','J61101', ...
    'B0030','B0032','J61100','J61101', ...
    'B0030','B0032','J61100','J61101', ...
    'B0034','B0034','B0034','B0034','B0034','B0034', ...
    'B0030','B0032','J61100','J61101','B0034'};

nTUs = numel(IDs);

%% =============================================================
% 2. Global parameters
% =============================================================

proms = unique(Promoters, 'stable');
rbss  = unique(RBSs,      'stable');
oris  = unique(ORIs,      'stable');

paramKeys = {};
paramLabels = {};

for i = 1:numel(proms)
    paramKeys{end+1}   = ['omega_' proms{i}]; %#ok<AGROW>
    paramLabels{end+1} = ['$\omega_{\mathrm{' proms{i} '}}$']; %#ok<AGROW>
end
idxEnd_omega = numel(paramKeys);

for i = 1:numel(rbss)
    paramKeys{end+1}   = ['kappa_' rbss{i}]; %#ok<AGROW>
    paramLabels{end+1} = ['$\kappa^{0}_{\mathrm{' rbss{i} '}}$']; %#ok<AGROW>
end
idxEnd_kappa = numel(paramKeys);

for i = 1:numel(rbss)
    paramKeys{end+1}   = ['rho_' rbss{i}]; %#ok<AGROW>
    paramLabels{end+1} = ['$\rho^{0}_{\mathrm{' rbss{i} '}}$']; %#ok<AGROW>
end
idxEnd_rho = numel(paramKeys);

for i = 1:numel(oris)
    paramKeys{end+1}   = ['N_' oris{i}]; %#ok<AGROW>
    paramLabels{end+1} = ['$N_{\mathrm{' oris{i} '}}$']; %#ok<AGROW>
end

nParams = numel(paramKeys);

%% =============================================================
% 3. Matrix: 35 ECs x 16 parameters
% =============================================================

M_full = zeros(nTUs, nParams);

for i = 1:nTUs
    prom = Promoters{i};
    rbs  = RBSs{i};
    ori  = ORIs{i};

    M_full(i, strcmp(paramKeys, ['omega_' prom])) = 1;
    M_full(i, strcmp(paramKeys, ['kappa_' rbs]))  = 1;
    M_full(i, strcmp(paramKeys, ['rho_' rbs]))    = 1;
    M_full(i, strcmp(paramKeys, ['N_' ori]))      = 1;
end

% Transposed layout for the figure:
%   rows    = parameters
%   columns = ECs
M_plot = M_full.';

%% =============================================================
% 4. Plot
% =============================================================

color_grey = [184 185 182]/255;
cmap = [1 1 1; color_grey];

fig_width_cm  = 42;
fig_height_cm = 26;

fs_ticks_x = 18;
fs_ticks_y = 23;
fs_axis    = 28;
fs_title   = 30;

fig = figure('Name','Combinatorial library Jacobian structure', ...
    'Color','w', ...
    'Units','centimeters', ...
    'Position',[1 1 fig_width_cm fig_height_cm], ...
    'PaperUnits','centimeters', ...
    'PaperPosition',[0 0 fig_width_cm fig_height_cm], ...
    'PaperSize',[fig_width_cm fig_height_cm], ...
    'Renderer','painters');

ax = axes(fig, ...
    'Units','normalized', ...
    'Position',[0.13 0.16 0.85 0.74]);

imagesc(ax, M_plot);
colormap(ax, cmap);
clim(ax, [0 1]);
hold(ax,'on');

axis(ax,'normal');

% Separators between libraries on X axis: L24 | L6 | L5.
for xb = [24, 30]
    line(ax, [xb+0.5 xb+0.5], [0.5 nParams+0.5], ...
        'Color', [0.20 0.20 0.20], 'LineWidth', 3.0);
end

% Separators between parameter families on Y axis.
for yb = [idxEnd_omega, idxEnd_kappa, idxEnd_rho]
    line(ax, [0.5 nTUs+0.5], [yb+0.5 yb+0.5], ...
        'Color', [0.45 0.45 0.45], 'LineWidth', 2.4);
end

ax.XTick = 1:nTUs;
ax.XTickLabel = escape_underscore_labels(IDs);
ax.XTickLabelRotation = 90;

ax.YTick = 1:nParams;
ax.YTickLabel = paramLabels;

ax.TickLength = [0 0];
ax.LineWidth = 1.4;
ax.Box = 'on';
ax.Layer = 'top';
ax.TickLabelInterpreter = 'latex';
ax.FontSize = fs_ticks_y;

xlabel(ax, 'ECs', ...
    'FontSize', fs_axis, ...
    'FontWeight','normal', ...
    'Interpreter','none');

ylabel(ax, 'Parameters', ...
    'FontSize', fs_axis, ...
    'FontWeight','normal', ...
    'Interpreter','none');

title(ax, '', ...
    'FontSize', fs_title, ...
    'FontWeight','bold', ...
    'Interpreter','none');

xlim(ax, [0.5 nTUs+0.5]);
ylim(ax, [0.5 nParams+0.5]);

% Add small library labels above the X axis blocks.
text(ax, 12.5, 0.35, 'L24', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontSize', 24, ...
    'FontWeight','bold', ...
    'Interpreter','none');

text(ax, 27.5, 0.35, 'L6', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontSize', 24, ...
    'FontWeight','bold', ...
    'Interpreter','none');

text(ax, 33.0, 0.35, 'L5', ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','bottom', ...
    'FontSize', 24, ...
    'FontWeight','bold', ...
    'Interpreter','none');

exportgraphics(fig, 'Combinatorial_Library_Jacobian_Structure.png', 'Resolution', 600);
exportgraphics(fig, 'Combinatorial_Library_Jacobian_Structure.svg', 'ContentType', 'vector');

end

function labels_out = escape_underscore_labels(labels_in)
    labels_out = labels_in;
    for ii = 1:numel(labels_in)
        labels_out{ii} = strrep(labels_in{ii}, '_', '\_');
    end
end
