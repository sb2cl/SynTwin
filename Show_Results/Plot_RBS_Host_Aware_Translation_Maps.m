%% Plot_RBS_Host_Aware_Translation_Maps.m
% Plot host-aware mappings from RBS intrinsic properties to effective
% translation quantities.
%
% Associated publication:
%   "Host-aware Identification of Intrinsic Gene Expression Biopart
%    Parameters using Combinatorial Libraries"
%
% Main-text association:
%   - Figure 4a: projection from intrinsic initiation capacity (IIC,
%     kappa^0) to the effective ribosome--RBS association constant
%     (ERAC = K_A^t/nu_t) as a function of growth rate.
%   - Figure 4b: predicted operational ranges of the effective translation
%     rate (ETR = K_A^t) for the five characterised RBS variants.
%
% The script does not rerun parameter estimation. It loads the complete L24
% and L6 result tensors, extracts stored digital-twin predictions and
% experimental summaries, and regenerates the published mappings.
%
% Inputs:
%   Estimation_Pi/L24_L1O_reduced_model/
%       Results_Tensor_Lib24_L1O_reduced_Wells.mat
%   Estimation_Pi/L6_L1O_reduced_model/
%       Results_Tensor_Lib6_L1O_reduced_Wells.mat
%
% Outputs in ./Figures:
%   Figure_4a_ERAC_projection.png/.svg/.pdf
%   Figure_4b_ETR_range_by_IIC.png/.svg/.pdf
%
% Requirements:
%   - MATLAB with scatteredInterpolant, exportgraphics, and local functions
%     in scripts.
%
% Usage:
%   Run Plot_RBS_Host_Aware_Translation_Maps from any working directory
%   while the script remains inside the SynTwin repository.
%
% See README_Show_Results.md for the complete Show_Results workflow.

clearvars;
close all;
clc;

% Use a consistent sans-serif font for all non-mathematical text.
% Mathematical symbols rendered with the LaTeX interpreter retain the
% appropriate mathematical typography.
set(groot, ...
    'defaultAxesFontName','Arial', ...
    'defaultTextFontName','Arial', ...
    'defaultLegendFontName','Arial');

%% ---- Portable package paths ----
ROOT = init_SynTwin();
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);
FIGURES_ROOT = fileparts(SCRIPT_DIR);

RBS_colors_hex = ['#A8DADC';'#F4A261';'#B5E48C';'#CDB4DB';'#FFE066'];
RBS_colors = hex2rgb_local(RBS_colors_hex);

%% ---- USER SETTINGS ----

required_files = { ...
    'Results_Tensor_Lib24_L1O_reduced_Wells.mat', ...
    'Results_Tensor_Lib6_L1O_reduced_Wells.mat'};

% Getting the data of Lib24 L1O reduced model: 
S24 = load(SynTwin_path('Estimation_Pi/L24_L1O_reduced_model','Results_Tensor_Lib24_L1O_reduced_Wells.mat'), ...
    'Results_Tensor_Lib24_L1O_reduced');
% Getting the estimated IIC (K/sigma) of RBS B0034 from Lib6: 
S6 = load(SynTwin_path('Estimation_Pi/L6_L1O_reduced_model','Results_Tensor_Lib6_L1O_reduced_Wells.mat'), ...
    'Results_Tensor_Lib6_L1O_reduced');

assert(isfield(S24,'Results_Tensor_Lib24_L1O_reduced'), ...
    'Variable Results_Tensor_Lib24_L1O_reduced is missing from %s.', required_files{1});
assert(isfield(S6,'Results_Tensor_Lib6_L1O_reduced'), ...
    'Variable Results_Tensor_Lib6_L1O_reduced is missing from %s.', required_files{2});

Results_Tensor_Lib24_L1O_reduced = S24.Results_Tensor_Lib24_L1O_reduced;
Results_Tensor_Lib6_L1O_reduced  = S6.Results_Tensor_Lib6_L1O_reduced;

% Constants (set to your model constants)
nu_max = 20.5156*60; % aa/min (as you stated; units cancel in normalization)
% If you want nu_max in min^-1 scaling, keep consistent with your KAt units.
% Since you normalize by nu_max*f, nu_max should match whatever you used when computing KAt.
% If KAt already includes nu_max in min^-1, keep nu_max in min^-1 as used in code.
% If unsure, set nu_max = 1 and interpret Z as KAt/f(s_n).

% kappa estimates (kappa = K0/sigma0) by RBS name:
kappa_map = struct();
kappa_map.B0030  = Results_Tensor_Lib24_L1O_reduced{1,1,1}.RBS_k0_sigma0_mean;
kappa_map.B0032  = Results_Tensor_Lib24_L1O_reduced{1,1,2}.RBS_k0_sigma0_mean;
kappa_map.B0034  = Results_Tensor_Lib6_L1O_reduced{1,1}.RBS_k0_sigma0_mean;
kappa_map.J61100 = Results_Tensor_Lib24_L1O_reduced{1,1,3}.RBS_k0_sigma0_mean;
kappa_map.J61101 = Results_Tensor_Lib24_L1O_reduced{1,1,4}.RBS_k0_sigma0_mean;

model_c.dm_c  =  0.2; %Mean degradation rate of non-ribosomal mRNA (1/min)


%  TU structs encode the RBS name in a field:
rbs_field_candidates = {'TU_RBS'};

% ---- Output settings for main 2D figures ----
save_main_figures = true;
fig_output_dir = fullfile(SCRIPT_DIR, 'Figures');
if save_main_figures && ~exist(fig_output_dir,'dir')
    mkdir(fig_output_dir);
end

% ---- Aesthetic settings for publication figures ----
% These settings follow the same publication-quality logic used in
% Plot_Results_Lib_reduced_model_v9.m: figures are created in centimeters,
% then exported to PNG/SVG. Increase Fig*SizeCm to create a larger canvas
% before final LaTeX scaling; increase font_* and marker_exp if labels or
% experimental points remain too small after reduction in the manuscript.
font_axis  = 28;
font_label = 34;
font_title = 28;
font_rbs   = 28;
line_iso   = 3.2;
line_width = 3.0;
marker_exp = 155;
marker_exp_3d = 130;
export_resolution = 600;
export_pause = 1.2;

% Publication canvas sizes. These are intentionally larger than final
% printed size: the figures are exported on a large canvas and then scaled
% down in LaTeX, preserving readable labels and visible experimental points.
% Increase/decrease the second number to tune vertical aspect ratio.
FigProjectionSizeCm = [64 38];
FigRangeSizeCm      = [60 30];
% Lighter colormap: white -> light blue -> green -> yellow.
map_light = make_light_iic_colormap(256);

%% ---- Extract stored predictions and experimental summaries ----
D1 = Results_Tensor_Lib24_L1O_reduced;
D2 = Results_Tensor_Lib6_L1O_reduced;

objs1 = find_objects_with_predictions(D1);
objs2 = find_objects_with_predictions(D2);

all_objs = [objs1; objs2];
if isempty(all_objs)
    error('No objects with a Synthesis_predictions field were found in the loaded .mat files.');
end

% Aggregate into arrays
MU = []; MU_EXP = []; 
KAPPA = [];  RBSNAME = strings(0,1); THETA = []; VARPHI = []; FS = []; 
KAT = []; KAT_EXP = []; 
Z = []; Z_EXP = [];

for i = 1:numel(all_objs)
    obj = all_objs{i};
    SP = obj.Synthesis_predictions;

    % Required fields check
    req = {'Mu_values','fs_pred_values','KA_t_pred_values'};
    for k = 1:numel(req)
        if ~isfield(SP, req{k})
            error('Missing field Synthesis_predictions.%s in object #%d.', req{k}, i);
        end
    end

    mu = SP.Mu_values(:);
    fs = SP.fs_pred_values(:);
    kat = SP.KA_t_pred_values(:);
    if ~isfield(obj,'Exp_Data')
        pi_exp = obj.Pi_mumax_pmax_global_mean;
        mu_exp = obj.Mu_mumax_pmax_global_mean;
    else
        pi_exp = obj.Exp_Data.Pi_mumax_pmax_global_mean;
        mu_exp = obj.Exp_Data.Mu_mumax_pmax_global_mean;
    end

    % Optional: load, varphi
    if isfield(SP,'load_pred_values'), theta = SP.load_pred_values(:); else, theta = nan(size(mu)); end
    if isfield(SP,'varphi_pred_values'), varphi = SP.varphi_pred_values(:); else, varphi = nan(size(mu)); end

    % Determine RBS name
    rbs_name = extract_rbs_name(obj, rbs_field_candidates);
    if strlength(rbs_name)==0
        error('Plot_RBS_Host_Aware_Translation_Maps:MissingRBSName', ...
            'TU_RBS was not found in result object #%d.', i);
    end

    if ~isfield(kappa_map, char(rbs_name))
        warning('Unknown RBS "%s" in object #%d. Skipping.', rbs_name, i);
        continue;
    end
    kappa = kappa_map.(char(rbs_name));

    % Compute Z = KAt / nu(sn)
    % nu(sn) = nu_max * f(s_n)
    % If you are uncertain about unit consistency, set nu_max = 1 above.
    nu = nu_max .* fs;
    z = kat ./ nu;
    pi_exp_interp = interpolate_pi_from_mu_log(mu_exp, pi_exp, mu, 'TolAbs', 1e-3);
    kat_exp = pi_exp_interp/obj.Gene_cn_mean/obj.Omega_mean*model_c.dm_c; 
    z_exp = kat_exp ./ nu;

    % Filter invalid points
    good = isfinite(mu) & isfinite(fs) & isfinite(kat)  & isfinite(kat_exp) & isfinite(z) & isfinite(z_exp) & (mu>0) & (fs>0) & (kat>0) & (z>0);

    MU     = [MU;     mu(good)];
    FS     = [FS;     fs(good)];
    KAT    = [KAT;    kat(good)];
    KAT_EXP    = [KAT_EXP;    kat_exp(good)];
    Z      = [Z;      z(good)];
    Z_EXP      = [Z_EXP;      z_exp(good)];
    THETA  = [THETA;  theta(good)];
    VARPHI = [VARPHI; varphi(good)];
    KAPPA  = [KAPPA;  repmat(kappa, sum(good), 1)];
    RBSNAME= [RBSNAME; repmat(string(rbs_name), sum(good), 1)];
end

if isempty(MU)
    error('No valid data points after filtering.');
end

%% FIGURE 4a (MAIN): 2D projection of the 3D mapping onto (mu, ERAC) with colour = kappa
% Axes: x = mu (min^-1), y = ERAC = K_A^t/nu(s_n), colour = kappa = K0/sigma0 (IIC)
% v3: lighter semi-transparent colour layer + larger fonts + PNG/SVG export.

maskF = isfinite(MU) & isfinite(KAPPA) & isfinite(Z);   % Z = ERAC pred = Kt/nu
[xF, yF, zF] = average_duplicate_xy(KAPPA(maskF), MU(maskF), Z(maskF));
Fz = scatteredInterpolant(xF, yF, zF, 'natural', 'nearest');

mu_min_plot = 0.0079;
mu_max_plot = 0.021;
kappa_min_plot = min(KAPPA(maskF)) * 0.8;
kappa_max_plot = max(KAPPA(maskF)) * 1.2;

mu_fine    = linspace(mu_min_plot, mu_max_plot, 180);
kappa_fine = logspace(log10(kappa_min_plot), log10(kappa_max_plot), 220);
[KKf, MMf] = meshgrid(kappa_fine, mu_fine);
ZZf = Fz(KKf, MMf);   % ERAC surface

fig3 = figure('Color','w','Name','ERAC projection', ...
    'Units','centimeters','Position',[1 1 FigProjectionSizeCm]);
ax3 = axes('Parent',fig3);
hold(ax3,'on');
set(fig3,'Renderer','opengl');
set(fig3,'PaperPositionMode','auto');

hs3 = surf(ax3, MMf, ZZf, zeros(size(ZZf)), KKf, ...
    'EdgeColor','none', 'FaceAlpha',0.72);
ax3 = ancestor(hs3,'axes');
view(ax3,2);
apply_colormap_safe(ax3,map_light);
set(ax3, 'YScale','log');
box(ax3,'on'); grid(ax3,'on');
set(ax3,'Layer','top','GridAlpha',0.18,'MinorGridAlpha',0.08);
ax3.FontName = 'Arial';

xlabel(ax3,'$\mu\; (\mathrm{min}^{-1})$', ...
    'FontName','Arial', ...
    'FontSize',font_label, ...
    'Interpreter','latex');
ylabel(ax3,'$\mathrm{ERAC}\; (\mathrm{molec}^{-1} )$', ...
    'FontName','Arial', ...
    'FontSize',font_label, ...
    'Interpreter','latex');

cb = colorbar(ax3);
cb.FontName = 'Arial';
cb.Label.String = '$\kappa^0\,  (\mathrm{molec}^{-1})$';
cb.Label.FontName = 'Arial';
cb.Label.FontSize = 1.2*font_label;
cb.Label.Interpreter = 'latex';
cb.FontSize = font_axis;

xlim(ax3,[mu_min_plot mu_max_plot]);
ylim(ax3,[1e-3 1]);

maskPts = isfinite(MU) & isfinite(KAPPA) & isfinite(Z_EXP);
scatter(ax3, MU(maskPts), Z_EXP(maskPts), marker_exp, KAPPA(maskPts), ...
    'filled', 'MarkerFaceAlpha',0.70, ...
    'MarkerEdgeAlpha',0.35, 'MarkerEdgeColor',[0.15 0.15 0.15]);

rbs_plot_order = ["B0030","B0034","B0032","J61100","J61101"];
label_cfg = rbs_label_config();
for nm = rbs_plot_order
    if isfield(kappa_map, char(nm))
        k_i = kappa_map.(char(nm));
        z_i = Fz(k_i*ones(size(mu_fine)), mu_fine);
        validLine = isfinite(z_i) & z_i > 0;
        if any(validLine)
            plot(ax3, mu_fine(validLine), z_i(validLine), 'k--', 'LineWidth',line_iso);
            cfg = label_cfg.(char(nm));
            z_lab = interp1(mu_fine(validLine), z_i(validLine), cfg.x, 'linear', 'extrap');
            if isfinite(z_lab) && z_lab > 0
                text(ax3, cfg.x, cfg.yfac*z_lab, char(nm), ...
                    'FontName', 'Arial', ...
                    'FontWeight', 'bold', 'FontSize', font_rbs, 'Color', 'k', ...
                    'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', ...
                    'Clipping', 'off');
            end
        end
    end
end

ax = ax3;
ax.FontName = 'Arial';
ax.XAxis.FontSize = font_axis;
ax.YAxis.FontSize = font_axis;
ax.XLabel.FontSize = font_label;
ax.YLabel.FontSize = font_label;

t3 = title(ax3,'Projection IIC ($\kappa^0$) $\rightarrow$ effective RBS--ribosome association (ERAC)', ...
    'FontName','Arial', ...
    'FontSize',1.2*font_title, ...
    'Interpreter','latex');
t3.FontWeight = 'bold';

uistack(hs3,'bottom');

if save_main_figures
    save_figure_publication(fig3, fig_output_dir, 'Figure_4a_ERAC_projection', export_resolution, export_pause);
end

%% FIGURE 4b (MAIN): ETR predicted range for each RBS IIC
% For each RBS, show the range of predicted ETR = K_A^t across the plotted
% mu interval. The x-axis is logarithmic because ETR spans orders of
% magnitude. The y-axis labels report both the RBS name and its IIC value.
fig5 = figure('Color','w','Name','ETR predicted range by IIC', ...
    'Units','centimeters','Position',[1 1 FigRangeSizeCm]);
ax5 = axes('Parent',fig5);
hold(ax5,'on');
set(fig5,'Renderer','opengl');
set(fig5,'PaperPositionMode','auto');

rbs_range_order_all = ["B0032","J61100","J61101","B0034","B0030"];
rbs_range_order = strings(0,1);
for rr0 = 1:numel(rbs_range_order_all)
    if isfield(kappa_map,char(rbs_range_order_all(rr0)))
        rbs_range_order(end+1,1) = rbs_range_order_all(rr0); %#ok<SAGROW>
    end
end

ypos = 1:numel(rbs_range_order);
ylabels = strings(numel(rbs_range_order),1);
xrange_all = [];
for rr = 1:numel(rbs_range_order)
    nm = rbs_range_order(rr);
    kv = kappa_map.(char(nm));
    idx = RBSNAME == nm & isfinite(KAT) & KAT > 0 & isfinite(MU) & ...
          MU >= mu_min_plot & MU <= mu_max_plot;
    if ~any(idx)
        idx = RBSNAME == nm & isfinite(KAT) & KAT > 0;
    end
    vals = KAT(idx);
    if isempty(vals)
        continue;
    end
    xmin = min(vals);
    xmax = max(vals);
    xmed = median(vals,'omitnan');
    xrange_all = [xrange_all; xmin; xmax]; %#ok<AGROW>
    C = rbs_color_from_name(nm, RBS_colors);
    plot(ax5, [xmin xmax], [ypos(rr) ypos(rr)], '-', 'Color', C, 'LineWidth', 8);
    plot(ax5, xmed, ypos(rr), 'o', 'MarkerSize', 11, ...
        'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'LineWidth',1.8);
    plot(ax5, [xmin xmax], [ypos(rr) ypos(rr)], '|', 'Color','k', 'LineWidth',2.0, 'MarkerSize',16);
    ylabels(rr) = sprintf('%s  ($\\kappa^0=%.3g$)', char(nm), kv);
end
set(ax5,'XScale','log');
if ~isempty(xrange_all)
    xrange_all = xrange_all(isfinite(xrange_all) & xrange_all > 0);
    if ~isempty(xrange_all)
        xlim_pair = [0.55*min(xrange_all) 2.20*max(xrange_all)];
        xlim(ax5,xlim_pair);

        % Increase the tick resolution on the logarithmic ETR axis.
        % Use 1-2-5 ticks per decade, plus decade ticks, so the reader can
        % see where each RBS range lies without overcrowding the axis.
        [xtick_vals, xtick_lbls] = log_ticks_125(xlim_pair(1), xlim_pair(2));
        ax5.XTick = xtick_vals;
        ax5.XTickLabel = xtick_lbls;
        ax5.XMinorTick = 'on';
        ax5.XMinorGrid = 'on';
        ax5.MinorGridAlpha = 0.12;
        ax5.GridAlpha = 0.28;
    end
end
grid(ax5,'on'); box(ax5,'on');
ylim(ax5,[0.5 numel(rbs_range_order)+0.5]);
yticks(ax5,ypos);
yticklabels(ax5,cellstr(ylabels));
xlabel(ax5,'$\mathrm{ETR}\; (\mathrm{min}^{-1})$', ...
    'FontName','Arial', ...
    'FontSize',font_label, ...
    'Interpreter','latex');
ylabel(ax5,'RBS IIC', ...
    'FontName','Arial', ...
    'FontSize',font_label, ...
    'Interpreter','latex');
set_axis_publication(ax5,font_axis);
ax5.YTickLabel = cellstr(ylabels);
ax5.TickLabelInterpreter = 'latex';
if save_main_figures
    save_figure_publication(fig5, fig_output_dir, 'Figure_4b_ETR_range_by_IIC', export_resolution, export_pause);
end

%% ---- DONE ----
disp('Host-aware RBS translation maps generated successfully.');

%% ===================== HELPER FUNCTIONS =====================

function assert_required_file(filePath)
%ASSERT_REQUIRED_FILE Fail with an informative message when input is absent.
    if ~isfile(filePath)
        error('Plot_RBS_Host_Aware_Translation_Maps:MissingInput', ...
            'Required data file not found: %s', filePath);
    end
end

function assert_required_variable(S,varName,filePath)
%ASSERT_REQUIRED_VARIABLE Verify the expected variable is stored in a MAT file.
    if ~isfield(S,varName)
        error('Plot_RBS_Host_Aware_Translation_Maps:MissingVariable', ...
            'Variable %s was not found in %s.', varName, filePath);
    end
end


function [ticks, labels] = log_ticks_125(xlo, xhi)
%LOG_TICKS_125 Readable 1-2-5 ticks for logarithmic axes.
% Returns ticks at 1, 2 and 5 times each power of ten inside [xlo,xhi],
% with compact labels suitable for publication figures.

    if ~isfinite(xlo) || ~isfinite(xhi) || xlo <= 0 || xhi <= 0 || xlo >= xhi
        ticks = [];
        labels = {};
        return;
    end

    e0 = floor(log10(xlo));
    e1 = ceil(log10(xhi));
    bases = [1 2 5];
    ticks = [];
    for e = e0:e1
        ticks = [ticks, bases.*10.^e]; %#ok<AGROW>
    end
    ticks = unique(ticks(ticks >= xlo & ticks <= xhi));

    % Ensure the plot does not end up with too many tick labels if the range
    % is unexpectedly wide. Prefer 1-2-5 spacing, then fall back to decades.
    if numel(ticks) > 14
        ticks_dec = 10.^(e0:e1);
        ticks_dec = ticks_dec(ticks_dec >= xlo & ticks_dec <= xhi);
        ticks = ticks_dec;
    end

    labels = cell(size(ticks));
    for ii = 1:numel(ticks)
        x = ticks(ii);
        if x >= 1000 || x < 0.01
            labels{ii} = sprintf('%.0e', x);
        elseif x >= 100
            labels{ii} = sprintf('%.0f', x);
        elseif x >= 10
            labels{ii} = sprintf('%.0f', x);
        elseif x >= 1
            labels{ii} = sprintf('%.1g', x);
        else
            labels{ii} = sprintf('%.2g', x);
        end
    end
end

function objs = find_objects_with_predictions(S)
% Return a cell array with all TU-like structs that contain a
% 'Synthesis_predictions' field. Works whether S is a scalar struct,
% a struct array (any size), or a nested struct containing such arrays.

objs = {};

% Case 1: S itself is a struct array of TU entries
if isstruct(S) && isfield(S,'Synthesis_predictions')
    % Return each element (so we can keep per-TU metadata such as TU_RBS)
    tmp = num2cell(S(:));
    objs = [objs; tmp(:)];
    return;
end

% Case 2: recursive descent into fields / elements
objs = descend(S);

    function out = descend(x)
        out = {};
        if isstruct(x)
            % If x is a struct array, visit each element
            for ii = 1:numel(x)
                out = [out; descend_struct_scalar(x(ii))]; %#ok<AGROW>
            end
        elseif iscell(x)
            for ii = 1:numel(x)
                out = [out; descend(x{ii})]; %#ok<AGROW>
            end
        end
    end

    function out = descend_struct_scalar(st)
        out = {};
        if isfield(st,'Synthesis_predictions')
            out{end+1,1} = st; %#ok<AGROW>
            return;
        end
        fns = fieldnames(st);
        for jj = 1:numel(fns)
            try
                v = st.(fns{jj});
            catch
                continue;
            end
            out = [out; descend(v)]; %#ok<AGROW>
        end
    end

end


function rbs_name = extract_rbs_name(obj, candidates)
% Try to extract an RBS name string from obj using candidate field names.
    rbs_name = "";
    for i = 1:numel(candidates)
        f = candidates{i};
        if isfield(obj, f)
            val = obj.(f);
            if ischar(val) || isstring(val)
                rbs_name = string(val);
                return;
            end
        end
    end
    % Sometimes nested: obj.Current_TU.RBS, etc.
    if isfield(obj, 'Current_TU')
        ct = obj.Current_TU;
        for i = 1:numel(candidates)
            f = candidates{i};
            if isstruct(ct) && isfield(ct, f)
                val = ct.(f);
                if ischar(val) || isstring(val)
                    rbs_name = string(val);
                    return;
                end
            end
        end
    end
end





function pi_interp = interpolate_pi_from_mu_log(mu_exp, pi_exp, mu_query, varargin)
%INTERPOLATE_PI_FROM_MU_LOG Interpola pi_exp(mu_exp) en mu_query SIN extrapolar.
%   Devuelve NaN fuera del rango "seguro" (rango de mu_exp +/- tolerancia).
%
% Usage:
%   pi_i = interpolate_pi_from_mu_log(mu_exp, pi_exp, mu_query);
%   pi_i = interpolate_pi_from_mu_log(mu_exp, pi_exp, mu_query, 'TolRel', 0.02);
%   pi_i = interpolate_pi_from_mu_log(mu_exp, pi_exp, mu_query, 'TolAbs', 1e-4);
%
% Notes:
% - Interpola en escala log(pi) si pi_exp>0.
% - Filtra NaNs y valores no positivos.
%
% David/ChatGPT

p = inputParser;
p.addParameter('TolRel', 0.00, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.addParameter('TolAbs', 0.00, @(x) isnumeric(x) && isscalar(x) && x>=0);
p.parse(varargin{:});
tolRel = p.Results.TolRel;
tolAbs = p.Results.TolAbs;

mu_exp   = mu_exp(:);
pi_exp   = pi_exp(:);
mu_query = mu_query(:);

% Filtrado de datos válidos
valid = isfinite(mu_exp) & isfinite(pi_exp) & (pi_exp > 0);
mu_e = mu_exp(valid);
pi_e = pi_exp(valid);

pi_interp = nan(size(mu_query));

if numel(mu_e) < 2
    return; % no hay soporte para interpolar
end

% Ordenar por mu
[mu_e, idx] = sort(mu_e);
pi_e = pi_e(idx);

% Rango seguro (sin extrapolación), con tolerancia opcional
mu_min = min(mu_e);
mu_max = max(mu_e);
mu_span = mu_max - mu_min;

tol = max(tolAbs, tolRel * mu_span);
mu_lo = mu_min - tol;
mu_hi = mu_max + tol;

inRange = isfinite(mu_query) & (mu_query >= mu_lo) & (mu_query <= mu_hi);

% Interpolación en log(pi)
logpi = log(pi_e);

% Importante: 'linear' + 'extrap' = NaN evita extrapolación dura
logpi_q = interp1(mu_e, logpi, mu_query(inRange), 'linear', nan);

pi_interp(inRange) = exp(logpi_q);

end

function [xu, yu, zu] = average_duplicate_xy(x, y, z)
%AVERAGE_DUPLICATE_XY Average duplicate scattered-interpolant support points.
% Prevents MATLAB duplicate-points warnings and makes averaging explicit.
    x = x(:); y = y(:); z = z(:);
    good = isfinite(x) & isfinite(y) & isfinite(z);
    x = x(good); y = y(good); z = z(good);
    key = [round(x,12), round(y,12)];
    [~,~,ic] = unique(key, 'rows');
    xu = accumarray(ic, x, [], @mean);
    yu = accumarray(ic, y, [], @mean);
    zu = accumarray(ic, z, [], @mean);
end

function cfg = rbs_label_config()
%RBS_LABEL_CONFIG Manual offsets to avoid label overlap in log-y plots.
% Labels are deliberately placed at slightly different x positions and
% multiplicative y offsets. No white label box is used, to keep the figure
% visually cleaner after manuscript scaling.
    cfg = struct();
    cfg.B0030  = struct('x',0.0100,'yfac',1.18);
    cfg.B0034  = struct('x',0.0100,'yfac',0.86);
    cfg.B0032  = struct('x',0.01035,'yfac',1.28);
    cfg.J61100 = struct('x',0.01035,'yfac',0.72);
    cfg.J61101 = struct('x',0.00985,'yfac',1.58);
end

function cmap = make_light_iic_colormap(n)
%MAKE_LIGHT_IIC_COLORMAP Publication colormap for IIC projection figures.
% Intermediate contrast: brighter than original dark map, but less washed
% out than v3, preserving readability of black iso-kappa lines and labels.
    if nargin < 1
        n = 256;
    end
    anchors = [ ...
        0.78 0.86 1.00;  % pale blue, not white
        0.42 0.62 1.00;  % medium blue
        0.20 0.78 0.78;  % teal
        0.58 0.84 0.36;  % green
        0.98 0.80 0.18;  % warm yellow
        1.00 0.66 0.05]; % orange-yellow high end
    x  = linspace(0,1,size(anchors,1));
    xi = linspace(0,1,n);
    cmap = interp1(x, anchors, xi, 'pchip');
    cmap = max(0,min(1,cmap));
end


function set_axis_publication(ax,font_size)
    ax.FontName = 'Arial';
    ax.FontSize = font_size;
    ax.LineWidth = 1.2;
    ax.TickDir = 'out';
    ax.TickLabelInterpreter = 'latex';
    ax.XLabel.FontSize = font_size + 2;
    ax.YLabel.FontSize = font_size + 2;
end

function save_figure_publication(fig,outdir,basename,resolution,export_pause)
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end

    drawnow;
    pause(export_pause);

    pngfile = fullfile(outdir,[basename '.png']);
    svgfile = fullfile(outdir,[basename '.svg']);
    pdffile = fullfile(outdir,[basename '.pdf']);

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
        exportgraphics(fig,svgfile, ...
            'ContentType','vector', ...
            'BackgroundColor','white');
    catch
        print(fig,svgfile,'-dsvg');
    end

    drawnow;
    pause(export_pause);

    try
        exportgraphics(fig,pdffile, ...
            'ContentType','vector', ...
            'BackgroundColor','white');
    catch
        print(fig,pdffile,'-dpdf','-bestfit');
    end

    fprintf('Saved figures:\n  %s\n  %s\n  %s\n', ...
        pngfile, svgfile, pdffile);
end

function apply_colormap_safe(ax,cmap)
%APPLY_COLORMAP_SAFE Apply colormap robustly to scalar axes/figure handles.
    if isempty(ax) || ~isscalar(ax) || ~isgraphics(ax)
        ax = gca;
    end
    try
        colormap(ax,cmap);
    catch
        % Fallback for MATLAB versions/contexts where axes colormap fails.
        fig = ancestor(ax,'figure');
        if isempty(fig) || ~isscalar(fig) || ~isgraphics(fig)
            colormap(cmap);
        else
            colormap(fig,cmap);
        end
    end
end

function rgb = hex2rgb_local(hex)
    if isstring(hex)
        hex = char(hex);
    end
    if iscell(hex)
        hex = char(string(hex));
    end
    if ischar(hex) && size(hex,1) > 1
        rgb = zeros(size(hex,1),3);
        for i = 1:size(hex,1)
            rgb(i,:) = hex2rgb_local(strtrim(hex(i,:)));
        end
        return;
    end
    h = char(hex);
    h = strtrim(h);
    if startsWith(h,'#')
        h = h(2:end);
    end
    if numel(h) ~= 6
        error('HEX color must have 6 characters.');
    end
    rgb = [hex2dec(h(1:2)), hex2dec(h(3:4)), hex2dec(h(5:6))]/255;
end

function C = rbs_color_from_name(nm,RBS_colors)
    nm = string(nm);
    switch nm
        case "B0030"
            idx = 1;
        case "B0032"
            idx = 2;
        case "B0034"
            idx = 3;
        case "J61100"
            idx = 4;
        case "J61101"
            idx = 5;
        otherwise
            idx = 1;
    end
    C = RBS_colors(idx,:);
end
