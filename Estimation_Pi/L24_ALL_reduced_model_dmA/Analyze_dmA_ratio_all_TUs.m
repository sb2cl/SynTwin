%% Analyze_dmA_ratio_all_TUs
% Analyze the turnover-to-resource ratio across Lib24 and Lib6 constructs.
%
% Computes
%
%   d_mA*(rho^0 + f_s(mu))/(kappa^0*varphi(mu))
%
% on the stored prediction grid and at construct-specific experimental growth
% rates. The script compares Lib24 and Lib6 for one selected d_mA value.
%
% IMPORTANT
%   This script requires both the Lib24 and Lib6 compiled result tensors for
%   the selected d_mA. The Lib6 tensor is produced by the corresponding Lib6
%   d_mA workflow and is not generated in this folder.
%
% Outputs:
%   Analysis_results/  MAT and CSV summaries.
%   Figures/           Per-library PNG and SVG plots.
%
% USAGE
%   Analyze_dmA_ratio_all_TUs
%
% See README.md for file placement and dependencies.

clear; clc;

%% =========================
% User options
% ==========================
dmA = 0.20;
rho0 = 0.02;
Use_mean = 'Wells';

SCRIPT_DIR = fileparts(mfilename('fullpath'));
RESULTS_DIR = SCRIPT_DIR;
ANALYSIS_DIR = fullfile(SCRIPT_DIR,'Analysis_results');
FIGURES_DIR = fullfile(SCRIPT_DIR,'Figures');
if ~exist(ANALYSIS_DIR,'dir'), mkdir(ANALYSIS_DIR); end
if ~exist(FIGURES_DIR,'dir'), mkdir(FIGURES_DIR); end

lineWidthCurve = 1.8;
lineWidthExp = 1.0;
markerSizeExp = 22;
baseFont = 10;
titleFont = 11;

%% =========================
% File names
% ==========================
dmA_tag = extractAfter(num2str(dmA), '.');
if isempty(dmA_tag)
    dmA_tag = '0';
end

% DIR_LIB6 = SCRIPT_DIR;
% DIR_LIB24 = resolve_existing_dir({ ...
%     fullfile(fileparts(SCRIPT_DIR), 'L24_ALL_reduced_model_dmA'), ...
%     fullfile(SCRIPT_DIR, '..', 'L24_ALL_reduced_model_dmA')});


file_lib24 = fullfile(RESULTS_DIR, ...
    "Results_Tensor_Lib24_ALL_reduced_model_" + Use_mean + "_dma_" + dmA_tag + ".mat");
file_lib6 = fullfile(RESULTS_DIR, ...
    "Results_Tensor_Lib6_ALL_reduced_model_" + Use_mean + "_dma_" + dmA_tag + ".mat");

if ~isfile(file_lib24)
    error('Lib24 results file not found: %s', file_lib24);
end
if ~isfile(file_lib6)
    error('Lib6 results file not found: %s', file_lib6);
end

S24 = load(file_lib24);
S6  = load(file_lib6);

if ~isfield(S24, 'Results_Tensor_Lib24_ALL_reduced')
    error('Variable Results_Tensor_Lib24_ALL_reduced not found in %s', file_lib24);
end
if ~isfield(S6, 'Results_Tensor_Lib6_ALL_reduced')
    error('Variable Results_Tensor_Lib6_ALL_reduced not found in %s', file_lib6);
end

Results_Tensor_Lib24_ALL_reduced = S24.Results_Tensor_Lib24_ALL_reduced;
Results_Tensor_Lib6_ALL_reduced  = S6.Results_Tensor_Lib6_ALL_reduced;

%% =========================
% Analyze Lib24 and Lib6
% ==========================
Lib24 = analyze_library_tensor(Results_Tensor_Lib24_ALL_reduced, dmA, rho0, 'Lib24');
Lib6  = analyze_library_tensor(Results_Tensor_Lib6_ALL_reduced,  dmA, rho0, 'Lib6');

%% =========================
% Build global output structure
% ==========================
RatioAnalysis = struct();
RatioAnalysis.meta.dmA = dmA;
RatioAnalysis.meta.rho0 = rho0;
RatioAnalysis.meta.Use_mean = Use_mean;
RatioAnalysis.meta.file_lib24 = file_lib24;
RatioAnalysis.meta.file_lib6 = file_lib6;
RatioAnalysis.Lib24 = Lib24;
RatioAnalysis.Lib6 = Lib6;

%% =========================
% Export summary tables
% ==========================
T24 = build_summary_table(Lib24);
T6  = build_summary_table(Lib6);

csv24 = fullfile(ANALYSIS_DIR, sprintf('RatioSummary_Lib24_dmA_%s.csv', strrep(sprintf('%.2f', dmA), '.', 'p')));
csv6  = fullfile(ANALYSIS_DIR, sprintf('RatioSummary_Lib6_dmA_%s.csv',  strrep(sprintf('%.2f', dmA), '.', 'p')));
writetable(T24, csv24);
writetable(T6,  csv6);

%% =========================
% Save MAT structure
% ==========================
mat_out = fullfile(ANALYSIS_DIR, sprintf('RatioAnalysis_dmA_%s.mat', strrep(sprintf('%.2f', dmA), '.', 'p')));
save(mat_out, 'RatioAnalysis', '-mat');

%% =========================
% Plots
% ==========================
fig24 = plot_library_ratios(Lib24, 'Lib24', baseFont, titleFont, lineWidthCurve, lineWidthExp, markerSizeExp);
export_figure(fig24, fullfile(FIGURES_DIR, sprintf('Fig_Ratio_AllTUs_Lib24_dmA_%s', strrep(sprintf('%.2f', dmA), '.', 'p'))));

fig6 = plot_library_ratios(Lib6, 'Lib6', baseFont, titleFont, lineWidthCurve, lineWidthExp, markerSizeExp);
export_figure(fig6, fullfile(FIGURES_DIR, sprintf('Fig_Ratio_AllTUs_Lib6_dmA_%s', strrep(sprintf('%.2f', dmA), '.', 'p'))));

fprintf('\nSaved:\n');
fprintf('  %s\n', mat_out);
fprintf('  %s\n', csv24);
fprintf('  %s\n', csv6);

%% =========================
% Local functions
% ==========================
function Lib = analyze_library_tensor(R, dmA, rho0, libName)

    sz = size(R);
    nd = ndims(R);

    entries = {};
    subs = {};
    for idx = 1:numel(R)
        if ~isempty(R{idx})
            entries{end+1} = R{idx}; %#ok<AGROW>
            s = cell(1, nd);
            [s{:}] = ind2sub(sz, idx);
            subs{end+1} = cell2mat(s); %#ok<AGROW>
        end
    end

    nTU = numel(entries);
    TU = repmat(struct(), nTU, 1);

    for t = 1:nTU
        X = entries{t};
        subidx = subs{t};

        mu_grid = get_field_flexible(X.Synthesis_predictions, {'Mu_values','mu_values','mu'});
        fs_grid = get_field_flexible(X.Synthesis_predictions, {'fs_pred_values','f_pred_values','fs_values'});
        varphi_grid = get_field_flexible(X.Synthesis_predictions, {'varphi_pred_values','phi_pred_values','varphi_values'});

        mu_grid = mu_grid(:);
        fs_grid = fs_grid(:);
        varphi_grid = varphi_grid(:);

        kappa = get_kappa_mean(X);
        ratio_grid = dmA .* (fs_grid + rho0) ./ (kappa .* varphi_grid);

        mu_exp = get_field_flexible(X, {'Mu_mumax_pmax_global_mean'});
        mu_exp = mu_exp(:);
        ratio_exp = interp1(mu_grid, ratio_grid, mu_exp, 'linear', 'extrap');

        TU(t).library = libName;
        TU(t).subindices = subidx;
        TU(t).TU_Name = get_optional_text_field(X, {'TU_Name'});
        TU(t).TU_Bioparts = get_optional_text_field(X, {'TU_Bioparts'});
        TU(t).TU_Ori = get_optional_text_field(X, {'TU_Ori'});
        TU(t).TU_Promoter = get_optional_text_field(X, {'TU_Promoter'});
        TU(t).TU_RBS = get_optional_text_field(X, {'TU_RBS'});

        TU(t).dmA = dmA;
        TU(t).rho0 = rho0;
        TU(t).kappa_mean = kappa;

        TU(t).mu_grid = mu_grid;
        TU(t).fs_grid = fs_grid;
        TU(t).varphi_grid = varphi_grid;
        TU(t).ratio_grid = ratio_grid;

        TU(t).mu_exp = mu_exp;
        TU(t).ratio_exp = ratio_exp;

        TU(t).ratio_grid_min = min(ratio_grid);
        TU(t).ratio_grid_max = max(ratio_grid);
        TU(t).ratio_exp_min = min(ratio_exp);
        TU(t).ratio_exp_max = max(ratio_exp);
        TU(t).ratio_exp_median = median(ratio_exp);
        TU(t).ratio_exp_mean = mean(ratio_exp);
        TU(t).fraction_exp_ratio_gt_1 = mean(ratio_exp > 1);
        TU(t).fraction_grid_ratio_gt_1 = mean(ratio_grid > 1);

        if all(ratio_exp > 1)
            TU(t).regime_exp = "turnover-dominated";
        elseif all(ratio_exp < 1)
            TU(t).regime_exp = "resource-dominated";
        else
            TU(t).regime_exp = "mixed";
        end
    end

    Lib = struct();
    Lib.name = libName;
    Lib.nTU = nTU;
    Lib.TU = TU;
end

function kappa = get_kappa_mean(X)
    candidates = {'RBS_k0_sigma0_mean','RBS_k0_mean','kappa_mean'};
    kappa = [];
    for i = 1:numel(candidates)
        if isfield(X, candidates{i})
            kappa = X.(candidates{i});
            break;
        end
    end
    if isempty(kappa)
        error('Could not find RBS mean parameter field.');
    end
    kappa = kappa(1);
end

function val = get_field_flexible(S, candidates)
    val = [];
    for i = 1:numel(candidates)
        if isfield(S, candidates{i})
            val = S.(candidates{i});
            return;
        end
    end
    error('Could not find any of these fields: %s', strjoin(candidates, ', '));
end

function txt = get_optional_text_field(S, candidates)
    txt = "";
    for i = 1:numel(candidates)
        if isfield(S, candidates{i})
            v = S.(candidates{i});
            if isstring(v)
                txt = v;
            elseif ischar(v)
                txt = string(v);
            elseif iscell(v)
                try
                    txt = string(strjoin(string(v), ', '));
                catch
                    txt = "";
                end
            else
                try
                    txt = string(v);
                catch
                    txt = "";
                end
            end
            if ~isempty(txt)
                txt = txt(1);
            end
            return;
        end
    end
end

function T = build_summary_table(Lib)
    nTU = Lib.nTU;

    Library   = strings(nTU,1);
    TU_Name   = strings(nTU,1);
    Ori       = strings(nTU,1);
    Promoter  = strings(nTU,1);
    RBS       = strings(nTU,1);
    KappaMean = nan(nTU,1);

    MuMin     = nan(nTU,1);
    MuMax     = nan(nTU,1);
    RatioExpMedian = nan(nTU,1);
    RatioExpMean   = nan(nTU,1);
    RatioExpMin    = nan(nTU,1);
    RatioExpMax    = nan(nTU,1);
    FracExpGT1     = nan(nTU,1);
    FracGridGT1    = nan(nTU,1);
    Regime         = strings(nTU,1);

    for t = 1:nTU
        X = Lib.TU(t);
        Library(t) = string(X.library);
        TU_Name(t) = string(X.TU_Name);
        Ori(t) = string(X.TU_Ori);
        Promoter(t) = string(X.TU_Promoter);
        RBS(t) = string(X.TU_RBS);
        KappaMean(t) = X.kappa_mean;

        MuMin(t) = min(X.mu_exp);
        MuMax(t) = max(X.mu_exp);
        RatioExpMedian(t) = X.ratio_exp_median;
        RatioExpMean(t) = X.ratio_exp_mean;
        RatioExpMin(t) = X.ratio_exp_min;
        RatioExpMax(t) = X.ratio_exp_max;
        FracExpGT1(t) = X.fraction_exp_ratio_gt_1;
        FracGridGT1(t) = X.fraction_grid_ratio_gt_1;
        Regime(t) = string(X.regime_exp);
    end

    T = table(Library, TU_Name, Ori, Promoter, RBS, KappaMean, ...
        MuMin, MuMax, RatioExpMedian, RatioExpMean, RatioExpMin, RatioExpMax, ...
        FracExpGT1, FracGridGT1, Regime);
end

function fig = plot_library_ratios(Lib, libName, baseFont, titleFont, lineWidthCurve, lineWidthExp, markerSizeExp)
    nTU = Lib.nTU;

    if strcmpi(libName, 'Lib24')
        nrows = 4; ncols = 6;
        figPos = [2 2 36 21];
    else
        nrows = 2; ncols = 3;
        figPos = [2 2 24 13];
    end

    fig = figure('Units','centimeters','Position',figPos,'Color','w');
    tiledlayout(nrows, ncols, 'TileSpacing','compact','Padding','compact');

    for t = 1:nTU
        X = Lib.TU(t);
        nexttile; hold on;

        plot(X.mu_grid, X.ratio_grid, '-', 'LineWidth', lineWidthCurve, 'Color', [0.10 0.45 0.80]);
        scatter(X.mu_exp, X.ratio_exp, markerSizeExp, 'o', ...
            'MarkerFaceColor', [0.10 0.45 0.80], ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.5);
        plot([min(X.mu_exp), max(X.mu_exp)], [1 1], '--', 'Color', [0.2 0.2 0.2], 'LineWidth', lineWidthExp);

        if max(X.mu_exp) > min(X.mu_exp)
            xpad = 0.07 * (max(X.mu_exp) - min(X.mu_exp));
        else
            xpad = 0.001;
        end
        xlim([min(X.mu_exp)-xpad, max(X.mu_exp)+xpad]);

        title(compose_title(X), 'Interpreter', 'none', 'FontSize', titleFont);
        xlabel('\mu (min^{-1})');
        ylabel('ratio');
        grid on; box on;
        set(gca, 'FontSize', baseFont, 'LineWidth', 1.0, 'Layer', 'top');
    end

    sgtitle(sprintf('%s: ratio = d_{mA}(\\rho^0 + f_s)/(\\kappa\\varphi) at TU-specific experimental \\mu ranges', libName), ...
        'FontSize', titleFont+1, 'FontWeight', 'bold');
end

function txt = compose_title(X)
    name = string(X.TU_Name);
    if strlength(name) == 0
        parts = [string(X.TU_Ori), string(X.TU_Promoter), string(X.TU_RBS)];
        parts = parts(parts ~= "");
        txt = strjoin(parts, ' | ');
    else
        txt = name;
    end
end

function export_figure(fig_handle, file_base)
    set(fig_handle, 'PaperPositionMode', 'auto');
    drawnow;
    exportgraphics(fig_handle, [file_base '.png'], 'Resolution', 300);
    exportgraphics(fig_handle, [file_base '.svg'], 'ContentType', 'vector');
end