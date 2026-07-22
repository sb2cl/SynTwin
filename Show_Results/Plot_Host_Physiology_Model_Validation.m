%% Plot_Host_Physiology_Model_Validation.m
% Validate the host-physiology model against literature datasets.
%
% Associated publication:
%   "Host-aware Identification of Intrinsic Gene Expression Biopart
%    Parameters using Combinatorial Libraries"
%
% Main-text association:
%   - Figure 2b: peptide elongation rate as a function of growth rate.
%   - Figure 2c: mature active ribosomal protein-mass fraction as a
%     function of growth rate.
%
% The script does not perform parameter estimation. It evaluates the fixed
% models and parameter values reported in Supplementary Section 3.5 against
% the literature datasets distributed in Generate_HEM/Third_party_data.
%
% Models:
%   Peptide elongation rate:
%       nu_t(mu) = gamma1*mu/(gamma2 + mu)
%
%   Mature active ribosomal protein-mass fraction:
%       Phi_R(mu) = a + b*mu
%
% Fixed parameters:
%       a           = 0.0660
%       b           = 7.1502 min
%       gamma1      = 25.9083 aa/s
%       gamma2      = 0.0092 min^-1
%       Phi_t*Phi_m = 0.72
%
% Inputs:
%   Generate_HEM/Third_party_data/
%       ecoli_peptide_elongation_rates_updated.xlsx
%       ecoli_ribosomal_mass_fractions_updated.xlsx
%
% Outputs in ./Figures:
%   Figure_2b_peptide_elongation.png/.svg
%   Figure_2c_ribosomal_mass_fraction.png/.svg
%
% Requirements:
%   - MATLAB R2020a or later.
%   - No optimization toolbox is required.
%
% Usage:
%   Run Plot_Host_Physiology_Model_Validation from any working directory
%   while the script remains inside the SynTwin repository.
%
% Data attribution:
%   The literature data are derived from the Chure and Cremer Flux Parity
%   compilation and upstream sources, including Bremer and Dennis. See the
%   package-level data-attribution documentation for provenance and reuse
%   conditions.
%
% See README_Show_Results.md for the complete Show_Results workflow.

clearvars;
close all;
clc;

% Restore the user's warning configuration even if execution stops early.
warning_state = warning;
warning_cleanup = onCleanup(@() warning(warning_state)); %#ok<NASGU>

% -------------------------------------------------------------------------
% Portable package paths
% -------------------------------------------------------------------------
% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin();
% Ensure this script folder is on the path (for local functions in this folder)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);
DATA_DIR = fullfile(ROOT, 'Generate_HEM', 'Third_party_data');
OUT_DIR = fullfile(SCRIPT_DIR, 'Figures');
if ~exist(OUT_DIR, 'dir')
    mkdir(OUT_DIR);
end

% -------------------------------------------------------------------------
% Load literature datasets
% -------------------------------------------------------------------------
nu_file = find_spreadsheet(DATA_DIR, 'ecoli_peptide_elongation_rates_updated');
phi_file = find_spreadsheet(DATA_DIR, 'ecoli_ribosomal_mass_fractions_updated');

[mu_nu_h, nu_exp] = load_peptide_elongation_data(nu_file);
[mu_phi_h, phi_reported, source_group] = load_ribosomal_fraction_data(phi_file);

% Literature growth rates are stored in h^-1; models use min^-1.
mu_nu = mu_nu_h / 60;
mu_phi = mu_phi_h / 60;

% The reported ribosomal fraction is scaled by the fitted active/mature
% fraction product, as described in Supplementary Section 3.5.
phi_t_phi_m = 0.72;
phi_active = phi_reported / phi_t_phi_m;

% Remove incomplete or non-finite rows without altering the source files.
valid_nu = isfinite(mu_nu) & isfinite(nu_exp);
mu_nu = mu_nu(valid_nu);
nu_exp = nu_exp(valid_nu);

valid_phi = isfinite(mu_phi) & isfinite(phi_active);
mu_phi = mu_phi(valid_phi);
phi_active = phi_active(valid_phi);
source_group = source_group(valid_phi);

% -------------------------------------------------------------------------
% Fixed published models
% -------------------------------------------------------------------------
a = 0.0660;
b = 7.1502;
gamma1 = 25.9083;  % aa/s
gamma2 = 0.0092;   % min^-1

nu_model = @(mu) gamma1 .* mu ./ (gamma2 + mu);
phi_model = @(mu) a + b .* mu;

mu_max = max([0.04; mu_nu(:); mu_phi(:)]);
mu_grid = linspace(0, mu_max, 500).';

% -------------------------------------------------------------------------
% Publication-oriented style
% -------------------------------------------------------------------------
font_name = 'Arial';
font_size = 14;
label_size = 17;
line_width = 2.6;
marker_size = 38;
model_color = [0.70 0.18 0.15];
data_color  = [0.15 0.15 0.15];

% -------------------------------------------------------------------------
% Figure 2b: peptide elongation rate
% -------------------------------------------------------------------------
fig_b = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 12.5 10.5], 'Name', 'Figure 2b');
ax_b = axes(fig_b);
hold(ax_b, 'on');

scatter(ax_b, mu_nu, nu_exp, marker_size, data_color, 'filled', ...
    'MarkerEdgeColor', 'none', 'DisplayName', 'Literature data');
plot(ax_b, mu_grid, nu_model(mu_grid), '-', 'Color', model_color, ...
    'LineWidth', line_width, 'DisplayName', 'Fitted model');

xlabel(ax_b, '$\mu\; (\mathrm{min}^{-1})$', 'Interpreter', 'latex', ...
    'FontSize', label_size);
ylabel(ax_b, '$\nu\; (\mathrm{aa}\,\mathrm{s}^{-1})$', ...
    'Interpreter', 'latex', 'FontSize', label_size);
format_axis(ax_b, font_name, font_size);
xlim(ax_b, [0 mu_max]);
ylim(ax_b, [0 max(22, 1.06 * max([nu_exp(:); nu_model(mu_grid)]))]);
legend(ax_b, 'Location', 'southeast', 'Box', 'off', ...
    'FontName', font_name, 'FontSize', font_size - 1);

export_publication_figure(fig_b, OUT_DIR, 'Figure_2b_peptide_elongation');

% -------------------------------------------------------------------------
% Figure 2c: mature active ribosomal protein-mass fraction
% -------------------------------------------------------------------------
fig_c = figure('Color', 'w', 'Units', 'centimeters', ...
    'Position', [2 2 12.5 10.5], 'Name', 'Figure 2c');
ax_c = axes(fig_c);
hold(ax_c, 'on');

% The sixth spreadsheet column is retained as a source-group indicator in
% the distributed dataset. When it contains more than one finite value, the
% groups are shown separately; otherwise all points use one marker style.
finite_groups = unique(source_group(isfinite(source_group)));
if numel(finite_groups) > 1
    group_a = source_group == finite_groups(1);
    group_b = ~group_a;
    scatter(ax_c, mu_phi(group_a), phi_active(group_a), marker_size, ...
        data_color, 'filled', 'MarkerEdgeColor', 'none', ...
        'DisplayName', 'Literature data');
    scatter(ax_c, mu_phi(group_b), phi_active(group_b), marker_size, ...
        data_color, 'filled', 'MarkerEdgeColor', 'none', ...
        'HandleVisibility', 'off');
else
    scatter(ax_c, mu_phi, phi_active, marker_size, data_color, 'filled', ...
        'MarkerEdgeColor', 'none', 'DisplayName', 'Literature data');
end

plot(ax_c, mu_grid, phi_model(mu_grid), '-', 'Color', model_color, ...
    'LineWidth', line_width, 'DisplayName', 'Fitted model');

xlabel(ax_c, '$\mu\; (\mathrm{min}^{-1})$', 'Interpreter', 'latex', ...
    'FontSize', label_size);
ylabel(ax_c, '$\Phi_\mathrm{R}$', 'Interpreter', 'latex', 'FontSize', label_size);
format_axis(ax_c, font_name, font_size);
xlim(ax_c, [0 mu_max]);
ylim(ax_c, [0 max(0.5, 1.06 * max([phi_active(:); phi_model(mu_grid)]))]);
legend(ax_c, 'Location', 'southeast', 'Box', 'off', ...
    'FontName', font_name, 'FontSize', font_size - 1);

export_publication_figure(fig_c, OUT_DIR, 'Figure_2c_ribosomal_mass_fraction');

fprintf('Figure 2b and Figure 2c outputs written to:\n  %s\n', OUT_DIR);

% =========================================================================
% Local functions
% =========================================================================
function file_path = find_spreadsheet(data_dir, stem)
    candidates = {
        fullfile(data_dir, [stem '.xlsx'])
        fullfile(data_dir, [stem '.xls'])
        fullfile(data_dir, stem)
    };

    file_path = '';
    for i = 1:numel(candidates)
        if exist(candidates{i}, 'file') == 2
            file_path = candidates{i};
            break;
        end
    end

    if isempty(file_path)
        error('Plot_Host_Physiology_Model_Validation:MissingDataFile', ...
            ['Could not find the spreadsheet "%s" in:\n  %s\n' ...
             'Expected extension: .xlsx or .xls.'], stem, data_dir);
    end
end

function [mu_h, nu] = load_peptide_elongation_data(file_path)
    opts = spreadsheetImportOptions('NumVariables', 2);
    opts.VariableTypes = {'double', 'double'};
    opts.Sheet = 'Data_nut';
    opts.DataRange = 'A3:B73';
    T = readtable(file_path, opts, 'UseExcel', false);
    A = table2array(T);
    mu_h = A(:,1);
    nu = A(:,2);
end

function [mu_h, phi, group] = load_ribosomal_fraction_data(file_path)
    opts = spreadsheetImportOptions('NumVariables', 6);
    opts.VariableTypes = {'double', 'double', 'char', 'char', 'char', 'double'};
    opts.Sheet = 'Data_phiR';
    opts.DataRange = 'A2:F272';
    T = readtable(file_path, opts, 'UseExcel', false);
    mu_h = table2array(T(:,1));
    phi = table2array(T(:,2));
    group = table2array(T(:,6));
end

function format_axis(ax, font_name, font_size)
    grid(ax, 'on');
    box(ax, 'on');
    ax.FontName = font_name;
    ax.FontSize = font_size;
    ax.LineWidth = 1.0;
    ax.TickDir = 'out';
    ax.TickLabelInterpreter = 'latex';
    ax.GridAlpha = 0.18;
    ax.MinorGridAlpha = 0.08;
    ax.XMinorGrid = 'on';
    ax.YMinorGrid = 'on';
    ax.Layer = 'top';
end

function export_publication_figure(fig, out_dir, base_name)
    drawnow;
    png_path = fullfile(out_dir, [base_name '.png']);
    svg_path = fullfile(out_dir, [base_name '.svg']);

    try
        exportgraphics(fig, png_path, 'Resolution', 600);
    catch
        print(fig, png_path, '-dpng', '-r600');
    end

    try
        exportgraphics(fig, svg_path, 'ContentType', 'vector');
    catch
        print(fig, svg_path, '-dsvg');
    end
end
