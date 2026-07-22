%% Show_Data_Exp_lib35
% Reproduce the experimental-data plots generated from the 35-construct library.
%
% Publication:
%   "Host-aware Identification of Intrinsic Gene Expression Biopart
%   Parameters using Combinatorial Libraries"
%
% Main-text output:
%   - Figure 2a: ordered construct synthesis rate, maximum growth rate,
%     and estimated construct-induced host burden.
%
% Additional outputs:
%   - Supplementary ordered synthesis/growth and synthesis/burden plots;
%   - Fano-factor and coefficient-of-variation summary;
%   - tiled Pi_p(mu) trajectories across all expression constructs.
%
% Inputs:
%   Experimental_data/ExpData_Tensor_lib35.mat
%
% Shared code:
%   Scripts_base/Plot_ExpData_Tensor_v10.m
%
% Requirements:
%   - MATLAB R2020a or later (uses tiledlayout and exportgraphics)
%   - No non-core MATLAB toolboxes are required by this script
%
% Usage:
%   Run this script from any working directory. Output files are written to
%   ./Figures/.
%
% Data provenance:
%   ExpData_Tensor_lib35.mat contains the processed experimental data and
%   summary statistics used to generate the plotted values. The full SynTwin
%   codebase and archived release are available at:
%   https://github.com/sb2cl/SynTwin
%   https://doi.org/10.5281/zenodo.18787107
%
% Notes:
%   The published figures were generated with the MATLAB implementation.
%   A Python reimplementation (Py_SynTwin) is under development.
%


clearvars;
close all;
clc;
warning('off','MATLAB:Axes:NegativeDataInLogAxis');
warning('off','MATLAB:plot:IgnoreImaginaryXYPart');
warning('off','MATLAB:print:ContentTypeImageSuggested');

% --- Portable project initialization (no absolute paths) ---
ROOT = init_SynTwin('experimental',true);
% Ensure this script folder is on the path (for local functions in this folder)
SCRIPT_DIR = fileparts(mfilename('fullpath'));
addpath(SCRIPT_DIR);


% Load the processed experimental tensor used for Figure 2a and the
% associated Supplementary figures.
% Getting the experimental data of Lib35: 
S = load(SynTwin_path('Experimental_Data','ExpData_Tensor_lib35.mat'));        % loads data of Lib30: ExpData_Tensor_lib30_micro

if ~isfield(S, 'Data_Exp_Tensor_lib35')
    error('Show_Data_Exp_lib35_v2:MissingVariable', ...
        ['The file ExpData_Tensor_lib35.mat must contain the variable ' ...
         'Data_Exp_Tensor_lib35.']);
end

outDir = fullfile(SCRIPT_DIR, 'Figures');

Plot_ExpData_Tensor_v10(S.Data_Exp_Tensor_lib35, 0, 0, 0, 0, ...
    'OutputDir', outDir, ...
    'SaveFigures', true, ...
    'MakeMainVariant', true, ...
    'LabelMode', 'ec_code', ...
    'UseVerticalLabels', true, ...
    'SupplementFontSize', 18, ...
    'MainFontSize', 16, ...
    'FigurePrefix', 'Lib35_Experimental');
