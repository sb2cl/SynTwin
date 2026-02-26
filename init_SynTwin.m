function ROOT = init_SynTwin(varargin)
% init_SynTwin  Initializes the SynTwin project (portable paths)
%
% Usage:
%   ROOT = init_SynTwin()
%   ROOT = init_SynTwin('experimental',true,'meigo',true,'bads',false,'all',false)
%
% Default behavior:
%   - Always adds:
%         Scripts_base
%         Generate_HEM/HEM_Surrogate
%   - Does NOT add Experimental_Data, MEIGO, or BADs unless requested
%   - Does NOT add the entire repository unless 'all' = true
%
% Optional flags:
%   'experimental' : adds Experimental_Data folder
%   'meigo'        : adds Matlab/MEIGO (including subfolders)
%   'bads'         : adds Matlab/bads-master (including subfolders)
%   'all'          : adds the entire project tree using genpath(ROOT)

p = inputParser;
p.addParameter('experimental', false, @(x)islogical(x) || isnumeric(x));
p.addParameter('meigo',       false, @(x)islogical(x) || isnumeric(x));
p.addParameter('bads',        false, @(x)islogical(x) || isnumeric(x));
p.addParameter('all',         false, @(x)islogical(x) || isnumeric(x));
p.parse(varargin{:});
opt = p.Results;

ROOT = SynTwin_root();   % Detect and cache project root

% --- Always-added base folders ---
add_if_exists(fullfile(ROOT,'Scripts_base'));
add_if_exists(fullfile(ROOT,'Generate_HEM','HEM_Surrogate'));

% --- Option: add entire project ---
if opt.all
    addpath(genpath(ROOT));
    return
end

% --- Optional folders ---
if opt.experimental
    add_if_exists(fullfile(ROOT,'Experimental_Data'));
end

if opt.meigo
    add_if_exists(fullfile(ROOT,'Matlab','MEIGO'));
    addpath(genpath(fullfile(ROOT,'Matlab','MEIGO')));
end

if opt.bads
    add_if_exists(fullfile(ROOT,'Matlab','bads-master'));
    addpath(genpath(fullfile(ROOT,'Matlab','bads-master')));
end

end


% ===================== Local helper =====================

function add_if_exists(folder_path)
% Adds folder to path if it exists
if exist(folder_path,'dir')
    addpath(folder_path);
else
    warning('init_SynTwin:MissingDir', ...
            'Directory not found: %s', folder_path);
end
end

