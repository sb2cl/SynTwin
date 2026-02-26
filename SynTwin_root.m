function ROOT = SynTwin_root()
% SynTwin_root  Returns absolute path to the "SynTwin" root folder.
%          Uses persistent caching for efficiency.
%
% The function walks upward from the caller script location
% until it finds a folder named "SynTwin".

persistent ROOT_CACHED

% Return cached root if valid
if ~isempty(ROOT_CACHED) && exist(ROOT_CACHED,'dir')
    ROOT = ROOT_CACHED;
    return
end

start_folder = guess_start_folder();
folder = start_folder;

while true
    [parent, name] = fileparts(folder);

    % Stop when we find folder named "SynTwin"
    if strcmp(name,'SynTwin')
        ROOT_CACHED = folder;
        ROOT = folder;
        return
    end

    % Stop if we reached filesystem root
    if strcmp(parent, folder)
        error(['Could not locate "SynTwin" root folder starting from: ' ...
               start_folder newline ...
               'Make sure the executed script is inside the project tree.']);
    end

    folder = parent;
end

end


% ===================== Helpers =====================

function start_folder = guess_start_folder()
% Attempts to detect the caller script location.
% Falls back to current working directory if needed.

try
    stack = dbstack('-completenames');
    if numel(stack) >= 2 && isfield(stack(2),'file') && ~isempty(stack(2).file)
        start_folder = fileparts(stack(2).file);
    else
        start_folder = pwd;
    end
catch
    start_folder = pwd;
end

end
