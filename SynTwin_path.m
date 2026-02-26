function full_path = SynTwin_path(varargin)
% SynTwin_path  Builds absolute paths inside SynTwin project.
%
% Examples:
%   p = SynTwin_path('Experimental_Data','file.mat')
%   p = SynTwin_path('Estimation_Pi','Problem1','Results','res.mat')

ROOT = SynTwin_root();
full_path = fullfile(ROOT, varargin{:});

end
