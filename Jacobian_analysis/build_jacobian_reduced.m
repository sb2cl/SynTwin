function [J_absolute, J_relative, paramNames, rowInfo] = build_jacobian_reduced(libraryIDs, TUs, ORIs, Promoters, RBSs, varargin)
%BUILD_JACOBIAN 
%   Construye un Jacobiano con la estructura sparse real de la librería.
%
%   [J_absolute, J_relative, paramNames, rowInfo] = build_jacobian_reduced(libraryIDs, TUs, ORIs, Promoters, RBSs, ...)
%
%   Name-value pairs relevantes:
%     'IncludePromoters' : cell array de promotores a tratar como parámetros
%     'IncludeRBS'       : cell array de RBS a tratar como parámetros
%     'IncludeORIs'      : cell array de ORIs a tratar como parámetros
%     'Sensitivities'    : struct array con campos
%                          .ID, .dPi_dN, .dPi_domega, .dPi_dthetaK

% ------------------- parse de opciones --------------------
p = inputParser;
addParameter(p,'IncludeORIs',     unique(ORIs,'stable'),      @(x)iscell(x));
addParameter(p,'IncludePromoters',unique(Promoters,'stable'),@(x)iscell(x));
addParameter(p,'IncludeRBS',      unique(RBSs,'stable'),      @(x)iscell(x));
addParameter(p,'Sensitivities', [], @(s)isstruct(s) || isempty(s));
parse(p,varargin{:});
opts = p.Results;

sensStruct = opts.Sensitivities;

% ------------------- seleccionar TUs ----------------------
[~, idxRows] = ismember(libraryIDs, TUs);
idxRows(idxRows==0) = [];  % por si acaso

oris_all  = unique(ORIs(idxRows),      'stable');
proms_all = unique(Promoters(idxRows), 'stable');
rbss_all  = unique(RBSs(idxRows),      'stable');

oris  = intersect(oris_all,  opts.IncludeORIs,      'stable');
proms = intersect(proms_all, opts.IncludePromoters, 'stable');
rbss  = intersect(rbss_all,  opts.IncludeRBS,       'stable');

% ------------------- lista de parámetros ----------------------
paramNames = {};

% Bloque 1: N_ori
for i=1:numel(oris)
    paramNames{end+1} = ['N_' oris{i}]; %#ok<AGROW>
end
idxEnd_N = numel(paramNames);

% Bloque 2: omega_prom
for i=1:numel(proms)
    paramNames{end+1} = ['omega_' proms{i}]; %#ok<AGROW>
end
idxEnd_omega = numel(paramNames);

% Bloque 3: thetaK_rbs
for i=1:numel(rbss)
    paramNames{end+1} = ['thetaK_' rbss{i}]; %#ok<AGROW>
end

nParams = numel(paramNames);

% ------------------- precontar nº filas ----------------------
nRows = 0;
for k = 1:numel(idxRows)
    tuID = libraryIDs{k};
    sens_k = find_sens_for_TU(sensStruct, tuID);
    nRows_k = local_sens_length(sens_k);
    nRows = nRows + nRows_k;
end


J_absolute = zeros(nRows, nParams);
J_relative = zeros(nRows, nParams);
rowInfo(nRows,1) = struct('ID','','localRow',0);

% ------------------- rellenar Jacobianos absoluto y relativo ----------------------
row = 0;

for kk = 1:numel(idxRows)
    iTU = idxRows(kk);
    tuID = TUs{iTU};

    ori  = ORIs{iTU};
    prom = Promoters{iTU};
    rbs  = RBSs{iTU};
    
    col_N      = find(strcmp(paramNames, ['N_' ori]));
    col_omega  = find(strcmp(paramNames, ['omega_' prom]));
    col_thetaK = find(strcmp(paramNames, ['thetaK_' rbs]));
    
    sens_k = find_sens_for_TU(sensStruct, tuID); %Includes absolute and relative sensitivities
    nRows_k = local_sens_length(sens_k);

    for r = 1:nRows_k
        row = row + 1;
         if ~isempty(col_N)      && isfield(sens_k,'dPi_dN_absolute')      && ~isempty(sens_k.dPi_dN_absolute)
            J_absolute(row,col_N)      = sens_k.dPi_dN_absolute(r);
        end
        if ~isempty(col_omega)  && isfield(sens_k,'dPi_domega_absolute')  && ~isempty(sens_k.dPi_domega_absolute)
            J_absolute(row,col_omega)  = sens_k.dPi_domega_absolute(r);
        end
        if ~isempty(col_thetaK) && isfield(sens_k,'dPi_dthetaK_absolute') && ~isempty(sens_k.dPi_dthetaK_absolute)
            J_absolute(row,col_thetaK) = sens_k.dPi_dthetaK_absolute(r);
        end
        % Relative next
        if ~isempty(col_N)      && isfield(sens_k,'dPi_dN_relative')      && ~isempty(sens_k.dPi_dN_relative)
            J_relative(row,col_N)      = sens_k.dPi_dN_relative(r);
        end
        if ~isempty(col_omega)  && isfield(sens_k,'dPi_domega_relative')  && ~isempty(sens_k.dPi_domega_relative)
            J_relative(row,col_omega)  = sens_k.dPi_domega_relative(r);
        end
        if ~isempty(col_thetaK) && isfield(sens_k,'dPi_dthetaK_relative') && ~isempty(sens_k.dPi_dthetaK_relative)
            J_relative(row,col_thetaK) = sens_k.dPi_dthetaK_relative(r);
        end
        rowInfo(row).ID       = tuID;
        rowInfo(row).localRow = r;
    end
   
end

end

% ---------- helpers internos ----------

function sens_k = find_sens_for_TU(sensStruct, tuID)
    sens_k = struct();
    if isempty(sensStruct)
        return;
    end
    ids = {sensStruct.ID};
    j = find(strcmp(ids, tuID), 1);
    if ~isempty(j)
        sens_k = sensStruct(j);
    else
        sens_k = struct();
    end
end

function nRows_k = local_sens_length(sens_k)
    if isempty(sens_k)
        nRows_k = 0;
        return;
    end
    fields = {'dPi_dN_absolute','dPi_domega_absolute','dPi_dthetaK_absolute'};
    nRows_k = 0;
    for f = 1:numel(fields)
        fn = fields{f};
        if isfield(sens_k,fn) && ~isempty(sens_k.(fn))
            nRows_k = max(nRows_k, numel(sens_k.(fn)));
        end
    end
end
