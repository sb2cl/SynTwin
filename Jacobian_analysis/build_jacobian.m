function [J, paramNames, rowInfo] = build_jacobian(libraryIDs, TUs, ORIs, Promoters, RBSs, varargin)
%BUILD_JACOBIAN 
%   Construye un Jacobiano con la estructura sparse real de la librería.
%   Puede operar en modo dummy o usar sensibilidades experimentales.
%
%   [J, paramNames, rowInfo] = build_jacobian(libraryIDs, TUs, ORIs, Promoters, RBSs, ...)
%
%   Name-value pairs relevantes:
%     'IncludePromoters' : cell array de promotores a tratar como parámetros
%     'IncludeRBS'       : cell array de RBS a tratar como parámetros
%     'IncludeORIs'      : cell array de ORIs a tratar como parámetros
%     'nMuDummy'         : nº de puntos de µ por TU (modo dummy)
%     'colinearRBS'      : true/false (modo dummy)
%     'Sensitivities'    : struct array con campos
%                          .ID, .dPi_domega, .dPi_dthetaK, .dPi_dthetaS, .dPi_dN

% ------------------- parse de opciones --------------------
p = inputParser;
addParameter(p,'IncludePromoters',unique(Promoters,'stable'),@(x)iscell(x));
addParameter(p,'IncludeRBS',      unique(RBSs,'stable'),      @(x)iscell(x));
addParameter(p,'IncludeORIs',     unique(ORIs,'stable'),      @(x)iscell(x));
addParameter(p,'nMuDummy',        8, @(x)isnumeric(x) && isscalar(x));
addParameter(p,'colinearRBS', true, @(x)islogical(x) && isscalar(x));
addParameter(p,'Sensitivities', [], @(s)isstruct(s) || isempty(s));
parse(p,varargin{:});
opts = p.Results;

useExperimental = ~isempty(opts.Sensitivities);
sensStruct = opts.Sensitivities;

% ------------------- seleccionar TUs ----------------------
[~, idxRows] = ismember(libraryIDs, TUs);
idxRows(idxRows==0) = [];  % por si acaso

proms_all = unique(Promoters(idxRows), 'stable');
rbss_all  = unique(RBSs(idxRows),      'stable');
oris_all  = unique(ORIs(idxRows),      'stable');

proms = intersect(proms_all, opts.IncludePromoters, 'stable');
rbss  = intersect(rbss_all,  opts.IncludeRBS,       'stable');
oris  = intersect(oris_all,  opts.IncludeORIs,      'stable');

% ------------------- lista de parámetros ----------------------
paramNames = {};

% Bloque 1: omega_prom
for i=1:numel(proms)
    paramNames{end+1} = ['omega_' proms{i}]; %#ok<AGROW>
end
idxEnd_omega = numel(paramNames);

% Bloque 2: thetaK_rbs
for i=1:numel(rbss)
    paramNames{end+1} = ['thetaK_' rbss{i}]; %#ok<AGROW>
end
idxEnd_thetaK = numel(paramNames);

% Bloque 3: thetaS_rbs
for i=1:numel(rbss)
    paramNames{end+1} = ['thetaS_' rbss{i}]; %#ok<AGROW>
end
idxEnd_thetaS = numel(paramNames);

% Bloque 4: N_ori
for i=1:numel(oris)
    paramNames{end+1} = ['N_' oris{i}]; %#ok<AGROW>
end

nParams = numel(paramNames);

% ------------------- precontar nº filas ----------------------
if useExperimental
    nRows = 0;
    for k = 1:numel(idxRows)
        tuID = libraryIDs{k};
        sens_k = find_sens_for_TU(sensStruct, tuID);
        nRows_k = local_sens_length(sens_k);
        nRows = nRows + nRows_k;
    end
else
    nTUs_sel = numel(idxRows);
    nRows = nTUs_sel * opts.nMuDummy;
end

J = zeros(nRows, nParams);
rowInfo(nRows,1) = struct('ID','','localRow',0);

% ------------------- rellenar Jacobiano ----------------------
row = 0;

for kk = 1:numel(idxRows)
    iTU = idxRows(kk);
    tuID = TUs{iTU};

    prom = Promoters{iTU};
    rbs  = RBSs{iTU};
    ori  = ORIs{iTU};

    col_omega  = find(strcmp(paramNames, ['omega_' prom]));
    col_thetaK = find(strcmp(paramNames, ['thetaK_' rbs]));
    col_thetaS = find(strcmp(paramNames, ['thetaS_' rbs]));
    col_N      = find(strcmp(paramNames, ['N_' ori]));

    if useExperimental
        sens_k = find_sens_for_TU(sensStruct, tuID);
        nRows_k = local_sens_length(sens_k);

        for r = 1:nRows_k
            row = row + 1;
            if ~isempty(col_omega)  && isfield(sens_k,'dPi_domega')  && ~isempty(sens_k.dPi_domega)
                J(row,col_omega)  = sens_k.dPi_domega(r);
            end
            if ~isempty(col_thetaK) && isfield(sens_k,'dPi_dthetaK') && ~isempty(sens_k.dPi_dthetaK)
                J(row,col_thetaK) = sens_k.dPi_dthetaK(r);
            end
            if ~isempty(col_thetaS) && isfield(sens_k,'dPi_dthetaS') && ~isempty(sens_k.dPi_dthetaS)
                J(row,col_thetaS) = sens_k.dPi_dthetaS(r);
            end
            if ~isempty(col_N)      && isfield(sens_k,'dPi_dN')      && ~isempty(sens_k.dPi_dN)
                J(row,col_N)      = sens_k.dPi_dN(r);
            end
            rowInfo(row).ID       = tuID;
            rowInfo(row).localRow = r;
        end

    else
        % Modo dummy
        nMu = opts.nMuDummy;
        base_thetaK = randn(nMu,1);
        if opts.colinearRBS
            base_thetaS = base_thetaK + 0.01*randn(nMu,1);
        else
            base_thetaS = randn(nMu,1);
        end
        base_omega = randn(nMu,1);
        base_N     = randn(nMu,1);

        for r = 1:nMu
            row = row + 1;
            if ~isempty(col_omega),  J(row,col_omega)  = base_omega(r); end
            if ~isempty(col_thetaK), J(row,col_thetaK) = base_thetaK(r); end
            if ~isempty(col_thetaS), J(row,col_thetaS) = base_thetaS(r); end
            if ~isempty(col_N),      J(row,col_N)      = base_N(r);     end

            rowInfo(row).ID       = tuID;
            rowInfo(row).localRow = r;
        end
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
    fields = {'dPi_domega','dPi_dthetaK','dPi_dthetaS','dPi_dN'};
    nRows_k = 0;
    for f = 1:numel(fields)
        fn = fields{f};
        if isfield(sens_k,fn) && ~isempty(sens_k.(fn))
            nRows_k = max(nRows_k, numel(sens_k.(fn)));
        end
    end
end
