function plot_jacobian_comparison_blocks()

%% =============================================================
% 1. Tabla completa de TUs (L24 + L30 + L5)
% =============================================================

IDs      = { ...
    'L24\_01','L24\_02','L24\_03','L24\_04', ...
    'L24\_05','L24\_06','L24\_07','L24\_08', ...
    'L24\_09','L24\_10','L24\_11','L24\_12', ...
    'L24\_13','L24\_14','L24\_15','L24\_16', ...
    'L24\_17','L24\_18','L24\_19','L24\_20', ...
    'L24\_21','L24\_22','L24\_23','L24\_24', ...
    'L30\_25','L30\_26','L30\_27','L30\_28','L30\_29','L30\_30', ...
    'L5\_1','L5\_2','L5\_4','L5\_5','L5\_3'};

ORIs     = { ...
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

RBSs     = { ...
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
% 2. Parámetros globales
% =============================================================

proms = unique(Promoters, 'stable');
rbss  = unique(RBSs,      'stable');
oris  = unique(ORIs,      'stable');

paramNames = {};
% omega_promoter
for i = 1:numel(proms)
    paramNames{end+1} = ['\omega_{' proms{i} '}']; %#ok<AGROW>
end
idxEnd_omega = numel(paramNames);

% thetaK_rbs, 
for i = 1:numel(rbss)
    paramNames{end+1} = ['\theta^K_{' rbss{i} '}']; %#ok<AGROW>
end
idxEnd_thetaK = numel(paramNames);
% thetaS_rbs
for i = 1:numel(rbss)
    paramNames{end+1} = ['\theta^S_{' rbss{i} '}']; %#ok<AGROW>
end
idxEnd_thetaS = numel(paramNames);

% N_ori
for i = 1:numel(oris)
    paramNames{end+1} = ['N_{' oris{i} '}']; %#ok<AGROW>
end

nParams = numel(paramNames);


%% =============================================================
% 3. Matriz global M_full(TUs x params)
% =============================================================

M_full = zeros(nTUs, nParams);
for i = 1:nTUs
    prom = Promoters{i};
    rbs  = RBSs{i};
    ori  = ORIs{i};
    M_full(i, strcmp(paramNames,  ['\omega_{' prom '}'])) = 1;
    M_full(i, strcmp(paramNames, ['\theta^K_{' rbs '}'])) = 1;
    M_full(i, strcmp(paramNames, ['\theta^S_{' rbs '}'])) = 1;
    M_full(i, strcmp(paramNames, ['N_{' ori '}']))      = 1;
end


%% =============================================================
% 4. Sublibrería J23100
% =============================================================

idxJ = find(strcmp(Promoters, 'J23100'));
IDs_J = IDs(idxJ);
ORIs_J = ORIs(idxJ);
Promoters_J = Promoters(idxJ);
RBSs_J = RBSs(idxJ);

proms_J = unique(Promoters_J, 'stable');
rbss_J  = unique(RBSs_J,      'stable');
oris_J  = unique(ORIs_J,      'stable');

paramNames_J = {};
for i=1:numel(proms_J)
    paramNames_J{end+1} = ['\omega_{' proms_J{i} '}']; %#ok<AGROW>
end
idxEnd_omega_J = numel(paramNames_J);
for i=1:numel(rbss_J)
    paramNames_J{end+1} = ['\theta^K_{' rbss_J{i} '}']; %#ok<AGROW>
end
idxEnd_thetaK_J = numel(paramNames_J);

for i=1:numel(rbss_J)
    paramNames_J{end+1} = ['\theta^S_{' rbss_J{i} '}']; %#ok<AGROW>
end
idxEnd_thetaS_J = numel(paramNames_J);

for i=1:numel(oris_J)
    paramNames_J{end+1} = ['N_{' oris_J{i} '}']; %#ok<AGROW>
end

nParams_J = numel(paramNames_J);
nTUs_J    = numel(idxJ);

M_J = zeros(nTUs_J, nParams_J);
for ii = 1:nTUs_J
    prom = Promoters_J{ii};
    rbs  = RBSs_J{ii};
    ori  = ORIs_J{ii};
    M_J(ii, strcmp(paramNames_J,  ['\omega_{' prom '}'])) = 1;
    M_J(ii, strcmp(paramNames_J, ['\theta^K_{' rbs '}'])) = 1;
    M_J(ii, strcmp(paramNames_J, ['\theta^S_{' rbs '}'])) = 1;
    M_J(ii, strcmp(paramNames_J, ['N_{' ori '}']))      = 1;
end

%% =============================================================
% 5. Plot comparativo con gris y líneas de bloques
% =============================================================

color_grey = [184 185 182]/255;   % #b8b9b6
cmap = [1 1 1; color_grey];

figure('Units','pixels','Position',[50 50 1500 600]);

% ---------- PANEL A: librería completa ----------
subplot(1,2,1);
imagesc(M_full);
colormap(cmap);
clim([0 1]);
hold on;

% Líneas verticales entre bloques de parámetros
colColor = [0.6 0.6 0.6]; % color de línea de bloque
lw = 1.25;
xBlocks = [idxEnd_omega, idxEnd_thetaK, idxEnd_thetaS];
for xb = xBlocks
    line([xb+0.5 xb+0.5], [0.5 nTUs+0.5], 'Color', colColor, 'LineWidth', lw);
end

% Líneas horizontales entre librerías: L24(1-24), L30(25-30), L5(31-35)
rowColor = [0.3 0.3 0.3];
lwRow = 2;
rowBounds = [24, 30];  % después de L24 y después de L30
for rb = rowBounds
    line([0.5 nParams+0.5], [rb+0.5 rb+0.5], 'Color', rowColor, 'LineWidth', lwRow);
end

ax = gca;
ax.XTick = 1:nParams;
ax.XTickLabel = paramNames;
ax.XTickLabelRotation = 90;
ax.YTick = 1:nTUs;
ax.YTickLabel = IDs;
ax.FontSize = 18;
set(gca,'TickLength',[0 0]);

xlabel('Parameters');
ylabel('TUs');
title('A. Full library (L24 + L30 + L5)');

% ---------- PANEL B: sublibrería J23100 ----------
subplot(1,2,2);
imagesc(M_J);
colormap(cmap);
clim([0 1]);
hold on;

% Líneas verticales entre bloques omega | thetaK | thetaS | N
xBlocks_J = [idxEnd_omega_J, idxEnd_thetaK_J, idxEnd_thetaS_J];
for xb = xBlocks_J
    line([xb+0.5 xb+0.5], [0.5 nTUs_J+0.5], 'Color', colColor, 'LineWidth', lw);
end

ax = gca;
ax.XTick = 1:nParams_J;
ax.XTickLabel = paramNames_J;
ax.XTickLabelRotation = 90;
ax.YTick = 1:nTUs_J;
ax.YTickLabel = IDs_J;
ax.FontSize = 20;
set(gca,'TickLength',[0 0]);

xlabel('Parameters');
ylabel('TUs');
title('B. J23100 sublibrary only');

% Guardar
print('-dpng','-r300','Jacobian_comparison_full_vs_J23100_blocks.png');

end
