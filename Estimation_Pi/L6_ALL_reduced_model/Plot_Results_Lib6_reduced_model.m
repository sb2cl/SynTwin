function Plot_Results_Lib6_reduced_model(Results_Lib,Lower_Bounds,Upper_Bounds,indices_plasmids,indices_promoters,indices_rbss,title_text,varargin)
%% Plot_Results_Lib6_reduced_model
% Generate publication-oriented plots for the Lib6 reduced model.
%
% The function produces:
%   - the B0034 kappa^0 estimate distribution;
%   - a 2-by-3 tiled comparison of experimental and predicted synthesis rates.
%
% The original seven positional arguments are retained for compatibility.
% Additional output, font, histogram, and figure-size settings are supplied
% through name-value arguments.
%
% The expected Lib6 tensor dimensions are 2 plasmids by 3 promoters by one
% RBS. The canonical Lib30 RBS index for B0034 is 3.

% Keep old signature variables for compatibility, even if not used below.
if nargin < 7 || isempty(title_text), title_text = ''; end %#ok<NASGU>
if nargin < 3, Upper_Bounds = []; end %#ok<NASGU>
if nargin < 2, Lower_Bounds = []; end %#ok<NASGU>

p = inputParser;
p.addParameter('OutputDir', fullfile(pwd,'Figures'), @(x)ischar(x)||isstring(x));
p.addParameter('SaveFigures', true, @(x)islogical(x)||isnumeric(x));
p.addParameter('FigurePrefix', 'Lib6_ALL_reduced', @(x)ischar(x)||isstring(x));
p.addParameter('BaseFontSize', 26, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('TileFontSize', 18, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('AxisLabelFontSize', 24, @(x)isnumeric(x)&&isscalar(x));
p.addParameter('HistogramBinMethod','auto', @(x)ischar(x)||isstring(x));
p.addParameter('ExportPause',2, @(x)isnumeric(x)&&isscalar(x));
% Figure sizes in centimeters: [width height]. These control aspect ratio.
p.addParameter('FigCopySizeCm',     [28 22], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('FigPromoterSizeCm', [72 24], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('FigRBSSizeCm',      [84 24], @(x)isnumeric(x)&&numel(x)==2);
p.addParameter('FigTilesSizeCm',    [76 48], @(x)isnumeric(x)&&numel(x)==2);
p.parse(varargin{:});
opt = p.Results;
opt.OutputDir = char(opt.OutputDir);
opt.FigurePrefix = char(opt.FigurePrefix);
opt.HistogramBinMethod = char(opt.HistogramBinMethod);
opt.FigCopySizeCm     = double(reshape(opt.FigCopySizeCm,1,2));
opt.FigPromoterSizeCm = double(reshape(opt.FigPromoterSizeCm,1,2));
opt.FigRBSSizeCm      = double(reshape(opt.FigRBSSizeCm,1,2));
opt.FigTilesSizeCm    = double(reshape(opt.FigTilesSizeCm,1,2));

if opt.SaveFigures && ~exist(opt.OutputDir,'dir')
    mkdir(opt.OutputDir);
end

% -------------------------------------------------------------------------
% Canonical colors. These are the authoritative colors used across figures.
% -------------------------------------------------------------------------
color_blue        = hex2rgb_local('#00008B');
color_grey        = hex2rgb_local('#4b4b4b');
color_grey_boira  = hex2rgb_local('#F2F2F2'); %#ok<NASGU>
color_grey_neutre = hex2rgb_local('#CCCCCC');
RBS_colors_hex      = ['#A8DADC';'#F4A261';'#B5E48C';'#CDB4DB';'#FFE066'];
Promoter_colors_hex = ['#264653';'#E76F51';'#2A9D8F';'#6A0572'];
RBS_colors      = hex2rgb_local(RBS_colors_hex);
Promoter_colors = hex2rgb_local(Promoter_colors_hex);

num_plasmids  = length(indices_plasmids);
num_promoters = length(indices_promoters);
num_rbss      = length(indices_rbss);


%% ESTIMATED PARAMETERS HISTOGRAMS

% -------------------------------------------------------------------------
% 3) RBS kappa^0 histograms
% -------------------------------------------------------------------------
fig1 = figure('Name','Lib6 RBS kappa histograms','Color','w', ...
    'Units','centimeters','Position',[1 1 opt.FigRBSSizeCm]);
tl1 = tiledlayout(fig1,1,num_rbss,'Padding','compact','TileSpacing','compact'); %#ok<NASGU>
for k = 1:num_rbss
    ax = nexttile;
    hold(ax,'on');
    [ii,jj] = find_result_with_field2(Results_Lib,'RBS_k0_sigma0_raw',1,1);
    Data_raw  = Results_Lib{ii,jj}.RBS_k0_sigma0_raw(:);
    Data_mean = Results_Lib{ii,jj}.RBS_k0_sigma0_mean;
    ridx = indices_rbss(k);
    C = RBS_colors(ridx,:);
    plot_histogram_counts(ax, Data_raw, C, opt.HistogramBinMethod);
    xline(ax,Data_mean,'Color','k','LineWidth',2.4);
    if isfield(Results_Lib{ii,jj},'TU_RBS')
        title(ax, Results_Lib{ii,jj}.TU_RBS, 'FontSize', opt.BaseFontSize, 'FontWeight','bold', 'Interpreter','none');
    end
    xlabel(ax,'$\kappa^0 = K^0/\sigma^0\; (\mathrm{molec}^{-1})$','FontSize',opt.BaseFontSize,'Interpreter','latex');
    if k == 1
        ylabel(ax,'Count','FontSize',opt.BaseFontSize,'Interpreter','latex');
    end
    set_axis_publication(ax,1.2*opt.BaseFontSize);
end
save_figure_if_requested(fig1,opt.OutputDir,[opt.FigurePrefix '_hist_RBS_kappa'],opt.SaveFigures,opt.ExportPause);

% -------------------------------------------------------------------------
% 2) Tiled experimental vs predicted Pi(mu)
% -------------------------------------------------------------------------
fig2 = figure('Name','Lib6 experimental vs predicted Pi(mu) tiles','Color','w', ...
    'Units','centimeters','Position',[1 1 opt.FigTilesSizeCm]);
tl2 = tiledlayout(fig2,num_plasmids,num_promoters, ...
    'Padding','compact','TileSpacing','compact'); %#ok<NASGU>

for i = 1:num_plasmids
    for j = 1:num_promoters
        for k = 1:num_rbss
            if isempty(Results_Lib{i,j,k})
                continue;
            end
            tile_num = j + num_promoters*(i-1) ;
            ax = nexttile(tile_num);
            hold(ax,'on');

            E = Results_Lib{i,j,k};
            C = color_from_result_or_rbs(E, RBS_colors, indices_rbss(k));
            ec_label = l6_code(i,j,k);

            [xmu, ypred_l, ypred_u, ypred_m] = extract_prediction_slices(E);
            [x_g,y_g,y_std] = extract_experimental_global(E);
            if isempty(xmu) || isempty(x_g)
                continue;
            end

            % Experimental mean +/- 2 std envelope first.
            % IMPORTANT: for the x-window we deliberately reproduce the
            % original S6.4 code convention:
            %     x_gn = x_g(end:-1:1);
            %     ax.XLim = [x_gn(1), x_gn(end)];
            % Do not sort x_g and do not use min/max over all values to define
            % the displayed range. This makes the visible mu interval match
            % the experimental trajectory interval used in the old figure.
            xg = x_g(:); yg = y_g(:); ys = y_std(:);
            x_gn = xg(end:-1:1);
            y_gn = yg(end:-1:1);
            ys_gn = ys(end:-1:1);
            yln = max(y_gn - 2*ys_gn, 0);
            yun = y_gn + 2*ys_gn;

            xlim_pair = [x_gn(1), x_gn(end)];
            if ~all(isfinite(xlim_pair)) || xlim_pair(1) == xlim_pair(2)
                xlim_pair = [xg(1), xg(end)];
            end
            if xlim_pair(1) > xlim_pair(2)
                xlim_pair = fliplr(xlim_pair);
            end

            fill(ax,[x_gn; flipud(x_gn)], [yln; flipud(yun)], color_grey, ...
                'FaceAlpha',0.20,'EdgeColor','none','HandleVisibility','off');
            plot(ax,xg,yg,'-','Color',color_grey_neutre,'LineWidth',1.6,'HandleVisibility','off');

            % Predicted uncertainty envelope. It is intentionally not
            % re-windowed by min/max; the axes below clip it to the exact
            % experimental interval, as in the original plotting code.
            fill(ax,[xmu; flipud(xmu)], [ypred_l; flipud(ypred_u)], C, ...
                'FaceAlpha',0.50,'EdgeColor','none','HandleVisibility','off');
            plot(ax,xmu,ypred_m,'-','Color',C,'LineWidth',2.4,'HandleVisibility','off');
            plot(ax,xmu,ypred_l,'-','Color',C,'LineWidth',0.9,'HandleVisibility','off');
            plot(ax,xmu,ypred_u,'-','Color',C,'LineWidth',0.9,'HandleVisibility','off');

            % EC label above the tile, matching the style used in the
            % experimental Lib35 supplementary figure.
            title(ax, ec_label, 'FontSize', opt.TileFontSize+2, ...
                'FontWeight','bold', 'Interpreter','none');

            grid(ax,'on');
            box(ax,'on');
            set_axis_publication(ax,1.2*opt.TileFontSize);
            ax.LineWidth = 0.9;

            % Apply the original experimental x-limits. Ticks are exactly
            % three values inside that visible interval.
            if ~all(isfinite(xlim_pair)) || xlim_pair(1) == xlim_pair(2)
                xlim_pair = [min(xmu,[],'omitnan') max(xmu,[],'omitnan')];
            end
            xlim(ax,xlim_pair);
            [xticks_use, xlabels_use] = three_axis_ticks_linear(xlim_pair(1),xlim_pair(2),'mu');
            ax.XTick = xticks_use;
            ax.XTickLabel = xlabels_use;

            % Y-limits must be computed after fixing the visible mu window.
            % Otherwise prediction/experimental values outside xlim can inflate
            % the y-axis and push the visible curves to the bottom of the tile.
            xlo = xlim_pair(1);
            xhi = xlim_pair(2);
            exp_visible = isfinite(x_gn) & x_gn >= xlo & x_gn <= xhi;
            pred_visible = isfinite(xmu) & xmu >= xlo & xmu <= xhi;

            y_candidates = [];
            if any(exp_visible)
                y_candidates = [y_candidates; yun(exp_visible); y_gn(exp_visible)]; %#ok<AGROW>
            end
            if any(pred_visible)
                y_candidates = [y_candidates; ypred_u(pred_visible); ypred_m(pred_visible)]; %#ok<AGROW>
            end
            y_candidates = y_candidates(isfinite(y_candidates));
            ymax_visible = max(y_candidates,[],'omitnan');
            if ~isfinite(ymax_visible) || ymax_visible <= 0
                ymax_visible = max([ypred_u(:); yun(:); y_gn(:)],[],'omitnan');
            end
            if ~isfinite(ymax_visible) || ymax_visible <= 0
                ymax_visible = 1;
            end

            % Add a modest headroom before computing the three readable ticks.
            [yticks_use, ylabels_use] = three_axis_ticks_linear(0,1.10*ymax_visible,'pi');
            ylim(ax,[yticks_use(1) yticks_use(end)]);
            ax.YTick = yticks_use;
            ax.YTickLabel = ylabels_use;

            ax.XTickLabelRotation = 35;
            ax.TickLabelInterpreter = 'none';

            if j == 1
                ylabel(ax,'$\Pi_p$','FontSize',opt.AxisLabelFontSize,'Interpreter','latex');
            end
            if i == num_plasmids
                xlabel(ax,'$\mu$','FontSize',opt.AxisLabelFontSize,'Interpreter','latex');
            end
        end
    end
end
save_figure_if_requested(fig2,opt.OutputDir,[opt.FigurePrefix '_Pi_vs_mu_tiles'],opt.SaveFigures,opt.ExportPause);

end % function Plot_Results_Lib6_reduced_model


% ========================================================================
function plot_histogram_counts(ax, data, color_rgb, bin_method)
    data = data(:);
    data = data(isfinite(data));
    histogram(ax,data,'BinMethod',bin_method, ...
        'FaceColor',color_rgb,'EdgeColor','none','FaceAlpha',0.95);
    grid(ax,'on');
    box(ax,'on');
end

% ========================================================================
function set_axis_publication(ax,font_size)
    ax.FontSize = font_size;
    ax.LineWidth = 1.1;
    ax.TickDir = 'out';
    ax.TickLabelInterpreter = 'latex';
    ax.XLabel.FontSize = font_size;
    ax.YLabel.FontSize = font_size;
end

% ========================================================================
function [ii,jj] = find_result_with_field2(R,field,ipref,jpref)
    sz = size(R);
    ii = min(ipref,sz(1)); jj = min(jpref,sz(2)); 
    if ~isempty(R{ii,jj}) && isfield(R{ii,jj},field)
        return;
    end
    for i = 1:sz(1)
        for j = 1:sz(2)
                if ~isempty(R{i,j}) && isfield(R{i,j},field)
                    ii = i; jj = j; 
                    return;
                end
        end
    end
    error('Could not find field "%s" in Results_Lib.',field);
end

% ========================================================================
function C = color_from_result_or_rbs(E,RBS_colors,rbs_index)
    C = RBS_colors(rbs_index,:);
    if isfield(E,'TU_color_code') && ~isempty(E.TU_color_code)
        try
            C = ensure_color_rgb(E.TU_color_code);
        catch
            % keep canonical RBS color fallback
        end
    end
end

% ========================================================================
function [xmu, yl, yu, ym] = extract_prediction_slices(E)
    xmu = []; yl = []; yu = []; ym = [];
    if ~isfield(E,'MC_mu_slices') || isempty(E.MC_mu_slices)
        return;
    end
    n = length(E.MC_mu_slices);
    xmu = nan(n,1); yl = nan(n,1); yu = nan(n,1); ym = nan(n,1);
    for q = 1:n
        xmu(q) = E.MC_mu_slices(q).Mu_slice;
        yl(q)  = E.MC_mu_slices(q).Pi_pred_q2p5;
        yu(q)  = E.MC_mu_slices(q).Pi_pred_q97p5;
        ym(q)  = E.MC_mu_slices(q).Pi_pred_q50;
    end
    valid = isfinite(xmu) & isfinite(yl) & isfinite(yu) & isfinite(ym);
    xmu = xmu(valid); yl = yl(valid); yu = yu(valid); ym = ym(valid);
    [xmu,ord] = sort(xmu);
    yl = yl(ord); yu = yu(ord); ym = ym(ord);
end

% ========================================================================
function [xg,yg,ys] = extract_experimental_global(E)
    xg = []; yg = []; ys = [];
    if isfield(E,'Mu_mumax_pmax_global_mean') && isfield(E,'Pi_mumax_pmax_global_mean')
        xg = E.Mu_mumax_pmax_global_mean(:);
        yg = E.Pi_mumax_pmax_global_mean(:);
        if isfield(E,'Pi_mumax_pmax_global_std')
            ys = E.Pi_mumax_pmax_global_std(:);
        else
            ys = zeros(size(yg));
        end
        n = min([numel(xg),numel(yg),numel(ys)]);
        xg = xg(1:n); yg = yg(1:n); ys = ys(1:n);
        valid = isfinite(xg) & isfinite(yg) & isfinite(ys);
        xg = xg(valid); yg = yg(valid); ys = ys(valid);
        % Preserve original ordering. The tile x-limits intentionally follow
        % the original code convention using x_g(end:-1:1), rather than
        % min/max over a sorted vector.
    end
end

% ========================================================================
function code = l6_code(i,j,k)
%L6_CODE Return EC code for library L6.
%
% INTERNAL tensor convention:
%   i = 1 -> pGreen
%   i = 2 -> pSC101
%
% Official L6/Table S1 convention:
%   L6_01--L6_3 -> pSC101
%   L6_4--L6_6 -> pGreen

    if i == 1
        iL = 2;   % pGreen
    elseif i == 2
        iL = 1;   % pSC101
    else
        error('Invalid plasmid index i.');
    end

    idx = (iL-1)*3 + j;
    code = sprintf('L6_%02d', idx);
end


% ========================================================================
function [ticks, labels] = three_mu_ticks(vmin,vmax)
    % Three readable ticks for growth-rate axes. Boundaries are rounded to
    % nearby 0.001 min^-1 values so the displayed range remains essentially
    % the experimental range while avoiding overcrowding.
    step = 0.001;
    lo = floor(vmin/step)*step;
    hi = ceil(vmax/step)*step;
    if hi <= lo
        lo = vmin; hi = vmax;
    end
    mid = round(((lo+hi)/2)/step)*step;
    if mid <= lo || mid >= hi
        mid = (lo+hi)/2;
    end
    ticks = unique([lo mid hi],'stable');
    labels = arrayfun(@(x)sprintf('%.3f',x),ticks,'UniformOutput',false);
end

% ========================================================================
function [ticks, labels] = three_pi_ticks(vmin,vmax)
    % Three readable ticks for Pi axes, avoiding scientific notation at zero.
    lo = max(0,vmin);
    hi = nice_axis_upper(vmax);
    if hi <= lo
        hi = vmax;
    end
    mid = nice_axis_middle(lo,hi);
    ticks = unique([lo mid hi],'stable');
    labels = arrayfun(@format_numeric_tick,ticks,'UniformOutput',false);
end

% ========================================================================
function hi = nice_axis_upper(vmax)
    if ~isfinite(vmax) || vmax <= 0
        hi = 1; return;
    end
    exponent = floor(log10(vmax));
    base = 10^exponent;
    candidates = [1 2 2.5 5 10] * base;
    idx = find(candidates >= vmax,1,'first');
    if isempty(idx)
        hi = 10*base;
    else
        hi = candidates(idx);
    end
end

% ========================================================================
function mid = nice_axis_middle(lo,hi)
    raw = (lo+hi)/2;
    step = nice_tick_step(raw);
    mid = round(raw/step)*step;
    if mid <= lo || mid >= hi
        mid = raw;
    end
end

% ========================================================================
function force_three_ticks(ax,which_axis,vmin,vmax,mode)
    if ~isfinite(vmin) || ~isfinite(vmax) || vmin == vmax
        return;
    end
    lo = min(vmin,vmax);
    hi = max(vmin,vmax);

    if strcmpi(mode,'mu')
        % Prefer readable 0.005 multiples for growth rate.
        step = 0.005;
        t0 = floor(lo/step)*step;
        t2 = ceil(hi/step)*step;
        if t0 < 0, t0 = 0; end
        if t2 <= t0
            t0 = lo; t2 = hi;
        end
        t1 = round(((t0+t2)/2)/step)*step;
        if t1 <= t0 || t1 >= t2
            t1 = (t0+t2)/2;
        end
        ticks = [t0 t1 t2];
        labels = arrayfun(@(x)sprintf('%.3f',x),ticks,'UniformOutput',false);
    else
        % Pi values: readable numbers spanning local data range.
        step = nice_tick_step((hi-lo)/2);
        t0 = floor(lo/step)*step;
        t2 = ceil(hi/step)*step;
        if t0 < 0, t0 = 0; end
        if t2 <= t0
            t0 = lo; t2 = hi;
        end
        t1 = round(((t0+t2)/2)/step)*step;
        if t1 <= t0 || t1 >= t2
            t1 = (t0+t2)/2;
        end
        ticks = [t0 t1 t2];
        labels = arrayfun(@format_numeric_tick,ticks,'UniformOutput',false);
    end

    if strcmpi(which_axis,'x')
        ax.XTick = ticks;
        ax.XTickLabel = labels;
    else
        ax.YTick = ticks;
        ax.YTickLabel = labels;
    end
end

% ========================================================================
function step = nice_tick_step(raw_step)
    if ~isfinite(raw_step) || raw_step <= 0
        step = 1;
        return;
    end
    exponent = floor(log10(raw_step));
    fraction = raw_step / 10^exponent;
    if fraction <= 1
        nice_fraction = 1;
    elseif fraction <= 2
        nice_fraction = 2;
    elseif fraction <= 2.5
        nice_fraction = 2.5;
    elseif fraction <= 5
        nice_fraction = 5;
    else
        nice_fraction = 10;
    end
    step = nice_fraction * 10^exponent;
end

% ========================================================================
function s = format_numeric_tick(v)
    if abs(v) < 1e-12
        s = '0';
        return;
    end
    av = abs(v);
    if abs(v-round(v)) < 1e-9 && av >= 1
        s = sprintf('%.0f',v);
    elseif av >= 100
        s = sprintf('%.0f',v);
    elseif av >= 10
        s = sprintf('%.0f',v);
    elseif av >= 1
        s = sprintf('%.1f',v);
    elseif av >= 0.01
        s = sprintf('%.3f',v);
    else
        s = sprintf('%.1e',v);
    end
end

% ========================================================================
function [ticks, labels] = three_axis_ticks_linear(vmin,vmax,mode)
    % Exactly three ticks spanning the actual plotted range: min, middle, max.
    % For mu, values are displayed with 3 decimals. For Pi, labels avoid
    % scientific notation and unnecessary decimals.
    if ~isfinite(vmin) || ~isfinite(vmax) || vmin == vmax
        ticks = [vmin vmax];
        labels = arrayfun(@format_numeric_tick,ticks,'UniformOutput',false);
        return;
    end
    lo = min(vmin,vmax);
    hi = max(vmin,vmax);
    if strcmpi(mode,'pi')
        lo = max(0,lo);
        hi = nice_axis_upper(hi);
    end
    mid = (lo + hi)/2;
    ticks = [lo mid hi];
    if strcmpi(mode,'mu')
        labels = arrayfun(@(x)sprintf('%.3f',x),ticks,'UniformOutput',false);
    else
        labels = arrayfun(@format_numeric_tick,ticks,'UniformOutput',false);
    end
end

% ========================================================================
function save_figure_if_requested(fig,outdir,basename,do_save,export_pause)
    if ~do_save
        return;
    end
    if nargin < 5 || isempty(export_pause), export_pause = 0.25; end
    if ~exist(outdir,'dir')
        mkdir(outdir);
    end
    drawnow;
    pause(export_pause);
    pngfile = fullfile(outdir,[basename '.png']);
    svgfile = fullfile(outdir,[basename '.svg']);
    try
        exportgraphics(fig,pngfile,'Resolution',300,'BackgroundColor','white');
    catch
        print(fig,pngfile,'-dpng','-r300');
    end
    drawnow;
    pause(export_pause);
    try
        exportgraphics(fig,svgfile,'ContentType','vector','BackgroundColor','white');
    catch
        print(fig,svgfile,'-dsvg');
    end
end

% ========================================================================
function rgb = ensure_color_rgb(c)
    if ischar(c) || isstring(c)
        rgb = hex2rgb_local(c);
        return;
    end
    if isnumeric(c)
        if numel(c) == 3
            rgb = double(reshape(c,1,3));
            if max(rgb) > 1
                rgb = rgb/255;
            end
            return;
        end
    end
    error('Unsupported color format.');
end

% ========================================================================
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
