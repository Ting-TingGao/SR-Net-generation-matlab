function stats = analyze_sr_network(Line_list, meta, varargin)
%Analysis for SR network
%
%   stats = analyze_sr_network(Line_list, meta, ...)
%   Computes and (optionally) plots the following:
%     • Adjacency matrix A (touch-based)
%     • Length per line L_t  + log-binned length distribution (bins_length)
%     • Thickness series (meta.v_hist) + running average + optional theory,
%       and log-binned thickness distribution (bins_thickness)
%     • Degree distribution (log-binned, fitted by curve_fit_log)
%
% OPTIONS (name-value):
%   'bins_degree'     30     # of logarithmic bins for degree distribution
%   'bins_length'     20     # of logarithmic bins for length distribution
%   'bins_thickness'  20     # of logarithmic bins for thickness distribution (v_hist)
%   'start_deg_fit'    1     First bin index used in degree fit
%   'start_len_fit'    1     First bin index used in length fit
%   'start_thk_fit'    1     First bin index used in thickness fit
%   'tolerance'     1e-8     Touching tolerance (point-to-line distance)
%   'show_figures'  true     Make figures
%   'save_plots'    false    If true, save PDFs with prefix
%   'save_prefix'     ''     Filename prefix when saving plots
%   'close_figures' false    Close figures at the end
%
% OUTPUT (stats struct):
%   .A                 (t x t) adjacency matrix
%   .degree            (t x 1) degree vector
%   .L_t               (t x 1) line lengths
%   .thickness         struct with fields:
%                        - v_hist: per-step thickness series (if available)
%                        - ave_thickness: running average of v_hist
%                        - theo: optional theoretical <lambda>(t) if (alpha,v0) valid
%   .degree_dist       struct: mid/pdf and fit (popt_log, pcov_log)
%   .length_dist       struct: mid/pdf and fit (popt_log, pcov_log)
%   .thickness_dist    struct: mid/pdf and fit (popt_log, pcov_log)
%   .alpha, .v0        copied from meta.opts if present (for record)
%
% Author: Ting-Ting Gao 

p = inputParser;
addParameter(p,'bins_degree',30,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'bins_length',20,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'bins_thickness',20,@(x)isnumeric(x)&&isscalar(x)&&x>=2);
addParameter(p,'start_deg_fit',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'start_len_fit',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'start_thk_fit',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'tolerance',1e-8,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'show_figures',true,@(x)islogical(x)&&isscalar(x));
addParameter(p,'save_plots',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'save_prefix','',@(x)ischar(x)||isstring(x));
addParameter(p,'close_figures',false,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:});
opt = p.Results;

% ---- Unpack meta (optional) ----
alpha_ = getfield_default(meta,'opts',struct());
alpha_ = getfield_default(alpha_,'alpha',NaN);
v0_    = getfield_default(meta,'opts',struct());
v0_    = getfield_default(v0_,'v0',NaN);

% ---- Basic sizes ----
if isempty(Line_list)
    t = 0;
else
    t = size(Line_list,1);
end

% ---- 1) Length per line ----
L_t = zeros(t,1);
for i = 1:t
    L_t(i) = hypot(Line_list(i,1)-Line_list(i,3), Line_list(i,2)-Line_list(i,4));
end

% ---- 2) Adjacency & degree (touching criterion) ----
A = zeros(t,t);
for i = 1:t
    for j = i+1:t
        % i's left/right edges vs j's midline
        isL = checkLineTouching(Line_list(i,5:6),  Line_list(i,7:8),  Line_list(j,1:2), Line_list(j,3:4), opt.tolerance);
        isR = checkLineTouching(Line_list(i,9:10), Line_list(i,11:12), Line_list(j,1:2), Line_list(j,3:4), opt.tolerance);
        % j's left/right edges vs i's midline (symmetric)
        isL2 = checkLineTouching(Line_list(j,5:6),  Line_list(j,7:8),  Line_list(i,1:2), Line_list(i,3:4), opt.tolerance);
        isR2 = checkLineTouching(Line_list(j,9:10), Line_list(j,11:12), Line_list(i,1:2), Line_list(i,3:4), opt.tolerance);
        if isL || isR || isL2 || isR2
            A(i,j) = 1; A(j,i) = 1;
        end
    end
end

degree_vec = sum(A,2);

% ---- 3) Degree distribution (log-binned) ----
nonzero_deg = degree_vec(degree_vec>0);
[deg_mid, deg_pdf] = get_distribution(nonzero_deg, opt.bins_degree);
fit_idx_d = max(1,opt.start_deg_fit):numel(deg_mid);
if ~isempty(fit_idx_d)
    xfit = deg_mid(fit_idx_d); yfit = deg_pdf(fit_idx_d);
    pos = (xfit>0) & (yfit>0) & isfinite(xfit) & isfinite(yfit);
    if any(pos)
        [poptD, pcovD] = curve_fit_log(xfit(pos), yfit(pos), 1);
    else
        poptD = [NaN NaN]; pcovD = [NaN NaN; NaN NaN];
    end
else
    poptD = [NaN NaN]; pcovD = [NaN NaN; NaN NaN];
end

% ---- 4) Length distribution (log-binned) ----
[len_mid, len_pdf] = get_distribution(L_t(L_t>0), opt.bins_length);
fit_idx_l = max(1,opt.start_len_fit):numel(len_mid);
if ~isempty(fit_idx_l)
    xfit = len_mid(fit_idx_l); yfit = len_pdf(fit_idx_l);
    pos = (xfit>0) & (yfit>0) & isfinite(xfit) & isfinite(yfit);
    if any(pos)
        [poptL, pcovL] = curve_fit_log(xfit(pos), yfit(pos), 1);
    else
        poptL = [NaN NaN]; pcovL = [NaN NaN; NaN NaN];
    end
else
    poptL = [NaN NaN]; pcovL = [NaN NaN; NaN NaN];
end

% ---- 5) Thickness distribution ----
th = struct('v_hist',[],'ave_thickness',[],'theo',[]);
if isfield(meta,'v_hist') && ~isempty(meta.v_hist)
    th.v_hist = meta.v_hist(:);
    th.ave_thickness = cumsum(th.v_hist) ./ (1:numel(th.v_hist))';
    if isfinite(alpha_) && isfinite(v0_) && abs(1-alpha_)>1e-9
        n = numel(th.v_hist);
        Tt = (1:n)';
        th.theo = (v0_/(1-alpha_) * Tt.^(1-alpha_) - v0_/(1-alpha_)) ./ Tt;
    end
end

[thk_mid, thk_pdf] = get_distribution(th.v_hist(th.v_hist>0), opt.bins_thickness);
fit_idx_t = max(1,opt.start_thk_fit):numel(thk_mid);
if ~isempty(fit_idx_t)
    xfit = thk_mid(fit_idx_t); yfit = thk_pdf(fit_idx_t);
    pos = (xfit>0) & (yfit>0) & isfinite(xfit) & isfinite(yfit);
    if any(pos)
        [poptT, pcovT] = curve_fit_log(xfit(pos), yfit(pos), 1);
    else
        poptT = [NaN NaN]; pcovT = [NaN NaN; NaN NaN];
    end
else
    poptT = [NaN NaN]; pcovT = [NaN NaN; NaN NaN];
end

% ---- 6) PLOTS ----
figs = struct();
if opt.show_figures
    % Degree distribution
    if ~isempty(deg_mid)
        figs.deg = figure('Name','Degree distribution');
        loglog(deg_mid, deg_pdf, 'o', 'MarkerSize', 6, 'MarkerFaceColor',[0.60,0.45,0.64], 'MarkerEdgeColor',[0.60,0.45,0.64]); hold on
        yfit = 10.^(poptD(1) + poptD(2)*log10(deg_mid(fit_idx_d)));
        loglog(deg_mid(fit_idx_d), yfit, '--', 'LineWidth', 1.5, 'Color',[0.41,0.24,0.45], 'DisplayName', sprintf('slope = %.3f', poptD(2)));
        legend('show'); xlabel('k'); ylabel('p(k)'); title('Degree distribution'); hold off
        maybe_print(figs.deg, opt.save_plots, opt.save_prefix, 'degree_dist');
    end

    % Length distribution
    if ~isempty(len_mid)
        figs.len = figure('Name','Length distribution');
        loglog(len_mid, len_pdf, 'o', 'MarkerSize', 6, 'MarkerFaceColor',[0.69,0.82,0.97], 'MarkerEdgeColor',[0.48,0.66,0.86]); hold on
        yfit = 10.^(poptL(1) + poptL(2)*log10(len_mid(fit_idx_l)));
        loglog(len_mid(fit_idx_l), yfit, '--', 'LineWidth', 1.5, 'Color',[0.18,0.64,0.76], 'DisplayName', sprintf('slope = %.3f', poptL(2)));
        legend('show'); xlabel('L'); ylabel('p(L)'); title('Line length distribution'); hold off
        maybe_print(figs.len, opt.save_plots, opt.save_prefix, 'length_dist');
    end

    % Thickness distribution
    if ~isempty(thk_mid)
        figs.thk = figure('Name','Thickness distribution');
        loglog(thk_mid, thk_pdf, 'o', 'MarkerSize', 6, 'MarkerFaceColor',[0.69,0.82,0.97], 'MarkerEdgeColor',[0.48,0.66,0.86]); hold on
        yfit = 10.^(poptT(1) + poptT(2)*log10(thk_mid(fit_idx_t)));
        loglog(thk_mid(fit_idx_t), yfit, '--', 'LineWidth', 1.5, 'Color',[0.86,0.71,0.47], 'DisplayName', sprintf('slope = %.3f', poptT(2)));
        legend('show'); xlabel('\lambda'); ylabel('p(\lambda)'); title('Thickness distribution'); hold off
        maybe_print(figs.thk, opt.save_plots, opt.save_prefix, 'thickness_dist');
    end
end

if opt.close_figures && ~isempty(fieldnames(figs))
    fns = fieldnames(figs);
    for k = 1:numel(fns), try close(figs.(fns{k})); end, end
end

% ---- Pack outputs ----
stats = struct();
stats.A           = A;
stats.degree      = degree_vec;
stats.L_t         = L_t;
stats.thickness   = th;
stats.degree_dist = struct('mid',deg_mid,'pdf',deg_pdf,'fit',struct('popt_log',poptD,'pcov_log',pcovD));
stats.length_dist = struct('mid',len_mid,'pdf',len_pdf,'fit',struct('popt_log',poptL,'pcov_log',pcovL));
stats.thickness_dist = struct('mid',thk_mid,'pdf',thk_pdf,'fit',struct('popt_log',poptT,'pcov_log',pcovT));
stats.alpha       = alpha_;
stats.v0          = v0_;

end

% ============================ Helpers ============================
function maybe_print(figH, do_save, prefix, tag)
% Save figure to PDF if requested; safe even if printing fails
if do_save && ~isempty(prefix)
    try
        fn = sprintf('%s_%s.pdf', char(prefix), tag);
        print(figH, fn, '-dpdf');
    catch ME
        warning('maybe_print:saveFailed','Failed to save figure %s: %s', tag, ME.message);
    end
end
end

function [x, y] = get_distribution(data_sequence, number_of_bins)
    % Log-binned PDF as provided
    if nargin < 2, number_of_bins = 20; end
    data_sequence = data_sequence(:)';
    if isempty(data_sequence)
        x = []; y = []; return; 
    end
    lower_bound = min(data_sequence);
    upper_bound = max(data_sequence);
    if lower_bound <= 0 || ~isfinite(lower_bound) || ~isfinite(upper_bound) || lower_bound==upper_bound
        x = []; y = []; return;
    end
    lower_bound = log10(lower_bound);
    upper_bound = log10(upper_bound) + 1;
    bins = logspace(lower_bound, upper_bound, number_of_bins);
    [y, ~] = histcounts(data_sequence, bins, 'Normalization', 'pdf');
    x = bins(2:end) - diff(bins) / 2.0;
    drop_indices = (y == 0) | ~isfinite(y) | ~isfinite(x);
    x(drop_indices) = [];
    y(drop_indices) = [];
end

function [popt_log, pcov_log] = curve_fit_log(xdata, ydata, intercept)
    % Fit log10(y) = a + b*log10(x)  (uses Curve Fitting Toolbox)
    xdata_log = log10(xdata);
    ydata_log = log10(ydata);
    if intercept == 1
        ft = fittype('a + b * x');
        opts = fitoptions('Method', 'NonlinearLeastSquares');
        opts.StartPoint = [1, 1];
        fitresult = fit(xdata_log(:), ydata_log(:), ft, opts);
        popt_log = coeffvalues(fitresult);
        pcov_log = confint(fitresult);
    else
        ft = fittype('b * x');
        opts = fitoptions('Method', 'NonlinearLeastSquares');
        opts.StartPoint = 1;
        fitresult = fit(xdata_log(:), ydata_log(:), ft, opts);
        popt_log = coeffvalues(fitresult);
        pcov_log = confint(fitresult);
    end
end

function v = getfield_default(S, f, default)
if isfield(S,f), v = S.(f); else, v = default; end
end

function isTouching = checkLineTouching(p1, p2, p3, p4, tolerance)
    if nargin < 5, tolerance = 1e-8; end
    if pointToLineDistance(p1, p3, p4) < tolerance, isTouching = true; return; end
    if pointToLineDistance(p2, p3, p4) < tolerance, isTouching = true; return; end
    if pointToLineDistance(p3, p1, p2) < tolerance, isTouching = true; return; end
    if pointToLineDistance(p4, p1, p2) < tolerance, isTouching = true; return; end
    isTouching = false;
end

function dist = pointToLineDistance(p, v1, v2)
    lineVec = v2 - v1; pointVec = p - v1; L = norm(lineVec);
    if L == 0, dist = norm(p - v1); return; end
    u = dot(pointVec, lineVec) / (L^2); u = max(0,min(1,u));
    proj = v1 + u*lineVec; dist = norm(p - proj);
end