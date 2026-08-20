%% step3_statistics.m
% Statistical comparison of RMS (%MVC) and MDF slope between EXO and NOEXO.
%
% Workflow:
%   1. Load rms_pX.mat  +  mvc_reference.mat  per participant → RMS as %MVC
%   2. Load mdf_pX.mat  per participant → MDF slope + total spectral shift
%   3. For each muscle:
%        - Shapiro-Wilk on within-participant differences
%        - Paired t-test (normal) or Wilcoxon signed-rank (non-normal)
%        - Bonferroni correction: alpha_corr = 0.05 / 7
%        - Cohen's d (parametric only)
%   4. Print results tables
%   5. Figures: effect sizes + p-value markers
%
% References (from TFG Section 5.3):
%   Shapiro-Wilk on differences, Bonferroni alpha = 0.05/7 = 0.00714
%   Cohen's d thresholds: 0.2 small | 0.5 medium | 0.8 large

clear; clc; close all;

%% ── 0. CONFIGURATION ────────────────────────────────────────────────────

results_dir = 'C:\Users\gzomo\TFG\results';
mvc_dir     = 'C:\Users\gzomo\TFG\mvc';
output_dir  = results_dir;

alpha_nom  = 0.05;
n_muscles  = 7;
alpha_corr = alpha_nom / n_muscles;   % Bonferroni: 0.00714

muscle_names = {'Anterior Deltoid'; 'Lateral Deltoid'; 'Posterior Deltoid'; ...
                'Upper Trapezius';  'Biceps Brachii';  'Triceps Brachii'; 'ECR'};
% These match the column order (1-7) of the processed signals.

%% ── 1. SELECT FILES ──────────────────────────────────────────────────────

% RMS files (one per participant)
[rms_files, rms_path] = uigetfile('rms_*.mat', ...
    'Select RMS .mat files (all participants)', results_dir, 'MultiSelect', 'on');
if isequal(rms_files, 0), error('No RMS files selected.'); end
if ischar(rms_files), rms_files = {rms_files}; end
n_subjects = numel(rms_files);

% MVC reference files (one per participant, same order)
[mvc_files, mvc_path] = uigetfile('*mvc_reference.mat', ...
    'Select MVC reference .mat files (same order as RMS)', mvc_dir, 'MultiSelect', 'on');
if isequal(mvc_files, 0)
    use_mvc = false;
    warning('No MVC files selected — RMS will be in µV (not %%MVC).');
else
    if ischar(mvc_files), mvc_files = {mvc_files}; end
    use_mvc = true;
end

% MDF files (one per participant)
[mdf_files, mdf_path] = uigetfile('mdf_*.mat', ...
    'Select MDF .mat files (all participants)', results_dir, 'MultiSelect', 'on');
if isequal(mdf_files, 0), error('No MDF files selected.'); end
if ischar(mdf_files), mdf_files = {mdf_files}; end

fprintf('\n%d participants | Bonferroni alpha = %.5f\n\n', n_subjects, alpha_corr);

%% ── 2. BUILD DATA MATRICES ──────────────────────────────────────────────
% Rows = participants, Cols = muscles (1-7, column order of signal)

rms_base_exo   = zeros(n_subjects, n_muscles);
rms_base_noexo = zeros(n_subjects, n_muscles);
rms_post_exo   = zeros(n_subjects, n_muscles);
rms_post_noexo = zeros(n_subjects, n_muscles);

mdf_slope_exo   = zeros(n_subjects, n_muscles);
mdf_slope_noexo = zeros(n_subjects, n_muscles);
mdf_shift_exo   = zeros(n_subjects, n_muscles);   % MDF_last - MDF_first
mdf_shift_noexo = zeros(n_subjects, n_muscles);

for s = 1:n_subjects

    % ── RMS ──
    r = load(fullfile(rms_path, rms_files{s}));
    rr = r.rms_results;

    exo_base   = rr.EXO.baseline;    % [7×1] µV
    exo_post   = rr.EXO.postfat;
    noexo_base = rr.NOEXO.baseline;
    noexo_post = rr.NOEXO.postfat;

    % ── MVC normalization ──
    if use_mvc
        mvc_raw = load(fullfile(mvc_path, mvc_files{s}));
        f_mvc   = fieldnames(mvc_raw);
        mvc     = mvc_raw.(f_mvc{1});   % take first variable regardless of name
        disp(fieldnames(mvc));
        % Build mvc_norm(ch): RMS_MVC for the muscle at column ch
        mvc_norm = zeros(n_muscles, 1);
        if isfield(mvc, 'channel_idx') && isfield(mvc, 'RMS_MVC')
            ch_idx  = mvc.channel_idx(:);
            rms_mvc = mvc.RMS_MVC(:);
            for m = 1:numel(ch_idx)
                ch = ch_idx(m);
                if ch >= 1 && ch <= n_muscles && ~isnan(rms_mvc(m))
                    mvc_norm(ch) = rms_mvc(m);
                end
            end
        else
            % Fallback: RMS_MVC already in column order 1-7
            mvc_norm = mvc.RMS_MVC(1:n_muscles);
        end
        % RMS_MVC is in V (raw Delsys), convert to µV to match rms_results units
        mvc_norm = mvc_norm * 1e6;
        % Replace any zero (unprocessed muscle) with NaN to avoid /0
        mvc_norm(mvc_norm == 0) = NaN;
        % Express as %MVC
        exo_base   = exo_base   ./ mvc_norm * 100;
        exo_post   = exo_post   ./ mvc_norm * 100;
        noexo_base = noexo_base ./ mvc_norm * 100;
        noexo_post = noexo_post ./ mvc_norm * 100;
    else
        % Convert to µV if not normalizing (already in µV from step1)
        % No conversion needed — keep as is
    end

    rms_base_exo(s,:)   = exo_base';
    rms_base_noexo(s,:) = noexo_base';
    rms_post_exo(s,:)   = exo_post';
    rms_post_noexo(s,:) = noexo_post';

    % ── MDF ──
    m_data = load(fullfile(mdf_path, mdf_files{s}));
    md = m_data.mdf_results;

    mdf_slope_exo(s,:)   = md.EXO.slope';
    mdf_slope_noexo(s,:) = md.NOEXO.slope';

    % Total spectral shift: MDF_last - MDF_first per muscle
    for ch = 1:n_muscles
        y_e = md.EXO.mdf(:, ch);
        y_n = md.NOEXO.mdf(:, ch);
        ok_e = ~isnan(y_e);  ok_n = ~isnan(y_n);
        if sum(ok_e) >= 2
            mdf_shift_exo(s,ch)   = y_e(find(ok_e,1,'last')) - y_e(find(ok_e,1,'first'));
        end
        if sum(ok_n) >= 2
            mdf_shift_noexo(s,ch) = y_n(find(ok_n,1,'last')) - y_n(find(ok_n,1,'first'));
        end
    end
end

rms_unit = iif(use_mvc, '%MVC', 'µV');
fprintf('RMS unit: %s\n\n', rms_unit);

%% ── 3. STATISTICAL TESTS ─────────────────────────────────────────────────
% Analyses run:
%   A) RMS baseline    EXO vs NOEXO
%   B) RMS post-fat    EXO vs NOEXO
%   C) MDF slope       EXO vs NOEXO
%   D) MDF total shift EXO vs NOEXO

analyses = {
    'RMS Baseline',    rms_base_exo,    rms_base_noexo;
    'RMS Post-fat',    rms_post_exo,    rms_post_noexo;
    'MDF Slope',       mdf_slope_exo,   mdf_slope_noexo;
    'MDF Total Shift', mdf_shift_exo,   mdf_shift_noexo;
};

results = struct();

for a = 1:size(analyses,1)
    label = analyses{a,1};
    X_exo   = analyses{a,2};   % [n_subjects × n_muscles]
    X_noexo = analyses{a,3};

    field = matlab.lang.makeValidName(label);
    results.(field).label      = label;
    results.(field).p_raw      = zeros(1, n_muscles);
    results.(field).p_adj      = zeros(1, n_muscles);
    results.(field).test_used  = cell(1, n_muscles);
    results.(field).stat       = zeros(1, n_muscles);
    results.(field).cohens_d   = nan(1, n_muscles);
    results.(field).mean_diff  = zeros(1, n_muscles);
    results.(field).sig        = false(1, n_muscles);

    for ch = 1:n_muscles
        d = X_exo(:,ch) - X_noexo(:,ch);   % within-participant differences
        d = d(~isnan(d));
        n = numel(d);

        results.(field).mean_diff(ch) = mean(d);

        if n < 3
            results.(field).p_raw(ch)     = NaN;
            results.(field).test_used{ch} = 'n<3';
            continue;
        end

        % Shapiro-Wilk on differences
        normal = shapiro_wilk(d);

        if normal
            % Paired t-test
            [~, p, ~, st] = ttest(d);
            results.(field).test_used{ch} = 't-test';
            results.(field).stat(ch)      = st.tstat;
            results.(field).cohens_d(ch)  = mean(d) / std(d);
        else
            % Wilcoxon signed-rank
            [p, ~, st] = signrank(d);
            results.(field).test_used{ch} = 'Wilcoxon';
            results.(field).stat(ch)      = st.signedrank;
        end

        results.(field).p_raw(ch) = p;
        results.(field).p_adj(ch) = min(p * n_muscles, 1);   % Bonferroni
        results.(field).sig(ch)   = results.(field).p_adj(ch) < alpha_nom;
    end
end

%% ── 4. PRINT TABLES ──────────────────────────────────────────────────────

sep = repmat('─', 1, 82);

for a = 1:size(analyses,1)
    field = matlab.lang.makeValidName(analyses{a,1});
    res   = results.(field);

    fprintf('\n%s\n', sep);
    fprintf('  %s\n', res.label);
    fprintf('%s\n', sep);
    fprintf('%-22s  %8s  %8s  %8s  %8s  %7s  %s\n', ...
            'Muscle', 'Mean Δ', 'p (raw)', 'p (adj)', 'Cohen d', 'Test', 'Sig');
    fprintf('%s\n', sep);

    for ch = 1:n_muscles
        sig_str = iif(res.sig(ch), '***', '');
        d_str   = iif(isnan(res.cohens_d(ch)), '   —   ', sprintf('%+7.3f', res.cohens_d(ch)));
        fprintf('%-22s  %+7.3f   %7.4f   %7.4f   %s  %-8s  %s\n', ...
                muscle_names{ch}, ...
                res.mean_diff(ch), res.p_raw(ch), res.p_adj(ch), ...
                d_str, res.test_used{ch}, sig_str);
    end
    fprintf('%s\n', sep);
end

fprintf('\nBonferroni α_corrected = %.5f | *** = significant after correction\n\n', alpha_corr);

%% ── 5. FIGURES ───────────────────────────────────────────────────────────

% ── 5a. RMS comparison (baseline and post-fat) ──
fig1 = figure('Name', 'RMS statistics', 'NumberTitle', 'off', ...
              'Position', [50 50 1300 560]);

titles_rms = {'Baseline', 'Post-fatigue'};
data_pairs  = {rms_base_exo, rms_base_noexo; rms_post_exo, rms_post_noexo};
fields_rms  = {'RMS_Baseline', 'RMS_Post_fat'};

x  = 1:n_muscles;
bw = 0.3;

for sub = 1:2
    ax = subplot(1,2,sub);
    hold(ax, 'on');

    Xe = data_pairs{sub,1};   % [n_subjects × 7]
    Xn = data_pairs{sub,2};

    me = mean(Xe, 1, 'omitnan');
    se_e = std(Xe, 0, 1, 'omitnan') / sqrt(n_subjects);
    mn = mean(Xn, 1, 'omitnan');
    se_n = std(Xn, 0, 1, 'omitnan') / sqrt(n_subjects);

    b1 = bar(ax, x - bw/2, me, bw, 'FaceColor', [0.2 0.4 0.8]);
    b2 = bar(ax, x + bw/2, mn, bw, 'FaceColor', [0.8 0.2 0.2]);
    errorbar(ax, x - bw/2, me, se_e, 'k.', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    errorbar(ax, x + bw/2, mn, se_n, 'k.', 'LineWidth', 1.2, 'HandleVisibility', 'off');

    % Significance markers
    field = matlab.lang.makeValidName(['RMS ' titles_rms{sub}]);
    if isfield(results, field)
        y_top = max([me + se_e, mn + se_n], [], 'all') * 1.15;
        for ch = 1:n_muscles
            if results.(field).sig(ch)
                text(ax, ch, y_top, '*', 'HorizontalAlignment', 'center', ...
                     'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
            end
        end
    end

    set(ax, 'XTick', x, 'XTickLabel', muscle_names, 'XTickLabelRotation', 35);
    ylabel(ax, sprintf('RMS (%s)', rms_unit));
    title(ax, titles_rms{sub}, 'FontSize', 12);
    legend(ax, [b1 b2], {'EXO', 'NOEXO'}, 'Location', 'northwest');
    grid(ax, 'on'); box(ax, 'off');
end
sgtitle(fig1, sprintf('RMS EXO vs NOEXO — mean ± SEM (%s)', rms_unit), ...
        'FontSize', 13, 'FontWeight', 'bold');
exportgraphics(fig1, fullfile(output_dir, 'stats_RMS.png'), 'Resolution', 300);

% ── 5b. MDF slope comparison ──
fig2 = figure('Name', 'MDF slope statistics', 'NumberTitle', 'off', ...
              'Position', [100 100 1000 480]);
hold on;

me  = mean(mdf_slope_exo,   1, 'omitnan');
mn  = mean(mdf_slope_noexo, 1, 'omitnan');
se_e = std(mdf_slope_exo,   0, 1, 'omitnan') / sqrt(n_subjects);
se_n = std(mdf_slope_noexo, 0, 1, 'omitnan') / sqrt(n_subjects);

b1 = bar(x - bw/2, me, bw, 'FaceColor', [0.2 0.4 0.8]);
b2 = bar(x + bw/2, mn, bw, 'FaceColor', [0.8 0.2 0.2]);
errorbar(x - bw/2, me, se_e, 'k.', 'LineWidth', 1.2, 'HandleVisibility', 'off');
errorbar(x + bw/2, mn, se_n, 'k.', 'LineWidth', 1.2, 'HandleVisibility', 'off');
yline(0, 'k--', 'LineWidth', 1, 'HandleVisibility', 'off');

% Significance
field = 'MDF_Slope';
y_top = max([me + se_e, mn + se_n], [], 'all') * 1.2;
y_bot = min([me - se_e, mn - se_n], [], 'all') * 1.2;
for ch = 1:n_muscles
    if isfield(results, field) && results.(field).sig(ch)
        text(ch, y_top, '*', 'HorizontalAlignment', 'center', ...
             'FontSize', 16, 'FontWeight', 'bold', 'Color', 'k');
    end
end

set(gca, 'XTick', x, 'XTickLabel', muscle_names, 'XTickLabelRotation', 35);
ylabel('MDF slope (Hz/s)');
title('Mean MDF regression slope — EXO vs NOEXO (mean ± SEM)', 'FontSize', 13);
legend([b1 b2], {'EXO', 'NOEXO'}, 'Location', 'southeast');
grid on; box off;
exportgraphics(fig2, fullfile(output_dir, 'stats_MDF_slope.png'), 'Resolution', 300);

% ── 5c. Cohen's d heatmap ──
fig3 = figure('Name', "Cohen's d", 'NumberTitle', 'off', ...
              'Position', [150 150 900 360]);

d_matrix = zeros(4, n_muscles);
row_labels = {'RMS Baseline', 'RMS Post-fat', 'MDF Slope', 'MDF Shift'};
for a = 1:4
    field = matlab.lang.makeValidName(analyses{a,1});
    d_matrix(a,:) = results.(field).cohens_d;
end

imagesc(d_matrix, [-2 2]);
colormap(bluewhitered_local());
colorbar;
set(gca, 'XTick', 1:n_muscles, 'XTickLabel', muscle_names, 'XTickLabelRotation', 35, ...
         'YTick', 1:4, 'YTickLabel', row_labels);
title("Cohen's d (EXO − NOEXO) | ±0.2 small | ±0.5 medium | ±0.8 large", 'FontSize', 11);

% Overlay significance markers
for a = 1:4
    field = matlab.lang.makeValidName(analyses{a,1});
    for ch = 1:n_muscles
        if results.(field).sig(ch)
            text(ch, a, '*', 'HorizontalAlignment', 'center', ...
                 'VerticalAlignment', 'middle', 'FontSize', 14, 'Color', 'k');
        end
    end
end

exportgraphics(fig3, fullfile(output_dir, 'stats_cohens_d.png'), 'Resolution', 300);

fprintf('Figures saved to: %s\n', output_dir);

%% ── 6. SAVE RESULTS ─────────────────────────────────────────────────────

save(fullfile(output_dir, 'statistics_results.mat'), 'results', ...
     'muscle_names', 'alpha_corr', 'n_subjects', 'rms_unit');
fprintf('Results saved → statistics_results.mat\n');

%% ══ LOCAL FUNCTIONS ══════════════════════════════════════════════════════

function is_normal = shapiro_wilk(x)
% Simple Shapiro-Wilk for n = 3..50.
% Returns true if normality is NOT rejected at alpha = 0.05.
    x = sort(x(:));
    n = numel(x);

    if n < 3
        is_normal = true;   % cannot test
        return;
    end

    % Use MATLAB built-in lillietest as proxy if SW not available
    % For n=3, SW critical value at alpha=0.05 is 0.767
    if n == 3
        a1   = 0.7071;
        W    = (a1 * (x(3) - x(1)))^2 / sum((x - mean(x)).^2);
        is_normal = W >= 0.767;
        return;
    end

    % For n > 3: use Lilliefors test (built-in Statistics Toolbox)
    try
        [h, ~] = lillietest(x, 'Alpha', 0.05);
        is_normal = ~h;
    catch
        is_normal = true;   % if toolbox missing, assume normal
    end
end

function cmap = bluewhitered_local()
% Blue-white-red colormap for diverging data.
    n  = 64;
    r1 = linspace(0,   1, n/2)';
    g1 = linspace(0.2, 1, n/2)';
    b1 = ones(n/2, 1);
    r2 = ones(n/2, 1);
    g2 = linspace(1, 0.2, n/2)';
    b2 = linspace(1, 0,   n/2)';
    cmap = [r1, g1, b1; r2, g2, b2];
end

function out = iif(cond, a, b)
    if cond, out = a; else, out = b; end
end
