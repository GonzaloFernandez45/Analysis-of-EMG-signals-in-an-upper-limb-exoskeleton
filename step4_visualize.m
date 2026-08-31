% =========================================================================
% step4_visualize.m
% =========================================================================
% Publication-ready figures for the TFG.
%
% FIGURES PRODUCED:
%   1. Baseline RMS %MVC  - EXO vs NOEXO (bars + individual paired dots)
%   2. Post-fatigue RMS %MVC - EXO vs NOEXO (bars + individual paired dots)
%   3. Delta MDF (Hz)     - EXO vs NOEXO (bars + individual paired dots)
%   4. MDF trajectory across the time-normalised fatigue block (0-100%)
%   5. Endurance          - n_reps + block_duration_s (bars + paired dots)
%
% INPUTS:
%   step3_results.mat  - output of step3_statistics.m
%   gente/ root folder - for Figure 4 (loads mdf5_*.mat per subject)
%
% Significance markers:
%   ***  = p < alpha_bonf  (Bonferroni-corrected)
%   dagger = p < 0.05 (nominal, not corrected)
%
% Compatibility: MATLAB R2016b+  (uses print() - no exportgraphics)
% =========================================================================

clear; clc; close all;
rng(1);   % reproducible jitter positions

% -------------------------------------------------------------------------
% Style
% -------------------------------------------------------------------------
COL_EXO   = [0.18 0.38 0.75];   % blue
COL_NOEXO = [0.78 0.18 0.18];   % red
COL_DOT   = [0.20 0.20 0.20];   % dark gray for individual dots
COL_LINE  = [0.70 0.70 0.70];   % light gray for connecting lines

DASH      = char(8211);   % en-dash (encoding-safe)
PM        = char(177);    % plus-minus sign (encoding-safe)
DAGGER    = char(8224);   % dagger symbol for nominal significance

set(groot, ...
    'defaultAxesFontName',  'Helvetica', ...
    'defaultAxesFontSize',   10, ...
    'defaultAxesLineWidth',  0.8, ...
    'defaultAxesTickDir',    'out', ...
    'defaultAxesBox',        'off', ...
    'defaultFigureColor',    'w');

muscle_names = {'AD','LD','PD','UT','BB','TB','ECR'};
muscle_long  = {'Ant. Deltoid','Lat. Deltoid','Post. Deltoid', ...
                'Upper Trap.','Biceps Br.','Triceps Br.','ECR'};
N_CH = 7;

% -------------------------------------------------------------------------
% 1. Load step3_results.mat
% -------------------------------------------------------------------------
results_dir = 'C:\Users\gzomo\TFG\results';
[res_fname, res_fpath] = uigetfile('step3_results.mat', ...
    'Select step3_results.mat', results_dir);
if isequal(res_fname, 0), error('No results file selected.'); end

R = load(fullfile(res_fpath, res_fname));
fprintf('Loaded: %s\n', fullfile(res_fpath, res_fname));

% Paired data matrices
rms_base_exo   = R.rms_base_exo;
rms_base_noexo = R.rms_base_noexo;
rms_post_exo   = R.rms_post_exo;
rms_post_noexo = R.rms_post_noexo;
dmdf_exo       = R.dmdf_exo;
dmdf_noexo     = R.dmdf_noexo;
nreps_exo      = R.nreps_exo;
nreps_noexo    = R.nreps_noexo;
dur_exo        = R.dur_exo;
dur_noexo      = R.dur_noexo;

N = size(rms_base_exo, 1);   % number of subjects

% Output folder: same as results file
output_dir = res_fpath;
fprintf('Figures will be saved to: %s\n', output_dir);

% =========================================================================
% FIGURE 1 - Baseline RMS %MVC
% =========================================================================
fig1 = plot_paired_bars(rms_base_exo, rms_base_noexo, ...
    R.res_base, muscle_long, ...
    ['Baseline RMS (%MVC) ' DASH ' EXO vs NOEXO'], 'RMS (%MVC)', ...
    R.ALPHA_BONF, DAGGER, COL_EXO, COL_NOEXO, COL_DOT, COL_LINE);

save_fig(fig1, fullfile(output_dir, 'fig1_baseline_rms.png'));
fprintf('Fig 1 saved.\n');

% =========================================================================
% FIGURE 2 - Post-fatigue RMS %MVC
% =========================================================================
fig2 = plot_paired_bars(rms_post_exo, rms_post_noexo, ...
    R.res_post, muscle_long, ...
    ['Post-fatigue RMS (%MVC) ' DASH ' EXO vs NOEXO'], 'RMS (%MVC)', ...
    R.ALPHA_BONF, DAGGER, COL_EXO, COL_NOEXO, COL_DOT, COL_LINE);

save_fig(fig2, fullfile(output_dir, 'fig2_postfatigue_rms.png'));
fprintf('Fig 2 saved.\n');

% =========================================================================
% FIGURE 3 - Delta MDF
% =========================================================================
fig3 = plot_paired_bars(dmdf_exo, dmdf_noexo, ...
    R.res_mdf, muscle_long, ...
    ['\DeltaMDF (Hz) ' DASH ' EXO vs NOEXO'], '\DeltaMDF (Hz)', ...
    R.ALPHA_BONF, DAGGER, COL_EXO, COL_NOEXO, COL_DOT, COL_LINE);

% Add horizontal reference line at 0
for ax_idx = 1:N_CH
    ax_tmp = subplot(1, N_CH, ax_idx);
    hold(ax_tmp, 'on');
    yline(0, '--', 'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off');
end

save_fig(fig3, fullfile(output_dir, 'fig3_delta_mdf.png'));
fprintf('Fig 3 saved.\n');

% =========================================================================
% FIGURE 4 - MDF trajectory normalized 0-100%
% =========================================================================
fprintf('\nLoading mdf5 files for trajectory figure...\n');

root_gente = uigetdir(fileparts(strtrim(res_fpath)), ...
    'Select root subjects folder (gente/) for MDF trajectory');
if ~isequal(root_gente, 0)
    fig4 = plot_mdf_trajectory(root_gente, muscle_names, COL_EXO, COL_NOEXO, PM);
    if ~isempty(fig4)
        save_fig(fig4, fullfile(output_dir, 'fig4_mdf_trajectory.png'));
        fprintf('Fig 4 saved.\n');
    end
else
    fprintf('Skipping Fig 4 (no folder selected).\n');
end

% =========================================================================
% FIGURE 5 - Endurance outcomes
% =========================================================================
fig5 = figure('Name', 'Fig5 Endurance', 'NumberTitle', 'off', ...
              'Position', [100 100 700 480]);

ok_r = ~isnan(nreps_exo) & ~isnan(nreps_noexo);
ok_d = ~isnan(dur_exo)   & ~isnan(dur_noexo);

outcomes   = {nreps_exo(ok_r), nreps_noexo(ok_r); ...
              dur_exo(ok_d),   dur_noexo(ok_d)};
titles_end = {'Number of repetitions', 'Block duration (s)'};
sig_bonf   = [false, false];   % endurance tested at alpha=0.05 (no Bonferroni)
sig_nom    = [R.res_end.nreps.sig,  R.res_end.duration.sig];

bw = 0.3;
for sp = 1:2
    ax = subplot(1, 2, sp);
    hold(ax, 'on');

    ve  = outcomes{sp, 1};
    vn  = outcomes{sp, 2};
    me  = mean(ve);  mn  = mean(vn);
    sem_e = std(ve) / sqrt(numel(ve));
    sem_n = std(vn) / sqrt(numel(vn));

    b1 = bar(ax, 1-bw/2, me, bw, 'FaceColor', COL_EXO,   'FaceAlpha', 0.75, 'EdgeColor', 'none');
    b2 = bar(ax, 1+bw/2, mn, bw, 'FaceColor', COL_NOEXO, 'FaceAlpha', 0.75, 'EdgeColor', 'none');
    errorbar(ax, 1-bw/2, me, sem_e, 'k.', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    errorbar(ax, 1+bw/2, mn, sem_n, 'k.', 'LineWidth', 1.2, 'HandleVisibility', 'off');

    % Individual paired dots
    jit = 0.04;
    for s_idx = 1:numel(ve)
        jx = (rand - 0.5) * jit;
        plot(ax, [1-bw/2+jx, 1+bw/2+jx], [ve(s_idx), vn(s_idx)], ...
             '-', 'Color', COL_LINE, 'LineWidth', 0.8, 'HandleVisibility', 'off');
        plot(ax, 1-bw/2+jx, ve(s_idx), 'o', 'MarkerFaceColor', COL_DOT, ...
             'MarkerEdgeColor', 'w', 'MarkerSize', 5, 'HandleVisibility', 'off');
        plot(ax, 1+bw/2+jx, vn(s_idx), 'o', 'MarkerFaceColor', COL_DOT, ...
             'MarkerEdgeColor', 'w', 'MarkerSize', 5, 'HandleVisibility', 'off');
    end

    % Significance marker
    all_vals = [ve; vn];
    y_range  = max(all_vals) - min(all_vals);
    y_top    = max(all_vals) + 0.12 * y_range;
    if sig_bonf(sp)
        sig_str = '***';
    elseif sig_nom(sp)
        sig_str = DAGGER;
    else
        sig_str = '';
    end
    if ~isempty(sig_str)
        y_line = max(all_vals) + 0.07 * y_range;
        plot(ax, [1-bw/2, 1+bw/2], [y_line, y_line], 'k-', ...
             'LineWidth', 0.8, 'HandleVisibility', 'off');
        text(ax, 1, y_top, sig_str, 'HorizontalAlignment', 'center', ...
             'FontSize', 13, 'FontWeight', 'bold', 'Color', 'k');
    end

    set(ax, 'XTick', [], 'XLim', [0.5 1.5]);
    ylabel(ax, titles_end{sp});
    title(ax, titles_end{sp}, 'FontSize', 10, 'FontWeight', 'bold');
    grid(ax, 'on'); ax.GridAlpha = 0.15;

    if sp == 1
        legend(ax, [b1 b2], {'EXO', 'NOEXO'}, 'Location', 'northeast', 'Box', 'off');
    end
end

sgtitle(fig5, ['Endurance outcomes ' DASH ' EXO vs NOEXO (mean ' PM ' SEM + individual)'], ...
        'FontSize', 11, 'FontWeight', 'bold');

save_fig(fig5, fullfile(output_dir, 'fig5_endurance.png'));
fprintf('Fig 5 saved.\n');

fprintf('\nAll figures saved to: %s\n', output_dir);

% =========================================================================
% LOCAL FUNCTIONS  (must be at end of script - MATLAB R2016b+)
% =========================================================================

function fig = plot_paired_bars(mat_exo, mat_noexo, res, mlong, ...
                                fig_title, y_label, alpha_bonf, dagger, ...
                                col_e, col_n, col_dot, col_line)
% Grouped bar chart per muscle with individual paired dots.
% Significance markers from res.p_test and alpha_bonf.
% *** = Bonferroni significant; dagger = nominal p < 0.05

    N_CH = size(mat_exo, 2);
    fig  = figure('Name', fig_title, 'NumberTitle', 'off', ...
                  'Position', [50 50 1500 440]);

    bw  = 0.32;
    jit = 0.06;

    for ch = 1:N_CH
        ax = subplot(1, N_CH, ch);
        hold(ax, 'on');

        e_ok  = mat_exo(:, ch);
        no_ok = mat_noexo(:, ch);
        valid = ~isnan(e_ok) & ~isnan(no_ok);
        ev    = e_ok(valid);
        nv    = no_ok(valid);
        n_ok  = sum(valid);

        if n_ok < 3
            title(ax, mlong{ch}, 'FontSize', 8);
            continue;
        end

        me    = mean(ev);    mn    = mean(nv);
        sem_e = std(ev) / sqrt(n_ok);
        sem_n = std(nv) / sqrt(n_ok);

        % Bars
        b1 = bar(ax, 1-bw/2, me, bw, 'FaceColor', col_e, 'FaceAlpha', 0.75, 'EdgeColor', 'none');
        b2 = bar(ax, 1+bw/2, mn, bw, 'FaceColor', col_n, 'FaceAlpha', 0.75, 'EdgeColor', 'none');

        % SEM error bars
        errorbar(ax, 1-bw/2, me, sem_e, 'k.', 'LineWidth', 1.2, 'HandleVisibility', 'off');
        errorbar(ax, 1+bw/2, mn, sem_n, 'k.', 'LineWidth', 1.2, 'HandleVisibility', 'off');

        % Individual paired dots + connecting lines
        for s = 1:n_ok
            jx = (rand - 0.5) * jit;
            plot(ax, [1-bw/2+jx, 1+bw/2+jx], [ev(s), nv(s)], ...
                 '-', 'Color', col_line, 'LineWidth', 0.8, 'HandleVisibility', 'off');
            plot(ax, 1-bw/2+jx, ev(s), 'o', 'MarkerFaceColor', col_dot, ...
                 'MarkerEdgeColor', 'w', 'MarkerSize', 4, 'HandleVisibility', 'off');
            plot(ax, 1+bw/2+jx, nv(s), 'o', 'MarkerFaceColor', col_dot, ...
                 'MarkerEdgeColor', 'w', 'MarkerSize', 4, 'HandleVisibility', 'off');
        end

        % Significance marker
        if isfield(res, 'p_test') && ~isnan(res.p_test(ch))
            p        = res.p_test(ch);
            all_vals = [ev; nv];
            y_range  = max(all_vals) - min(all_vals);
            y_top    = max(all_vals) + 0.12 * y_range;
            if p < alpha_bonf
                sig_str = '***';
            elseif p < 0.05
                sig_str = dagger;
            else
                sig_str = '';
            end
            if ~isempty(sig_str)
                y_line = max(all_vals) + 0.07 * y_range;
                plot(ax, [1-bw/2, 1+bw/2], [y_line, y_line], 'k-', ...
                     'LineWidth', 0.8, 'HandleVisibility', 'off');
                text(ax, 1, y_top, sig_str, 'HorizontalAlignment', 'center', ...
                     'FontSize', 12, 'FontWeight', 'bold', 'Color', 'k');
            end
        end

        set(ax, 'XTick', [], 'XLim', [0.5 1.5]);
        title(ax, mlong{ch}, 'FontSize', 8, 'FontWeight', 'bold');
        grid(ax, 'on'); ax.GridAlpha = 0.15;

        if ch == 1
            ylabel(ax, y_label, 'FontSize', 9);
            legend(ax, [b1 b2], {'EXO', 'NOEXO'}, 'Location', 'best', 'Box', 'off', 'FontSize', 7);
        end
    end

    sgtitle(fig, fig_title, 'FontSize', 11, 'FontWeight', 'bold');
end


function fig = plot_mdf_trajectory(root, mnames, col_e, col_n, pm)
% Load mdf5_*.mat from all subject subfolders and plot mean +/- SEM MDF
% trajectory across the time-normalised fatigue block (0-100%).
% SEM is computed point-wise using only valid (non-NaN) observations.
% No extrapolation: NaN is kept where data is absent.

    N_GRID = 101;   % normalised time grid: 0%, 1%, ..., 100%
    N_CH   = numel(mnames);
    grid_t = linspace(0, 100, N_GRID);

    traj_exo   = {};
    traj_noexo = {};

    entries   = dir(root);
    is_subdir = [entries.isdir] & ~ismember({entries.name}, {'.', '..'});
    subj_dirs = entries(is_subdir);

    for di = 1:numel(subj_dirs)
        subj_folder = fullfile(root, subj_dirs(di).name);
        mdf5_d = dir(fullfile(subj_folder, 'results', 'results', 'mdf5_*.mat'));

        for fi = 1:numel(mdf5_d)
            S = load(fullfile(mdf5_d(fi).folder, mdf5_d(fi).name));

            if isfield(S, 'condition')
                cond = upper(strtrim(S.condition));
            elseif ~isempty(regexpi(mdf5_d(fi).name, '_NOEXO'))
                cond = 'NOEXO';
            else
                cond = 'EXO';
            end

            mdf_r = S.mdf_reps;   % [n_reps x 7]
            n_r   = size(mdf_r, 1);
            if n_r < 3, continue; end

            t_orig     = linspace(0, 100, n_r);
            interp_mat = nan(N_GRID, N_CH);
            for ch = 1:N_CH
                y  = mdf_r(:, ch);
                ok = ~isnan(y);
                if sum(ok) < 3, continue; end
                % No extrapolation: grid points outside observed range stay NaN
                interp_mat(:, ch) = interp1(t_orig(ok), y(ok), grid_t, 'linear');
            end

            if strcmp(cond, 'EXO')
                traj_exo{end+1} = interp_mat; %#ok<AGROW>
            else
                traj_noexo{end+1} = interp_mat; %#ok<AGROW>
            end
        end
    end

    if isempty(traj_exo) || isempty(traj_noexo)
        fprintf('  Not enough data for trajectory figure.\n');
        fig = [];
        return;
    end

    n_e = numel(traj_exo);
    n_n = numel(traj_noexo);
    arr_e = nan(n_e, N_GRID, N_CH);
    arr_n = nan(n_n, N_GRID, N_CH);
    for i = 1:n_e, arr_e(i,:,:) = traj_exo{i}; end
    for i = 1:n_n, arr_n(i,:,:) = traj_noexo{i}; end

    % Point-wise mean and SEM (valid N per grid point per channel)
    mu_e  = squeeze(nanmean(arr_e, 1));   % [N_GRID x N_CH]
    mu_n  = squeeze(nanmean(arr_n, 1));
    nv_e  = squeeze(sum(~isnan(arr_e), 1));
    nv_n  = squeeze(sum(~isnan(arr_n), 1));
    sem_e = squeeze(nanstd(arr_e, 0, 1)) ./ sqrt(max(nv_e, 1));
    sem_n = squeeze(nanstd(arr_n, 0, 1)) ./ sqrt(max(nv_n, 1));

    fig = figure('Name', 'Fig4 MDF Trajectory', 'NumberTitle', 'off', ...
                 'Position', [60 60 1500 600]);

    for ch = 1:N_CH
        ax = subplot(2, 4, ch);
        hold(ax, 'on');

        % Shaded SEM band - skip NaN points
        x_fill = [grid_t, fliplr(grid_t)];
        y_up_e  = mu_e(:,ch)' + sem_e(:,ch)';
        y_lo_e  = mu_e(:,ch)' - sem_e(:,ch)';
        y_up_n  = mu_n(:,ch)' + sem_n(:,ch)';
        y_lo_n  = mu_n(:,ch)' - sem_n(:,ch)';

        valid_e = ~any(isnan([y_up_e; y_lo_e]), 1);
        valid_n = ~any(isnan([y_up_n; y_lo_n]), 1);

        if any(valid_e)
            fill(ax, [grid_t(valid_e), fliplr(grid_t(valid_e))], ...
                 [y_up_e(valid_e), fliplr(y_lo_e(valid_e))], ...
                 col_e, 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end
        if any(valid_n)
            fill(ax, [grid_t(valid_n), fliplr(grid_t(valid_n))], ...
                 [y_up_n(valid_n), fliplr(y_lo_n(valid_n))], ...
                 col_n, 'FaceAlpha', 0.20, 'EdgeColor', 'none', 'HandleVisibility', 'off');
        end

        plot(ax, grid_t, mu_e(:,ch), '-', 'Color', col_e, 'LineWidth', 1.8, 'DisplayName', 'EXO');
        plot(ax, grid_t, mu_n(:,ch), '-', 'Color', col_n, 'LineWidth', 1.8, 'DisplayName', 'NOEXO');

        title(ax, mnames{ch}, 'FontSize', 9, 'FontWeight', 'bold');
        xlabel(ax, 'Fatigue block (%)', 'FontSize', 8);
        if ch == 1 || ch == 5
            ylabel(ax, 'MDF (Hz)', 'FontSize', 8);
        end
        grid(ax, 'on'); ax.GridAlpha = 0.12;
        if ch == 1
            legend(ax, 'Location', 'southwest', 'Box', 'off', 'FontSize', 7);
        end
    end

    subplot(2, 4, 8); axis off;

    sgtitle(fig, ['MDF trajectory across the time-normalised fatigue block (0-100%)' ...
                  ' ' char(8211) ' mean ' pm ' SEM, EXO vs NOEXO'], ...
            'FontSize', 10, 'FontWeight', 'bold');
end


function save_fig(fig, fpath)
% Save figure as PNG (300 dpi). R2018b compatible (no exportgraphics).
    print(fig, fpath, '-dpng', '-r300');
end
