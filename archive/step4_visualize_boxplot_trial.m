% =========================================================================
% step4_visualize_boxplot_trial.m - Trial: boxplot version of Chapter 6 figures
% =========================================================================
% Loads step3_results.mat and plots the same four outcome blocks as
% step4_visualize.m, but as boxplots (box + median + whiskers + a dot
% for the mean) instead of bar + paired dots + connecting lines.
%
% No block-level title (that goes in the table/figure caption in the
% document instead). Significance markers per muscle:
%   *** = Bonferroni-significant (p < alpha_bonf)
%   dagger = nominally significant (p < 0.05, not Bonferroni)
%   * (endurance only) = significant at uncorrected alpha = 0.05
%
% This is a TRIAL script, separate from the frozen step4_visualize.m,
% to compare both styles before deciding which one goes in the TFG.
%
% Output: 4 PNGs (300 dpi) saved next to step3_results.mat.
%
% Compatibility: MATLAB R2016b+ (boxplot requires Statistics Toolbox)
% =========================================================================

clear; clc;

[fname, fpath] = uigetfile('*.mat', 'Select step3_results.mat', 'C:\Users\gzomo\TFG\results');
if isequal(fname, 0), error('No file selected.'); end

S = load(fullfile(fpath, fname));
muscle_names = S.muscle_names;
n_muscles    = numel(muscle_names);

% -------------------------------------------------------------------------
% Figure 1 - Baseline RMS (%MVC)
% -------------------------------------------------------------------------
plot_block(S.rms_base_noexo, S.rms_base_exo, muscle_names, 'Baseline RMS (%MVC)', ...
    logical(S.res_base.sig_bonf), logical(S.res_base.sig_nom), ...
    fullfile(fpath, 'boxplot_baseline_rms.png'));

% -------------------------------------------------------------------------
% Figure 2 - Post-fatigue RMS (%MVC)
% -------------------------------------------------------------------------
plot_block(S.rms_post_noexo, S.rms_post_exo, muscle_names, 'Post-fatigue RMS (%MVC)', ...
    logical(S.res_post.sig_bonf), logical(S.res_post.sig_nom), ...
    fullfile(fpath, 'boxplot_postfatigue_rms.png'));

% -------------------------------------------------------------------------
% Figure 3 - Delta MDF (Hz)
% -------------------------------------------------------------------------
plot_block(S.dmdf_noexo, S.dmdf_exo, muscle_names, 'Delta MDF (Hz)', ...
    logical(S.res_mdf.sig_bonf), logical(S.res_mdf.sig_nom), ...
    fullfile(fpath, 'boxplot_delta_mdf.png'));

% -------------------------------------------------------------------------
% Figure 4 - Endurance (n_reps, block duration) - uncorrected alpha = 0.05
% -------------------------------------------------------------------------
fig_end = figure('Name', 'Endurance', 'Position', [100 100 700 420]);

subplot(1,2,1);
plot_pair(S.nreps_noexo, S.nreps_exo, 'n_{reps}', false, logical(S.res_end.nreps.sig));
ylabel('n_{reps}');

subplot(1,2,2);
plot_pair(S.dur_noexo, S.dur_exo, 'Block duration (s)', false, logical(S.res_end.duration.sig));
ylabel('Duration (s)');

print(fig_end, fullfile(fpath, 'boxplot_endurance.png'), '-dpng', '-r300');
fprintf('Saved: %s\n', fullfile(fpath, 'boxplot_endurance.png'));

fprintf('\nDone. 4 boxplot PNGs saved in %s\n', fpath);

% =========================================================================
% LOCAL FUNCTIONS (must be at end of script - MATLAB R2016b+)
% =========================================================================

function plot_block(mat_noexo, mat_exo, muscle_names, ylab, sig_bonf, sig_nom, save_path)
    n_muscles = numel(muscle_names);
    fig = figure('Name', ylab, 'Position', [50 50 1500 320]);
    for ch = 1:n_muscles
        subplot(1, n_muscles, ch);
        plot_pair(mat_noexo(:,ch), mat_exo(:,ch), muscle_names{ch}, sig_bonf(ch), sig_nom(ch));
        if ch == 1, ylabel(ylab); end
    end
    print(fig, save_path, '-dpng', '-r300');
    fprintf('Saved: %s\n', save_path);
end

function plot_pair(vec_noexo, vec_exo, ttl, is_bonf, is_nom)
    data   = [vec_noexo; vec_exo];
    groups = [ones(numel(vec_noexo),1); 2*ones(numel(vec_exo),1)];
    boxplot(data, groups, 'Labels', {'NOEXO','EXO'});
    hold on;
    plot(1, mean(vec_noexo, 'omitnan'), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);
    plot(2, mean(vec_exo,   'omitnan'), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 6);

    if is_bonf || is_nom
        all_vals = [vec_noexo; vec_exo];
        y_max = max(all_vals);
        y_min = min(all_vals);
        y_rng = y_max - y_min;
        y_sig = y_max + 0.08 * y_rng;

        line([1 2], [y_sig y_sig], 'Color', 'k', 'LineWidth', 1);
        if is_bonf
            sig_txt = '***';
        else
            sig_txt = char(8224); % dagger, nominal only
        end
        text(1.5, y_sig, sig_txt, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'bottom', 'FontSize', 12);
        ylim([y_min - 0.05*y_rng, y_sig + 0.15*y_rng]);
    end

    hold off;
    title(ttl);
    grid on;
end
