% =========================================================================
% step2_mdf_window_comparison.m
% =========================================================================
% Compares three spectral window strategies for MDF estimation:
%
%   A - Current (v5): ECR onset/offset + 10% edge trim (middle 80%)
%   B - Fixed centred: ±0.5 s window around rep centre (1.0 s total)
%   C - Muscle-specific: window centred on each muscle's RMS peak
%
% Input: one or more mdf5_[name].mat files (from step2_mdf_v5.m) +
%        the corresponding pre_[name].mat files (for raw signal).
%
% Selection criterion (non-circular, per ChatGPT debate Aug 2026):
%   - % of reps with valid segment (no out-of-bounds, no NaN MDF)
%   - CV of MDF across reps per channel (lower = more stable estimation)
%   - Visual coincidence with active muscle windows (manual inspection)
%   Strategy that minimises CV AND maximises valid% wins.
%   If A ≈ B → keep A for simplicity.
%
% After computing A/B/C, the script also shows:
%   - Delta MDF comparison (as sensitivity, NOT as selection criterion)
%   - MDF rep-by-rep trajectory for each strategy (overlay plot)
%
% Compatibility: MATLAB R2016b+
% =========================================================================

clear; clc;

% -------------------------------------------------------------------------
% Parameters
% -------------------------------------------------------------------------
Fs          = 2148;
low_cut     = 20;
high_cut    = 450;
nfft        = 1024;
trim_pct_A  = 0.10;          % Strategy A: trim 10% each side
half_win_B  = round(0.5*Fs); % Strategy B: ±0.5 s from centre (1.0 s total)
peak_win_C  = round(0.4*Fs); % Strategy C: ±0.4 s around muscle RMS peak

muscle_names = {'AD','LD','PD','UT','BB','TB','ECR'};
n_muscles    = numel(muscle_names);

% -------------------------------------------------------------------------
% Select mdf5 .mat files (can select multiple)
% -------------------------------------------------------------------------
[fnames, fpath] = uigetfile('mdf5_*.mat', ...
    'Select one or more mdf5_*.mat files', 'C:\Users\gzomo\TFG\gente', ...
    'MultiSelect', 'on');
if isequal(fnames, 0), error('No files selected.'); end
if ischar(fnames), fnames = {fnames}; end

% -------------------------------------------------------------------------
% Locate corresponding pre_*.mat files
% -------------------------------------------------------------------------
pre_dir = fpath;
if contains(fpath, 'results')
    pre_dir = fileparts(strtrim(fpath));
end

fprintf('\nFound %d mdf5 file(s). Looking for pre_*.mat in: %s\n', numel(fnames), pre_dir);

% -------------------------------------------------------------------------
% Storage for group-level comparison
% -------------------------------------------------------------------------
n_files = numel(fnames);
all_cv  = nan(n_files, 3, n_muscles);   % [file, strategy, channel]
all_pct = nan(n_files, 3, n_muscles);   % % valid reps
all_delta = nan(n_files, 3, n_muscles); % delta MDF (sensitivity only)

for fi = 1:n_files
    fname_mdf5 = fnames{fi};
    [~, base] = fileparts(fname_mdf5);
    base_name = strrep(base, 'mdf5_', '');

    fprintf('\n--- Processing: %s ---\n', base_name);

    % Load mdf5 (contains onset/offset from ECR detection)
    S5 = load(fullfile(fpath, fname_mdf5));
    onset_idx  = S5.onset_idx;
    offset_idx = S5.offset_idx;
    n_reps     = S5.n_reps;
    mdf_A      = S5.mdf_reps;   % Strategy A already computed

    % Load corresponding pre_ file for raw signal
    pre_name = ['pre_' base_name '.mat'];
    pre_path = fullfile(pre_dir, pre_name);
    if ~exist(pre_path, 'file')
        % Try inside a 'results' sibling
        pre_path = fullfile(pre_dir, '..', pre_name);
    end
    if ~exist(pre_path, 'file')
        fprintf('  WARNING: cannot find %s — skipping signal-dependent strategies B and C.\n', pre_name);
        continue;
    end

    Spre   = load(pre_path);
    bp_fat = Spre.fatigue.filtered;   % bandpass-filtered [M x 7], Volts
    t_fat  = Spre.fatigue.t;

    % Linear envelope for strategy C (muscle RMS peak detection)
    [b_lp, a_lp] = butter(4, 8/(Fs/2), 'low');
    env_fat = zeros(size(bp_fat));
    for ch = 1:n_muscles
        env_fat(:,ch) = filtfilt(b_lp, a_lp, abs(bp_fat(:,ch)));
    end

    % ------------------------------------------------------------------
    % Strategy B: fixed window ±half_win_B around rep centre
    % ------------------------------------------------------------------
    mdf_B = nan(n_reps, n_muscles);
    for k = 1:n_reps
        rep_center = round((onset_idx(k) + offset_idx(k)) / 2);
        i1 = max(1, rep_center - half_win_B);
        i2 = min(size(bp_fat,1), rep_center + half_win_B);
        seg_idx = i1 : i2;
        if numel(seg_idx) < round(0.5*Fs)
            continue;   % segment too short
        end
        for ch = 1:n_muscles
            seg = bp_fat(seg_idx, ch);
            seg = seg - mean(seg);
            [pxx, f] = pburg(seg, 3, nfft, Fs);
            mask  = f >= low_cut & f <= high_cut;
            cum_p = cumsum(pxx(mask));
            f_m   = f(mask);
            mdf_B(k,ch) = interp1(cum_p, f_m, cum_p(end)/2, 'linear', 'extrap');
        end
    end

    % ------------------------------------------------------------------
    % Strategy C: window centred on each muscle's envelope peak
    % ------------------------------------------------------------------
    mdf_C = nan(n_reps, n_muscles);
    for k = 1:n_reps
        rep_range = onset_idx(k) : offset_idx(k);   % full rep range
        for ch = 1:n_muscles
            [~, pk_local] = max(env_fat(rep_range, ch));
            pk_idx = rep_range(1) + pk_local - 1;
            i1 = max(1, pk_idx - peak_win_C);
            i2 = min(size(bp_fat,1), pk_idx + peak_win_C);
            seg_idx = i1 : i2;
            if numel(seg_idx) < round(0.3*Fs)
                continue;
            end
            seg = bp_fat(seg_idx, ch);
            seg = seg - mean(seg);
            [pxx, f] = pburg(seg, 3, nfft, Fs);
            mask  = f >= low_cut & f <= high_cut;
            cum_p = cumsum(pxx(mask));
            f_m   = f(mask);
            mdf_C(k,ch) = interp1(cum_p, f_m, cum_p(end)/2, 'linear', 'extrap');
        end
    end

    % ------------------------------------------------------------------
    % Quality metrics: CV and % valid per strategy and channel
    % ------------------------------------------------------------------
    strat_labels = {'A (ECR+trim10%)', 'B (fixed 1s)', 'C (muscle peak)'};
    mdf_all = {mdf_A, mdf_B, mdf_C};

    fprintf('\n  CV of MDF across reps (lower = more stable estimation):\n');
    fprintf('  %-5s  %16s  %16s  %16s\n', 'Musc.', strat_labels{1}, strat_labels{2}, strat_labels{3});
    for ch = 1:n_muscles
        cv_row = nan(1,3);
        pv_row = nan(1,3);
        for s = 1:3
            y = mdf_all{s}(:,ch);
            valid = ~isnan(y);
            pv_row(s) = 100 * mean(valid);
            if sum(valid) >= 3
                cv_row(s) = std(y(valid)) / mean(y(valid)) * 100;
            end
            all_cv(fi, s, ch)  = cv_row(s);
            all_pct(fi, s, ch) = pv_row(s);
        end
        fprintf('  %-5s  %6.1f%% (CV%5.1f)  %6.1f%% (CV%5.1f)  %6.1f%% (CV%5.1f)\n', ...
                muscle_names{ch}, pv_row(1), cv_row(1), pv_row(2), cv_row(2), pv_row(3), cv_row(3));
    end

    % ------------------------------------------------------------------
    % Delta MDF per strategy (sensitivity analysis — NOT selection criterion)
    % ------------------------------------------------------------------
    n_edge = min(5, floor(n_reps/4));
    fprintf('\n  Delta MDF per strategy (sensitivity — first %d vs last %d reps):\n', n_edge, n_edge);
    fprintf('  %-5s  %20s  %20s  %20s\n', 'Musc.', strat_labels{1}, strat_labels{2}, strat_labels{3});
    for ch = 1:n_muscles
        d_row = nan(1,3);
        for s = 1:3
            y = mdf_all{s}(:,ch);
            vi = find(~isnan(y));
            if numel(vi) >= 2*n_edge
                d_row(s) = mean(y(vi(end-n_edge+1:end))) - mean(y(vi(1:n_edge)));
                all_delta(fi, s, ch) = d_row(s);
            end
        end
        fprintf('  %-5s  %+20.1f  %+20.1f  %+20.1f\n', muscle_names{ch}, d_row(1), d_row(2), d_row(3));
    end

    % ------------------------------------------------------------------
    % Plot: MDF trajectories for all 3 strategies (one figure per file)
    % ------------------------------------------------------------------
    t_rep = S5.t_rep;
    fig = figure('Name', ['A/B/C comparison: ' base_name], 'Position', [60 30 1200 900]);
    col_abc = {[0.2 0.5 0.85], [0.85 0.4 0.1], [0.1 0.7 0.3]};
    for ch = 1:n_muscles
        subplot(n_muscles, 1, ch);
        hold on;
        for s = 1:3
            y = mdf_all{s}(:,ch);
            plot(t_rep, y, 'o-', 'Color', col_abc{s}, 'MarkerSize', 3, ...
                 'LineWidth', 0.8, 'DisplayName', strat_labels{s});
        end
        ylabel('Hz', 'FontSize', 7);
        title(muscle_names{ch}, 'FontSize', 8, 'FontWeight', 'bold');
        grid on; hold off;
        set(gca, 'FontSize', 7);
        if ch == 1, legend('Location', 'northeast', 'FontSize', 7); end
    end
    xlabel('Relative time (s)');
    sgtitle(sprintf('MDF strategy comparison A/B/C — %s', base_name), 'FontSize', 10);

    % Save figure
    results_dir = fullfile(fpath);
    png_path = fullfile(results_dir, ['abc_' base_name '.png']);
    fig.Position = [60 30 1200 900];
    print(fig, png_path, '-dpng', '-r150');
    fprintf('\n  Figure saved: %s\n', png_path);
end

% =========================================================================
% Group-level summary
% =========================================================================
if n_files > 1
    fprintf('\n\n=== GROUP SUMMARY: mean CV across all files ===\n');
    fprintf('  (lower CV = more stable MDF estimation)\n\n');
    strat_labels = {'A (ECR+trim)', 'B (fixed 1s)', 'C (muscle pk)'};
    fprintf('  %-5s  %14s  %14s  %14s\n', 'Musc.', strat_labels{1}, strat_labels{2}, strat_labels{3});
    for ch = 1:n_muscles
        cv_mean = squeeze(nanmean(all_cv(:,:,ch), 1));
        pct_mean = squeeze(nanmean(all_pct(:,:,ch), 1));
        fprintf('  %-5s  %6.1f%%(CV%4.1f)  %6.1f%%(CV%4.1f)  %6.1f%%(CV%4.1f)\n', ...
                muscle_names{ch}, pct_mean(1), cv_mean(1), ...
                pct_mean(2), cv_mean(2), pct_mean(3), cv_mean(3));
    end

    % Winner per channel
    fprintf('\n  Recommended strategy by channel (lowest mean CV):\n');
    for ch = 1:n_muscles
        cv_row = squeeze(nanmean(all_cv(:,:,ch), 1));
        [~, best] = min(cv_row);
        fprintf('  %-5s  → %s  (CV %.1f%%)\n', muscle_names{ch}, strat_labels{best}, cv_row(best));
    end

    fprintf('\n  NOTE: Use CV as primary metric. If A and B are within ~1%% CV,\n');
    fprintf('        keep A for simplicity and reproducibility.\n');
    fprintf('        Delta MDF comparison above is sensitivity-only — do NOT\n');
    fprintf('        choose strategy based on which produces stronger fatigue signal.\n');
end

fprintf('\nDone.\n');
