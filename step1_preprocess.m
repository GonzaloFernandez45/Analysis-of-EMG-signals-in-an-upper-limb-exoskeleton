% =========================================================================
% step1_preprocess.m - Loading, filtering, block selection and RMS %MVC
% =========================================================================
% Single preprocessing step. For each CSV file (EXO or NOEXO):
%
%   1. Load the raw CSV (Format A or B, auto-detected)
%   2. Apply standard filtering (notch 50Hz + bandpass 20-450Hz)
%   3. User selects the 3 blocks with 6 clicks on the signal:
%        clicks 1-2 -> baseline (start and end)
%        clicks 3-4 -> fatigue block (start and end)
%        clicks 5-6 -> post-fatigue (start and end)
%   4. Load mvc_reference.mat and compute RMS %MVC for baseline and post-fatigue
%   5. Save everything to results/pre_[name].mat
%
% Output - pre_[name].mat:
%   .fname            source CSV filename
%   .condition        'EXO' | 'NOEXO'
%   .Fs               2148
%   .muscle_names     {'AD','LD','PD','UT','BB','TB','ECR'}
%   .t                full time vector [N x 1] (s)
%   .filtered         full filtered signal [N x 7] (V)
%   .envelope         full envelope [N x 7] (V)
%   .baseline         struct: .idx .t .filtered .envelope
%   .fatigue          struct: .idx .t .filtered .envelope
%   .postfat          struct: .idx .t .filtered .envelope
%   .rms_baseline     [1 x 7] uV
%   .rms_postfat      [1 x 7] uV
%   .mvc_uv           [1 x 7] uV
%   .rms_base_mvc     [1 x 7] %MVC
%   .rms_post_mvc     [1 x 7] %MVC
%
% Compatibility: MATLAB R2016b+
% =========================================================================

clear; clc;

% -- Parameters -------------------------------------------------------------
Fs       = 2148;
low_cut  = 20;   high_cut = 450;
bp_order = 4;    notch_f  = 50;
env_cut  = 8;

muscle_names = {'AD','LD','PD','UT','BB','TB','ECR'};
n_muscles    = numel(muscle_names);

% -- Filters ----------------------------------------------------------------
[b_n,  a_n]  = butter(2,          [notch_f-2, notch_f+2]/(Fs/2), 'stop');
[b_bp, a_bp] = butter(bp_order/2, [low_cut, high_cut]/(Fs/2),    'bandpass');
[b_lp, a_lp] = butter(4,          env_cut/(Fs/2),                 'low');

% -- Load CSV ---------------------------------------------------------------
[fname, fpath] = uigetfile('*.csv', 'Select EMG file (EXO or NOEXO)', 'C:\Users\gzomo\TFG');
if isequal(fname, 0), error('No file selected.'); end

fprintf('Loading %s ...\n', fname);
raw     = read_eu_csv(fullfile(fpath, fname));
t       = raw(:, 1);
emg_raw = raw(:, 2:end);
fprintf('  %d samples | %.1f s | %d channels\n', numel(t), t(end), size(emg_raw,2));

% Detect condition from filename
if ~isempty(regexpi(fname, 'NOEXO'))
    condition = 'NOEXO';
elseif ~isempty(regexpi(fname, 'EXO'))
    condition = 'EXO';
else
    condition = input('Condition not detected. Enter EXO or NOEXO: ', 's');
end
fprintf('  Condition: %s\n', condition);

% -- Filtering --------------------------------------------------------------
fprintf('Filtering...\n');
filtered = zeros(size(emg_raw));
envelope = zeros(size(emg_raw));
for ch = 1:n_muscles
    sig = filtfilt(b_n,  a_n,  emg_raw(:,ch));
    sig = filtfilt(b_bp, a_bp, sig);
    filtered(:,ch) = sig;
    envelope(:,ch) = filtfilt(b_lp, a_lp, abs(sig));
end

% -- Global activity signal for visualisation -------------------------------
env_norm = nan(size(envelope));
for ch = 1:n_muscles
    e   = envelope(:,ch);
    rng = max(e) - min(e);
    if rng > 1e-10
        env_norm(:,ch) = (e - min(e)) / rng;
    end
end
global_act = nanmean(env_norm, 2);

% -- Block selection (6 clicks) ---------------------------------------------
fig = figure('Name', sprintf('Block selection %s', fname), ...
             'Position', [50 50 1400 600]);
plot(t, global_act, 'k', 'LineWidth', 1.2);
xlabel('Time (s)'); ylabel('Global activity (norm.)');
grid on; hold on;

labels      = {'Baseline START','Baseline END','Fatigue START','Fatigue END', ...
               'Post-fatigue START','Post-fatigue END'};
line_colors = {'b','b','r','r','g','g'};
clicks      = zeros(6,1);

for k = 1:6
    title(sprintf('Click %d/6: %s', k, labels{k}), 'FontSize', 12);
    [xc, ~] = ginput(1);
    clicks(k) = xc;
    xline(xc, line_colors{k}, labels{k}, 'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
    drawnow;
end

% Convert time to sample index
idx6 = zeros(6,1);
for k = 1:6
    [~, idx6(k)] = min(abs(t - clicks(k)));
end
if ~issorted(idx6)
    warning('Clicks out of order sorting automatically.');
    idx6 = sort(idx6);
end

i_base_s = idx6(1); i_base_e = idx6(2);
i_fat_s  = idx6(3); i_fat_e  = idx6(4);
i_post_s = idx6(5); i_post_e = idx6(6);

% Shade the three blocks
y_top = 1.05;
patch([t(i_base_s) t(i_base_e) t(i_base_e) t(i_base_s)], [0 0 y_top y_top], ...
      'b', 'FaceAlpha', 0.12, 'EdgeColor', 'none');
patch([t(i_fat_s)  t(i_fat_e)  t(i_fat_e)  t(i_fat_s)],  [0 0 y_top y_top], ...
      'r', 'FaceAlpha', 0.12, 'EdgeColor', 'none');
patch([t(i_post_s) t(i_post_e) t(i_post_e) t(i_post_s)], [0 0 y_top y_top], ...
      'g', 'FaceAlpha', 0.12, 'EdgeColor', 'none');
title('Selection complete - close figure to continue', ...
      'FontSize', 12, 'Color', [0 0.5 0]);
drawnow; pause(1); close(fig);

fprintf('\nBlocks selected:\n');
fprintf('  Baseline:     %.1f - %.1f s\n', t(i_base_s), t(i_base_e));
fprintf('  Fatigue:      %.1f - %.1f s  (%.1f s)\n', t(i_fat_s), t(i_fat_e), t(i_fat_e)-t(i_fat_s));
fprintf('  Post-fatigue: %.1f - %.1f s\n', t(i_post_s), t(i_post_e));

% -- Segment blocks ---------------------------------------------------------
baseline.idx      = [i_base_s, i_base_e];
baseline.t        = t(i_base_s:i_base_e);
baseline.filtered = filtered(i_base_s:i_base_e, :);
baseline.envelope = envelope(i_base_s:i_base_e, :);

fatigue.idx       = [i_fat_s, i_fat_e];
fatigue.t         = t(i_fat_s:i_fat_e);
fatigue.filtered  = filtered(i_fat_s:i_fat_e, :);
fatigue.envelope  = envelope(i_fat_s:i_fat_e, :);

postfat.idx       = [i_post_s, i_post_e];
postfat.t         = t(i_post_s:i_post_e);
postfat.filtered  = filtered(i_post_s:i_post_e, :);
postfat.envelope  = envelope(i_post_s:i_post_e, :);

% -- RMS per block ----------------------------------------------------------
rms_baseline = sqrt(mean(baseline.filtered .^ 2)) * 1e6;   % [1 x 7] uV
rms_postfat  = sqrt(mean(postfat.filtered  .^ 2)) * 1e6;

% -- MVC normalisation ------------------------------------------------------
[mvc_fname, mvc_fpath] = uigetfile('*.mat', 'Select MVC reference file (.mat)', 'C:\Users\gzomo\TFG\mvc');
if isequal(mvc_fname, 0), error('No MVC file selected.'); end

mvc_data = load(fullfile(mvc_fpath, mvc_fname));
mvc_ref  = mvc_data.mvc_reference;

ch_idx        = double(mvc_ref.channel_idx);
[~, sort_ord] = sort(ch_idx);
mvc_uv        = mvc_ref.RMS_MVC(sort_ord)' * 1e6;   % [1 x 7] uV

fprintf('\nMVC loaded: %s\n', mvc_fname);
fprintf('  %-5s  %8s\n', 'Musc.', 'MVC(uV)');
for ch = 1:n_muscles
    fprintf('  %-5s  %8.2f\n', muscle_names{ch}, mvc_uv(ch));
end

rms_base_mvc = (rms_baseline ./ mvc_uv) * 100;   % %MVC
rms_post_mvc = (rms_postfat  ./ mvc_uv) * 100;

fprintf('\n--- RMS %%MVC ---\n');
fprintf('  %-5s  %10s  %12s\n', 'Musc.', 'Baseline', 'Post-fatigue');
for ch = 1:n_muscles
    fprintf('  %-5s  %9.1f%%  %11.1f%%\n', muscle_names{ch}, rms_base_mvc(ch), rms_post_mvc(ch));
end

% -- Bar chart %MVC ---------------------------------------------------------
[~, base_name] = fileparts(fname);
fig_rms = figure('Name', ['RMS %MVC' base_name]);
bar_data = [rms_base_mvc; rms_post_mvc]';
b = bar(1:n_muscles, bar_data, 'grouped');
b(1).FaceColor = [0.2 0.5 0.8];
b(2).FaceColor = [0.8 0.3 0.3];
set(gca, 'XTickLabel', muscle_names, 'XTick', 1:n_muscles);
ylabel('%MVC'); grid on;
legend({'Baseline','Post-fatigue'}, 'Location', 'northeast');
title(['RMS by block ' base_name]);

% -- Save -------------------------------------------------------------------
results_dir = fullfile(fpath, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

mat_path = fullfile(results_dir, ['pre_' base_name '.mat']);
save(mat_path, ...
    'fname', 'condition', 'Fs', 'muscle_names', ...
    't', 'filtered', 'envelope', ...
    'baseline', 'fatigue', 'postfat', ...
    'rms_baseline', 'rms_postfat', 'mvc_uv', ...
    'rms_base_mvc', 'rms_post_mvc');
fprintf('\nMat saved: %s\n', mat_path);

png_path = fullfile(results_dir, ['pre_' base_name '_rms.png']);
fig_rms.Position = [100 100 900 500];
print(fig_rms, png_path, '-dpng', '-r300');
fprintf('PNG saved: %s\n', png_path);
