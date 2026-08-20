 %% step0_preprocessing.m
% Loads a Delsys CSV (EXO or NOEXO), filters the signal, and lets the user
% manually mark the three experimental blocks.
%
% Protocol structure:
%   Baseline (5 reps) → Rest → Fatigue (reps to failure) → Rest → Post-fatigue (5 reps)
%
% User interaction:
%   6 clicks on the full signal:
%     click 1 → baseline start
%     click 2 → baseline end
%     click 3 → fatigue start
%     click 4 → fatigue end
%     click 5 → post-fatigue start
%     click 6 → post-fatigue end
%
% Output (.mat):
%   emg_struct
%     .source_file      char
%     .condition        'EXO' | 'NOEXO'
%     .Fs               2148
%     .muscle_names     cell{7,1}
%     .channel_idx      double[7]
%     .time             full time vector
%     .filtered         [N × 7]  bandpass-filtered signal
%     .envelope         [N × 7]  linear envelope
%     .blocks.baseline  struct   .idx [start end]  .time  .filtered  .envelope
%     .blocks.fatigue   struct   .idx  .time  .filtered  .envelope
%                                .windows  struct array  (for MDF)
%     .blocks.postfat   struct   .idx  .time  .filtered  .envelope

clear; clc; close all;

%% ── 0. CONFIGURATION ────────────────────────────────────────────────────

data_dir   = 'C:\Users\gzomo\TFG\raw';
output_dir = 'C:\Users\gzomo\TFG\processed';

% Signal parameters
Fs       = 2148;
low_cut  = 20;
high_cut = 450;
env_cut  = 8;
notch_f  = 50;

% Muscle / channel mapping (col index 1–7 in CSV after time column)
%   col 1 → sensor 1  Anterior Deltoid
%   col 2 → sensor 2  Lateral Deltoid
%   col 3 → sensor 4  Posterior Deltoid  (sensor 3 absent)
%   col 4 → sensor 5  Upper Trapezius
%   col 5 → sensor 6  Biceps Brachii
%   col 6 → sensor 7  Triceps Brachii
%   col 7 → sensor 8  ECR
muscle_names = {'Anterior Deltoid'; 'Lateral Deltoid'; 'Posterior Deltoid'; ...
                'Upper Trapezius'; 'Biceps Brachii'; 'Triceps Brachii'; 'ECR'};
channel_idx  = (1:7)';   % columns 1–7 already in order

% Fatigue block windowing (for MDF)
win_sec  = 5;    % window length [s]  — 5 s gives good spectral resolution at 2148 Hz
overlap  = 0.5;  % fractional overlap (0.5 = 50%)

%% ── 1. SELECT FILE AND DETECT CONDITION ─────────────────────────────────

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

[input_file, input_path] = uigetfile('*.csv', 'Select EXO or NOEXO CSV', data_dir);
if isequal(input_file, 0), error('No file selected.'); end
file_path = fullfile(input_path, input_file);   % use path from dialog, not hardcoded dir

% Detect condition from filename
if ~isempty(regexpi(input_file, 'NOEXO'))
    condition = 'NOEXO';
elseif ~isempty(regexpi(input_file, 'EXO'))
    condition = 'EXO';
else
    condition = 'UNKNOWN';
    warning('Could not detect condition from filename. Set manually if needed.');
end
fprintf('File: %s  |  Condition: %s\n', input_file, condition);

%% ── 2. LOAD CSV ─────────────────────────────────────────────────────────
% Delsys exports European format: delimiter = ';', decimal separator = ','
% textscan handles DecimalSeparator reliably across all MATLAB versions.

fprintf('Loading...\n');
raw     = read_eu_csv(file_path);
t       = raw(:, 1);
emg_raw = raw(:, 2:end);   % [N × 7]
N       = size(emg_raw, 1);
fprintf('  %d samples  |  %.1f s  |  %d channels\n', N, t(end), size(emg_raw, 2));

%% ── 3. FILTER ───────────────────────────────────────────────────────────

fprintf('Filtering...\n');

% Notch 50 Hz
bw = (notch_f / (Fs/2)) / 35;
[bn, an] = butter(2, [notch_f-bw/2, notch_f+bw/2] / (Fs/2), 'stop');

% Bandpass 20–450 Hz (order 4 applied as 2nd-order zerophase = filtfilt)
[b_bp, a_bp] = butter(2, [low_cut, high_cut] / (Fs/2), 'bandpass');

% Envelope LP 8 Hz
[b_env, a_env] = butter(4, env_cut / (Fs/2), 'low');

filtered = zeros(N, 7);
envelope = zeros(N, 7);
for ch = 1:7
    s             = filtfilt(bn,   an,   emg_raw(:, ch));
    s             = filtfilt(b_bp, a_bp, s);
    filtered(:,ch) = s;
    envelope(:,ch) = filtfilt(b_env, a_env, abs(s));
end
fprintf('Done.\n');

%% ── 4. MANUAL BLOCK SELECTION ───────────────────────────────────────────
% Show all 7 channels (offset for clarity) and ask for 6 clicks.

colors = lines(7);
offset_scale = max(envelope(:)) * 1.5;   % vertical spacing

fig = figure('Name', sprintf('Block selection — %s | %s', input_file, condition), ...
             'NumberTitle', 'off', 'Position', [50 50 1400 650]);
hold on;

for ch = 1:7
    plot(t, envelope(:, ch)*1e6 + (ch-1)*offset_scale*1e6, ...
         'Color', colors(ch,:), 'LineWidth', 0.8);
    text(t(end)*1.002, (ch-1)*offset_scale*1e6, muscle_names{ch}, ...
         'FontSize', 8, 'Color', colors(ch,:));
end

xlabel('Time (s)');
ylabel('Envelope (µV, offset per channel)');
title({'Click 6 points in order:'; ...
       '1=Baseline start  2=Baseline end  3=Fatigue start  4=Fatigue end  5=Post-fat start  6=Post-fat end'}, ...
       'FontSize', 12);
grid on;

% Draw vertical guide lines as user clicks
block_labels = {'Baseline↑', 'Baseline↓', 'Fatigue↑', 'Fatigue↓', 'PostFat↑', 'PostFat↓'};
block_colors = {'g','g','r','r','b','b'};

clicks = zeros(6, 1);
for k = 1:6
    title({sprintf('Click %d/6: %s', k, block_labels{k}); ...
           '(zoom/pan first if needed, then click)'}, 'FontSize', 12);
    [xc, ~] = ginput(1);
    clicks(k) = xc;
    xline(xc, block_colors{k}, block_labels{k}, 'LineWidth', 1.5, 'LabelVerticalAlignment', 'top');
    drawnow;
end

% Convert time clicks → sample indices
idx = zeros(6, 1);
for k = 1:6
    [~, idx(k)] = min(abs(t - clicks(k)));
end

% Validate order
if ~issorted(idx)
    warning('Clicks are not in ascending time order — sorting automatically.');
    idx = sort(idx);
end

i_base_s  = idx(1);  i_base_e  = idx(2);
i_fat_s   = idx(3);  i_fat_e   = idx(4);
i_post_s  = idx(5);  i_post_e  = idx(6);

% Shade the three blocks
y_top = (7 * offset_scale * 1e6) * 1.05;
patch([t(i_base_s) t(i_base_e) t(i_base_e) t(i_base_s)], [0 0 y_top y_top], ...
      'g', 'FaceAlpha', 0.12, 'EdgeColor', 'none');
patch([t(i_fat_s) t(i_fat_e) t(i_fat_e) t(i_fat_s)], [0 0 y_top y_top], ...
      'r', 'FaceAlpha', 0.12, 'EdgeColor', 'none');
patch([t(i_post_s) t(i_post_e) t(i_post_e) t(i_post_s)], [0 0 y_top y_top], ...
      'b', 'FaceAlpha', 0.12, 'EdgeColor', 'none');

title('Block selection complete — close figure to continue', 'FontSize', 12, 'Color', [0 0.5 0]);
drawnow;

fprintf('\nBlocks selected:\n');
fprintf('  Baseline:     t = %.1f – %.1f s  (%d samples)\n', t(i_base_s), t(i_base_e), i_base_e-i_base_s+1);
fprintf('  Fatigue:      t = %.1f – %.1f s  (%d samples)\n', t(i_fat_s),  t(i_fat_e),  i_fat_e-i_fat_s+1);
fprintf('  Post-fatigue: t = %.1f – %.1f s  (%d samples)\n', t(i_post_s), t(i_post_e), i_post_e-i_post_s+1);

pause(1);  % let user inspect before close
close(fig);

%% ── 5. SEGMENT BLOCKS ───────────────────────────────────────────────────

blocks = struct();

% Baseline
blocks.baseline.idx      = [i_base_s, i_base_e];
blocks.baseline.time     = t(i_base_s:i_base_e);
blocks.baseline.filtered = filtered(i_base_s:i_base_e, :);
blocks.baseline.envelope = envelope(i_base_s:i_base_e, :);

% Post-fatigue
blocks.postfat.idx      = [i_post_s, i_post_e];
blocks.postfat.time     = t(i_post_s:i_post_e);
blocks.postfat.filtered = filtered(i_post_s:i_post_e, :);
blocks.postfat.envelope = envelope(i_post_s:i_post_e, :);

% Fatigue block (full)
blocks.fatigue.idx      = [i_fat_s, i_fat_e];
blocks.fatigue.time     = t(i_fat_s:i_fat_e);
blocks.fatigue.filtered = filtered(i_fat_s:i_fat_e, :);
blocks.fatigue.envelope = envelope(i_fat_s:i_fat_e, :);

%% ── 6. WINDOW FATIGUE BLOCK FOR MDF ─────────────────────────────────────
% Non-overlapping windows of win_sec seconds (overlap parameter available
% but left at 0 for simplicity — change if needed).

win_samp  = round(win_sec * Fs);
step_samp = round(win_samp * (1 - overlap));
fat_sig   = blocks.fatigue.filtered;   % [M × 7]
M         = size(fat_sig, 1);

win_starts = 1 : step_samp : (M - win_samp + 1);
num_wins   = numel(win_starts);

windows(num_wins) = struct('t_center', [], 'filtered', []);
for w = 1:num_wins
    ws = win_starts(w);
    we = ws + win_samp - 1;
    windows(w).t_center  = t(i_fat_s + round((ws+we)/2) - 1);
    windows(w).filtered  = fat_sig(ws:we, :);   % [win_samp × 7]
end

blocks.fatigue.windows    = windows;
blocks.fatigue.win_sec    = win_sec;
blocks.fatigue.win_overlap = overlap;

fprintf('  Fatigue block split into %d windows of %.0f s (%.0f%% overlap)\n', ...
        num_wins, win_sec, overlap*100);

%% ── 7. SAVE ─────────────────────────────────────────────────────────────

[~, stem, ~] = fileparts(input_file);

emg_struct = struct();
emg_struct.source_file  = input_file;
emg_struct.condition    = condition;
emg_struct.Fs           = Fs;
emg_struct.muscle_names = muscle_names;
emg_struct.channel_idx  = channel_idx;
emg_struct.time         = t;
emg_struct.filtered     = filtered;
emg_struct.envelope     = envelope;
emg_struct.blocks       = blocks;

out_file = fullfile(output_dir, ['processed_' stem '.mat']);
save(out_file, 'emg_struct');
fprintf('\nSaved → %s\n', out_file);
