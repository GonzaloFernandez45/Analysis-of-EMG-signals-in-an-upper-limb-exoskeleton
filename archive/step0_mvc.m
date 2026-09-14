%% step0_mvc.m
% Extracts the MVC reference (RMS_MVC) from MVC recordings.
%
% Two modes (auto-detected by number of files selected):
%
%   MULTI-FILE mode (select 7 files):
%     One CSV per muscle, named with abbreviation: UT · AD · LD · PD · BB · TB · ECR
%     Only the target channel is used; detected automatically from filename.
%     Click 4 points per file: S1 E1 S2 E2 (or same contraction twice).
%
%   SINGLE-FILE mode (select 1 file):
%     All 7 muscles recorded sequentially in one CSV.
%     For each muscle: target channel shown in blue, others in gray.
%     Click 4 points: S1 E1 S2 E2.
%
% Output — mvc_reference.mat:
%   .muscle_names   cell{7,1}
%   .channel_idx    double[7]
%   .RMS_MVC        double[7]   max RMS across contractions [V]
%   .RMS_each       double[7,2] RMS per contraction [V]

clear; clc; close all;

%% ── 0. CONFIGURATION ────────────────────────────────────────────────────

data_dir   = 'C:\Users\gzomo\TFG\raw';
output_dir = 'C:\Users\gzomo\TFG\mvc';

Fs        = 2148;
low_cut   = 20;
high_cut  = 450;
bp_order  = 4;
env_cut   = 8;
notch_f   = 50;
rms_win_s = 0.5;

% ── Muscle / channel mapping ──────────────────────────────────────────────
%   col 1 → sensor 1  Anterior Deltoid
%   col 2 → sensor 2  Lateral Deltoid
%   col 3 → sensor 4  Posterior Deltoid  (sensor 3 absent)
%   col 4 → sensor 5  Upper Trapezius
%   col 5 → sensor 6  Biceps Brachii
%   col 6 → sensor 7  Triceps Brachii
%   col 7 → sensor 8  ECR

muscle_names = {'Upper Trapezius'; 'Anterior Deltoid'; 'Lateral Deltoid'; ...
                'Posterior Deltoid'; 'Biceps Brachii'; 'Triceps Brachii'; 'ECR'};
channel_idx  = [4; 1; 2; 3; 5; 6; 7];
%               UT  AD  LD  PD  BB  TB  ECR

% ── Filename keyword → muscle index ──────────────────────────────────────
keyword_map = {
    '_UT',   1;  'UT',  1;
    '_AD',   2;  'AD',  2;
    '_LD',   3;  'LD',  3;
    '_PD',   4;  'PD',  4;
    '_BB',   5;  'BB',  5;
    '_TB',   6;  'TB',  6;
    'ECR',   7;
};

num_muscles      = numel(muscle_names);
num_contractions = 2;

if ~exist(output_dir, 'dir'), mkdir(output_dir); end

%% ── 1. SELECT FILES ──────────────────────────────────────────────────────

[files, fpath] = uigetfile('*.csv', ...
    'Select 1 MVC file (single) or 7 files (one per muscle)', data_dir, ...
    'MultiSelect', 'on');
if isequal(files, 0), error('No files selected.'); end
if ischar(files), files = {files}; end

num_files  = numel(files);
multi_mode = num_files > 1;
fprintf('Mode: %s  (%d file(s) selected)\n', ...
        iif(multi_mode, 'MULTI-FILE', 'SINGLE-FILE'), num_files);

%% ── 2. FILTERS ───────────────────────────────────────────────────────────

bw = (notch_f / (Fs/2)) / 35;
[bn, an] = butter(2, [notch_f-bw/2, notch_f+bw/2] / (Fs/2), 'stop');
[b_bp, a_bp] = butter(bp_order/2, [low_cut, high_cut] / (Fs/2), 'bandpass');
[b_env,a_env] = butter(4, env_cut / (Fs/2), 'low');

filter_all = @(raw) apply_filters(raw, bn, an, b_bp, a_bp, b_env, a_env);

%% ── 3A. MULTI-FILE MODE ──────────────────────────────────────────────────

RMS_each     = NaN(num_muscles, num_contractions);
sel_idx      = zeros(num_muscles, num_contractions, 2);
source_files = cell(num_muscles, 1);

if multi_mode

    for fi = 1:num_files
        fname     = files{fi};
        file_path = fullfile(fpath, fname);

        % Detect muscle
        muscle_idx = detect_muscle(fname, keyword_map);
        if muscle_idx == 0
            warning('Cannot detect muscle from: %s — skipping.', fname);
            continue;
        end

        ch   = channel_idx(muscle_idx);
        name = muscle_names{muscle_idx};
        fprintf('\n[%d/%d] %s  →  %s (col %d)\n', fi, num_files, fname, name, ch);
        source_files{muscle_idx} = fname;

        raw     = read_eu_csv(file_path);
        t       = raw(:,1);
        emg_raw = raw(:,2:end);
        [filtered, envelope] = filter_all(emg_raw);

        [RMS_each(muscle_idx,:), sel_idx(muscle_idx,:,:)] = ...
            pick_contractions(t, filtered, envelope, ch, name, fi, num_files, ...
                              num_contractions, rms_win_s, Fs);

        fprintf('  contrac 1: %6.1f µV  |  contrac 2: %6.1f µV\n', ...
                RMS_each(muscle_idx,1)*1e6, RMS_each(muscle_idx,2)*1e6);
    end

%% ── 3B. SINGLE-FILE MODE ─────────────────────────────────────────────────

else
    fname     = files{1};
    file_path = fullfile(fpath, fname);
    fprintf('Loading: %s\n', fname);
    source_files{1} = fname;

    raw     = read_eu_csv(file_path);
    t       = raw(:,1);
    emg_raw = raw(:,2:end);
    [filtered, envelope] = filter_all(emg_raw);

    fprintf('  %d samples  |  %.1f s\n', size(emg_raw,1), t(end));
    fprintf('For each muscle: click S1 E1 S2 E2 (4 clicks in order).\n\n');

    for m = 1:num_muscles
        ch   = channel_idx(m);
        name = muscle_names{m};

        [RMS_each(m,:), sel_idx(m,:,:)] = ...
            pick_contractions(t, filtered, envelope, ch, name, m, num_muscles, ...
                              num_contractions, rms_win_s, Fs);

        fprintf('  %-22s  contrac 1: %6.1f µV  |  contrac 2: %6.1f µV\n', ...
                name, RMS_each(m,1)*1e6, RMS_each(m,2)*1e6);
    end
end

%% ── 4. RMS_MVC ───────────────────────────────────────────────────────────

RMS_MVC = max(RMS_each, [], 2);

fprintf('\n── RMS_MVC summary ──────────────────────────────────────────\n');
for m = 1:num_muscles
    if isnan(RMS_MVC(m))
        fprintf('  %-22s  NOT PROCESSED\n', muscle_names{m});
    else
        fprintf('  %-22s  RMS_MVC = %7.2f µV\n', muscle_names{m}, RMS_MVC(m)*1e6);
    end
end

%% ── 5. SAVE ─────────────────────────────────────────────────────────────

tok = regexp(files{1}, '(p\d+)', 'tokens');
participant = tok{1}{1};

mvc_reference = struct();
mvc_reference.muscle_names   = muscle_names;
mvc_reference.channel_idx    = channel_idx;
mvc_reference.RMS_MVC        = RMS_MVC;
mvc_reference.RMS_each       = RMS_each;
mvc_reference.sel_idx        = sel_idx;
mvc_reference.source_files   = source_files;
mvc_reference.Fs             = Fs;
mvc_reference.rms_win_s      = rms_win_s;

out_file = fullfile(output_dir, sprintf('%s_mvc_reference.mat', participant));
save(out_file, 'mvc_reference');
fprintf('\nSaved → %s\n', out_file);

%% ══ LOCAL FUNCTIONS ══════════════════════════════════════════════════════

function [filtered, envelope] = apply_filters(emg_raw, bn, an, b_bp, a_bp, b_env, a_env)
    filtered = zeros(size(emg_raw));
    envelope = zeros(size(emg_raw));
    for channel_idx = 1:size(emg_raw, 2)
        signal = filtfilt(bn, an, emg_raw(:,channel_idx));
        signal = filtfilt(b_bp, a_bp, signal);
        filtered(:,channel_idx) = signal;
        envelope(:,channel_idx) = filtfilt(b_env, a_env, abs(signal));
    end
end

function muscle_idx = detect_muscle(fname, keyword_map)
    muscle_idx = 0;
    for k = 1:size(keyword_map, 1)
        if contains(fname, keyword_map{k,1})
            muscle_idx = keyword_map{k,2};
            return;
        end
    end
end

function [rms_row, idx_row] = pick_contractions(t, filtered, ~, ch, name, ...
                                                 fi, n_total, n_contrac, rms_win_s, Fs)
    % Calcular sliding window RMS sobre toda la señal del canal objetivo
    win_samp = round(rms_win_s * Fs);
    sig      = filtered(:, ch);
    n_steps  = length(sig) - win_samp + 1;
    rms_full = zeros(n_steps, 1);
    for w = 1:n_steps
        rms_full(w) = sqrt(mean(sig(w:w+win_samp-1).^2));
    end
    t_rms = t(1:n_steps);

    % Mostrar solo el canal objetivo (RMS ya calculado)
    fig = figure('Name', sprintf('MVC — %s', name), ...
                 'NumberTitle', 'off', 'Position', [80 80 1300 520]);
    plot(t_rms, rms_full * 1e6, 'b', 'LineWidth', 1.5);
    xlabel('Time (s)'); ylabel('Sliding RMS (µV)');
    title(sprintf('[%d/%d]  %s  — click: S1  E1  S2  E2', fi, n_total, name), 'FontSize', 13);
    grid on;

    % Usuario selecciona rango de cada contracción: S1 E1 S2 E2
    [x_clicks, ~] = ginput(n_contrac * 2);
    x_clicks = sort(x_clicks);

    idx_row    = zeros(n_contrac, 2);
    rms_row    = NaN(1, n_contrac);
    colors_sel = {[0.2 0.7 0.2], [0.8 0.3 0.1]};
    hold on;
    y_max = max(rms_full) * 1e6;

    for k = 1:n_contrac
        [~, i1] = min(abs(t_rms - x_clicks(2*k-1)));
        [~, i2] = min(abs(t_rms - x_clicks(2*k)));
        idx_row(k,:) = [i1, i2];
        rms_row(k)   = max(rms_full(i1:i2));   % pico RMS dentro del rango

        patch([t_rms(i1) t_rms(i2) t_rms(i2) t_rms(i1)], [0 0 y_max y_max], ...
              colors_sel{k}, 'FaceAlpha', 0.25, 'EdgeColor', 'none');
        text(mean(t_rms([i1 i2])), y_max * 0.88, ...
             sprintf('RMS=%.1f µV', rms_row(k)*1e6), ...
             'HorizontalAlignment', 'center', 'FontSize', 10, ...
             'Color', colors_sel{k}, 'FontWeight', 'bold');
    end

    drawnow; pause(0.5); close(fig);
end

function out = iif(cond, a, b)
    if cond, out = a; else, out = b; end
end

% read_eu_csv está en read_eu_csv.m (archivo separado en la misma carpeta)
