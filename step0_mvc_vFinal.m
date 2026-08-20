%% step0_mvc_vFinal.m
% Extracts the MVC reference by searching for the peak activation of each
% muscle across ALL MVC files of the subject, not just its own file.
%
% Improvement over v2: peak is located using a smoothed signal (robust to
% artefacts) but measured on the raw sliding RMS (recovers the true peak).
%
% Workflow:
%   1. Select all MVC CSV files for the subject
%   2. Filter all and compute peak RMS of every channel in every file
%   3. Display summary table: row = muscle, column = file
%   4. Automatically select the maximum per muscle
%   5. User can accept or manually adjust muscle by muscle
%   6. Save mvc_reference.mat compatible with step1_preprocess.m
%
% Compatibility: MATLAB R2016b+

clear; clc; close all;

%% -- CONFIGURATION ----------------------------------------------------------
Fs        = 2148;
low_cut   = 20;   high_cut = 450;
bp_order  = 4;    notch_f  = 50;
rms_win_s = 0.5;  % sliding RMS window [s]
smooth_s  = 0.5;  % smoothing window for peak localisation [s]

% Muscles in column order (col 1=AD, 2=LD, ... 7=ECR)
muscle_names = {'AD','LD','PD','UT','BB','TB','ECR'};
channel_idx  = (1:7)';
n_muscles    = numel(muscle_names);

% Filename keyword -> muscle index map
keyword_map = {
    'AD', 1;
    'LD', 2;
    'PD', 3;
    'UT', 4;
    'BB', 5;
    'TB', 6;
    'ECR', 7;
};

%% -- FILTERS ----------------------------------------------------------------
bw = (notch_f / (Fs/2)) / 35;
[bn,   an]   = butter(2,          [notch_f-bw/2, notch_f+bw/2]/(Fs/2), 'stop');
[b_bp, a_bp] = butter(bp_order/2, [low_cut, high_cut]/(Fs/2),           'bandpass');

%% -- SELECT FILES -----------------------------------------------------------
[files, fpath] = uigetfile('*.csv', ...
    'Select ALL MVC files for this subject (multi-select)', ...
    'C:\Users\gzomo\TFG', ...
    'MultiSelect', 'on');
if isequal(files, 0), error('No files selected.'); end
if ischar(files), files = {files}; end
n_files = numel(files);

% Extract short label from each filename
file_labels = cell(n_files, 1);
for fi = 1:n_files
    [~, fn, ~] = fileparts(files{fi});
    lbl = fn;
    for k = 1:size(keyword_map, 1)
        if ~isempty(strfind(fn, keyword_map{k,1}))
            lbl = keyword_map{k,1};
            break;
        end
    end
    file_labels{fi} = lbl;
end

fprintf('\n%d files loaded:\n', n_files);
for fi = 1:n_files
    fprintf('  %d. %s\n', fi, files{fi});
end

%% -- PROCESS ALL FILES ------------------------------------------------------
% peak_matrix(muscle, file) = peak RMS in V
peak_matrix = NaN(n_muscles, n_files);
rms_data    = cell(n_files, 1);   % stored for manual mode
t_data      = cell(n_files, 1);

win_samp    = round(rms_win_s * Fs);
smooth_samp = round(smooth_s  * Fs);
half_win    = round(0.5 * Fs);   % search window around detected peak

for fi = 1:n_files
    fname = files{fi};
    fprintf('\nProcessing [%d/%d]: %s ...\n', fi, n_files, fname);

    raw     = read_eu_csv(fullfile(fpath, fname));
    t       = raw(:, 1);
    emg_raw = raw(:, 2:end);

    % Filter all channels
    filtered = zeros(size(emg_raw));
    for ch = 1:n_muscles
        s = filtfilt(bn,   an,   emg_raw(:, ch));
        s = filtfilt(b_bp, a_bp, s);
        filtered(:, ch) = s;
    end

    % Vectorised sliding RMS (movmean on squared signal)
    rms_all = sqrt(movmean(filtered .^ 2, win_samp));   % [N x 7]

    % Locate peak with smoothed signal, measure on raw sliding RMS
    for ch = 1:n_muscles
        rms_smooth    = movmean(rms_all(:, ch), smooth_samp);
        [~, peak_loc] = max(rms_smooth);
        i1            = max(1, peak_loc - half_win);
        i2            = min(size(rms_all, 1), peak_loc + half_win);
        peak_matrix(ch, fi) = max(rms_all(i1:i2, ch));
    end

    rms_data{fi} = rms_all;
    t_data{fi}   = t;

    fprintf('  Peaks (uV): ');
    fprintf('%s', strjoin( ...
        arrayfun(@(ch) sprintf('%s=%.1f', muscle_names{ch}, peak_matrix(ch,fi)*1e6), ...
                 1:n_muscles, 'UniformOutput', false), '  '));
    fprintf('\n');
end

%% -- AUTOMATIC SELECTION: maximum per muscle across all files ---------------
[RMS_MVC, best_file_idx] = max(peak_matrix, [], 2);   % [7x1]

%% -- SUMMARY TABLE ----------------------------------------------------------
col_w = 9;
sep   = repmat('-', 1, 8 + n_files * col_w + 22);

fprintf('\n%s\n', sep);
fprintf('  MVC SUMMARY - peak RMS per muscle and file (uV)\n');
fprintf('  [value] = selected as MVC\n');
fprintf('%s\n', sep);

% Header
fprintf('  %-6s', 'Musc.');
for fi = 1:n_files
    fprintf('%*s', col_w, file_labels{fi});
end
fprintf('   |  MVC(uV)  Source\n');
fprintf('%s\n', sep);

% Rows
for m = 1:n_muscles
    fprintf('  %-6s', muscle_names{m});
    for fi = 1:n_files
        val_str = sprintf('%.1f', peak_matrix(m, fi) * 1e6);
        if fi == best_file_idx(m)
            val_str = ['[' val_str ']'];
        end
        fprintf('%*s', col_w, val_str);
    end
    fprintf('   |  %7.1f   %s\n', RMS_MVC(m)*1e6, file_labels{best_file_idx(m)});
end
fprintf('%s\n\n', sep);

%% -- FIGURE: grouped bar chart per muscle -----------------------------------
fig_sum = figure('Name', 'MVC Summary - peak per muscle and file', ...
                 'Position', [80 80 1000 480]);
bar_data = peak_matrix' * 1e6;   % [n_files x n_muscles]
b_h = bar(1:n_muscles, bar_data', 'grouped');
set(gca, 'XTick', 1:n_muscles, 'XTickLabel', muscle_names);
ylabel('Peak RMS (uV)');
title('Peak RMS MVC per muscle and file - [auto-selected maximum]');
legend(file_labels, 'Location', 'northeast', 'FontSize', 8);
grid on;

% Mark the selected maximum with an asterisk
hold on;
for m = 1:n_muscles
    n_groups = n_files;
    x_pos    = m + (best_file_idx(m) - (n_groups+1)/2) * (0.8/n_groups);
    y_pos    = RMS_MVC(m) * 1e6 * 1.05;
    text(x_pos, y_pos, '*', 'HorizontalAlignment', 'center', ...
         'FontSize', 10, 'Color', [0.8 0 0]);
end
hold off;

%% -- OPTIONAL MANUAL ADJUSTMENT ---------------------------------------------
fprintf('Adjust any muscle manually?\n');
fprintf('  Enter muscle number (1=AD, 2=LD, 3=PD, 4=UT, 5=BB, 6=TB, 7=ECR)\n');
fprintf('  or press Enter to accept all: ');
resp = input('', 's');

while ~isempty(resp)
    m = str2double(resp);
    if isnan(m) || m < 1 || m > n_muscles || floor(m) ~= m
        fprintf('  Invalid value. Enter a number between 1 and 7.\n');
    else
        m = round(m);
        fprintf('\n  Muscle: %s\n', muscle_names{m});
        fprintf('  Select the file to use:\n');
        for fi = 1:n_files
            fprintf('    %d. %s  (auto peak: %.1f uV)\n', fi, files{fi}, peak_matrix(m,fi)*1e6);
        end
        fi_sel = input('  File number: ');

        if isnan(fi_sel) || fi_sel < 1 || fi_sel > n_files
            fprintf('  Invalid file.\n');
        else
            fi_sel = round(fi_sel);
            % Show sliding RMS for the selected channel and file
            fig_man = figure('Name', sprintf('Manual - %s in %s', ...
                             muscle_names{m}, file_labels{fi_sel}), ...
                             'Position', [120 120 1200 400]);
            plot(t_data{fi_sel}, rms_data{fi_sel}(:,m)*1e6, 'b', 'LineWidth', 1.2);
            xlabel('Time (s)'); ylabel('Sliding RMS (uV)');
            title(sprintf('%s - click 4 points: start1 end1 start2 end2', muscle_names{m}), ...
                  'FontSize', 12);
            grid on;

            fprintf('  Click 4 points: S1 E1 S2 E2\n');
            [x_clicks, ~] = ginput(4);
            x_clicks = sort(x_clicks);

            rms_ch     = rms_data{fi_sel}(:, m);
            t_fi       = t_data{fi_sel};
            rms_manual = NaN(1, 2);
            hold on;
            colors_k = {[0.1 0.6 0.1], [0.8 0.3 0.0]};
            for k = 1:2
                [~, i1] = min(abs(t_fi - x_clicks(2*k-1)));
                [~, i2] = min(abs(t_fi - x_clicks(2*k)));
                rms_manual(k) = max(rms_ch(i1:i2));
                patch([t_fi(i1) t_fi(i2) t_fi(i2) t_fi(i1)], ...
                      [0 0 max(rms_ch)*1e6*1.1 max(rms_ch)*1e6*1.1], ...
                      colors_k{k}, 'FaceAlpha', 0.25, 'EdgeColor','none');
                text(mean(t_fi([i1 i2])), max(rms_ch(i1:i2))*1e6*1.05, ...
                     sprintf('%.1f uV', rms_manual(k)*1e6), ...
                     'HorizontalAlignment','center','FontSize',10,'Color',colors_k{k});
            end
            hold off; drawnow; pause(0.5); close(fig_man);

            new_mvc = max(rms_manual);
            fprintf('  %s -> %.1f uV (manual, from %s)\n', ...
                    muscle_names{m}, new_mvc*1e6, file_labels{fi_sel});
            RMS_MVC(m)       = new_mvc;
            best_file_idx(m) = fi_sel;
        end
    end

    fprintf('  Another muscle (number) or Enter to finish: ');
    resp = input('', 's');
end

%% -- FINAL SUMMARY ----------------------------------------------------------
fprintf('\n==============================================\n');
fprintf('  FINAL MVC\n');
fprintf('==============================================\n');
for m = 1:n_muscles
    fprintf('  %-4s  %7.1f uV  (from %s)\n', ...
            muscle_names{m}, RMS_MVC(m)*1e6, file_labels{best_file_idx(m)});
end
fprintf('==============================================\n\n');

%% -- SAVE -------------------------------------------------------------------
tok = regexp(files{1}, 'p\d+', 'match');
if isempty(tok)
    participant = input('Participant ID not detected. Enter it (e.g. p1): ', 's');
else
    participant = tok{1};
end

mvc_reference = struct();
mvc_reference.muscle_names   = muscle_names;    % {'AD','LD','PD','UT','BB','TB','ECR'}
mvc_reference.channel_idx    = channel_idx;      % [1;2;3;4;5;6;7]
mvc_reference.RMS_MVC        = RMS_MVC;          % [7x1] V  <- normalisation denominator
mvc_reference.peak_matrix    = peak_matrix;      % [7xn_files] all peaks
mvc_reference.best_file_idx  = best_file_idx;
mvc_reference.source_files   = files;
mvc_reference.file_labels    = file_labels;
mvc_reference.Fs             = Fs;
mvc_reference.rms_win_s      = rms_win_s;

mvc_dir  = 'C:\Users\gzomo\TFG\mvc';
out_file = fullfile(mvc_dir, sprintf('%s_mvc_reference.mat', participant));
save(out_file, 'mvc_reference');
fprintf('Saved -> %s\n', out_file);
