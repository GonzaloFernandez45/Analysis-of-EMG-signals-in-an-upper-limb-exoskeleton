% =========================================================================
% step1_rms.m — RMS de baseline y post-fatiga
% =========================================================================
% Carga el CSV raw, aplica el mismo filtrado que step2_mdf_v2.m,
% el usuario selecciona los bloques baseline y post-fatiga con ginput,
% calcula el RMS medio de cada bloque por canal, normaliza por MVC → %MVC,
% y guarda en results/rms1_[nombre].mat
%
% Output: rms1_p[N]_[EXO/NOEXO].mat
% =========================================================================

clear; clc;

% ── Parámetros ────────────────────────────────────────────────────────────
Fs       = 2148;
low_cut  = 20;   high_cut = 450;
bp_order = 4;    notch_f  = 50;
muscle_names = {'AD','LD','PD','UT','BB','TB','ECR'};
n_muscles    = numel(muscle_names);

% ── Cargar CSV ────────────────────────────────────────────────────────────
[fname, fpath] = uigetfile('*.csv', 'Selecciona el archivo EMG (EXO o NOEXO)');
if isequal(fname, 0), error('No se seleccionó archivo.'); end

fprintf('Cargando %s ...\n', fname);
raw     = read_eu_csv(fullfile(fpath, fname));
t_all   = raw(:, 1);
emg_all = raw(:, 2:end);
fprintf('  %d muestras | %.1f s | %d canales\n', numel(t_all), t_all(end), size(emg_all,2));

% ── Filtrado ──────────────────────────────────────────────────────────────
[b_n,  a_n]  = butter(2,          [notch_f-2, notch_f+2]/(Fs/2), 'stop');
[b_bp, a_bp] = butter(bp_order/2, [low_cut, high_cut]/(Fs/2),    'bandpass');
[b_lp, a_lp] = butter(4,          8/(Fs/2),                       'low');

emg_bp  = zeros(size(emg_all));
emg_env = zeros(size(emg_all));
for ch = 1:n_muscles
    sig = emg_all(:, ch);
    sig = filtfilt(b_n,  a_n,  sig);
    sig = filtfilt(b_bp, a_bp, sig);
    emg_bp(:, ch)  = sig;
    emg_env(:, ch) = filtfilt(b_lp, a_lp, abs(sig));
end

% ── Señal global de actividad ─────────────────────────────────────────────
env_norm = nan(size(emg_env));
for ch = 1:n_muscles
    e   = emg_env(:, ch);
    rng = max(e) - min(e);
    if rng > 1e-10
        env_norm(:, ch) = (e - min(e)) / rng;
    end
end
global_act = nanmean(env_norm, 2);

% ── Selección de bloques ──────────────────────────────────────────────────
figure('Name', 'Selección de bloques — baseline y post-fatiga');
plot(t_all, global_act, 'k', 'LineWidth', 1.2);
xlabel('Tiempo (s)'); ylabel('Actividad global (norm.)');
grid on; hold on;

% Baseline (azul)
title('Click en INICIO y FIN del bloque BASELINE (azul)');
fprintf('\nClick en INICIO del baseline...\n');
[t_base_start, ~] = ginput(1); xline(t_base_start, 'b--', 'LineWidth', 1.5);
fprintf('Click en FIN del baseline...\n');
[t_base_end, ~]   = ginput(1); xline(t_base_end,   'b--', 'LineWidth', 1.5);
fprintf('Baseline: %.1f – %.1f s (%.1f s)\n', t_base_start, t_base_end, t_base_end - t_base_start);

% Post-fatiga (rojo)
title('Click en INICIO y FIN del bloque POST-FATIGA (rojo)');
fprintf('Click en INICIO del post-fatiga...\n');
[t_post_start, ~] = ginput(1); xline(t_post_start, 'r--', 'LineWidth', 1.5);
fprintf('Click en FIN del post-fatiga...\n');
[t_post_end, ~]   = ginput(1); xline(t_post_end,   'r--', 'LineWidth', 1.5);
fprintf('Post-fatiga: %.1f – %.1f s (%.1f s)\n', t_post_start, t_post_end, t_post_end - t_post_start);
hold off;

% ── Recortar bloques ──────────────────────────────────────────────────────
idx_base = t_all >= t_base_start & t_all <= t_base_end;
idx_post = t_all >= t_post_start & t_all <= t_post_end;
bp_base  = emg_bp(idx_base, :);
bp_post  = emg_bp(idx_post, :);

% ── RMS por canal por bloque ──────────────────────────────────────────────
rms_baseline    = sqrt(mean(bp_base.^2)) * 1e6;   % µV [1x7]
rms_postfatigue = sqrt(mean(bp_post.^2)) * 1e6;   % µV [1x7]

fprintf('\n--- RMS (µV) ---\n');
fprintf('  %-5s  %10s  %12s\n', 'Musc.', 'Baseline', 'Post-fatiga');
for ch = 1:n_muscles
    fprintf('  %-5s  %10.2f  %12.2f\n', muscle_names{ch}, rms_baseline(ch), rms_postfatigue(ch));
end

% ── Normalización por MVC ─────────────────────────────────────────────────
% Cargar el archivo mvc_reference.mat generado por step0_mvc.m
[mvc_fname, mvc_fpath] = uigetfile('*.mat', 'Selecciona el archivo MVC reference (.mat)');
if isequal(mvc_fname, 0), error('No se seleccionó archivo MVC.'); end

mvc_data = load(fullfile(mvc_fpath, mvc_fname));
mvc_ref  = mvc_data.mvc_reference;

% Reordenar por channel_idx → orden AD(1),LD(2),PD(3),UT(4),BB(5),TB(6),ECR(7)
ch_idx        = double(mvc_ref.channel_idx);
[~, sort_ord] = sort(ch_idx);
mvc_uv        = mvc_ref.RMS_MVC(sort_ord)' * 1e6;   % [1x7] µV

fprintf('\nMVC cargado desde: %s\n', mvc_fname);
fprintf('  %-5s  %8s\n', 'Musc.', 'MVC (µV)');
for ch = 1:n_muscles
    fprintf('  %-5s  %8.2f\n', muscle_names{ch}, mvc_uv(ch));
end

rms_base_mvc = (rms_baseline    ./ mvc_uv) * 100;   % %MVC
rms_post_mvc = (rms_postfatigue ./ mvc_uv) * 100;   % %MVC

fprintf('\n--- %%MVC ---\n');
fprintf('  %-5s  %10s  %12s\n', 'Musc.', 'Baseline', 'Post-fatiga');
for ch = 1:n_muscles
    fprintf('  %-5s  %9.1f%%  %11.1f%%\n', muscle_names{ch}, rms_base_mvc(ch), rms_post_mvc(ch));
end

% ── Plot barras ───────────────────────────────────────────────────────────
[~, base_name, ~] = fileparts(fname);
fig_rms = figure('Name', ['RMS %%MVC — ' base_name]);
bar_data = [rms_base_mvc; rms_post_mvc]';
b = bar(1:n_muscles, bar_data, 'grouped');
b(1).FaceColor = [0.2 0.5 0.8];
b(2).FaceColor = [0.8 0.3 0.3];
set(gca, 'XTickLabel', muscle_names, 'XTick', 1:n_muscles);
ylabel('%MVC');
legend({'Baseline', 'Post-fatiga'}, 'Location', 'northeast');
title(['RMS por bloque — ' base_name]);
grid on;

% ── Guardar ───────────────────────────────────────────────────────────────
results_dir = fullfile(fpath, 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
    fprintf('\nCarpeta creada: %s\n', results_dir);
end

mat_path = fullfile(results_dir, ['rms1_' base_name '.mat']);
save(mat_path, ...
    'rms_baseline', 'rms_postfatigue', ...
    'rms_base_mvc', 'rms_post_mvc',    ...
    'mvc_uv', 'muscle_names',          ...
    't_base_start', 't_base_end',      ...
    't_post_start', 't_post_end',      ...
    'fname');
fprintf('Mat guardado: %s\n', mat_path);

png_path = fullfile(results_dir, ['rms1_' base_name '_bars.png']);
fig_rms.Position = [100 100 900 500];
exportgraphics(fig_rms, png_path, 'Resolution', 300);
fprintf('PNG guardado: %s\n', png_path);
