% =========================================================================
% compare_exo_noexo.m — Comparación visual EXO vs NOEXO canal a canal
% =========================================================================
% Carga dos CSVs (EXO y NOEXO) del mismo sujeto, aplica el filtrado
% estándar y muestra el envelope RMS de cada canal en la misma figura,
% EXO en azul y NOEXO en rojo.
% =========================================================================

clear; clc;

% ── Parámetros ─────────────────────────────────────────────────────────────
Fs       = 2148;
low_cut  = 20;   high_cut = 450;
bp_order = 4;    notch_f  = 50;
env_cut  = 8;    % LP envolvente

muscle_names = {'AD','LD','PD','UT','BB','TB','ECR'};
n_muscles    = numel(muscle_names);

% ── Filtros ────────────────────────────────────────────────────────────────
[b_n,  a_n]  = butter(2,          [notch_f-2, notch_f+2]/(Fs/2), 'stop');
[b_bp, a_bp] = butter(bp_order/2, [low_cut, high_cut]/(Fs/2),    'bandpass');
[b_lp, a_lp] = butter(4,          env_cut/(Fs/2),                 'low');

% ── Cargar archivos ────────────────────────────────────────────────────────
[f_exo,  p_exo]  = uigetfile('*.csv', 'Selecciona archivo EXO');
if isequal(f_exo,0), error('Cancelado.'); end

[f_noexo, p_noexo] = uigetfile('*.csv', 'Selecciona archivo NOEXO', p_exo);
if isequal(f_noexo,0), error('Cancelado.'); end

function [t, env] = load_and_filter(fpath, fname, b_n, a_n, b_bp, a_bp, b_lp, a_lp, n_muscles)
    raw = read_eu_csv(fullfile(fpath, fname));
    t   = raw(:,1);
    emg = raw(:,2:end);
    env = zeros(size(emg));
    for ch = 1:n_muscles
        sig = filtfilt(b_n,  a_n,  emg(:,ch));
        sig = filtfilt(b_bp, a_bp, sig);
        env(:,ch) = filtfilt(b_lp, a_lp, abs(sig)) * 1e6;  % µV
    end
end

fprintf('Cargando y filtrando EXO...\n');
[t_exo,  env_exo]  = load_and_filter(p_exo,  f_exo,  b_n, a_n, b_bp, a_bp, b_lp, a_lp, n_muscles);
fprintf('Cargando y filtrando NOEXO...\n');
[t_noexo, env_noexo] = load_and_filter(p_noexo, f_noexo, b_n, a_n, b_bp, a_bp, b_lp, a_lp, n_muscles);

% ── Nombre base para el título ─────────────────────────────────────────────
[~, base_exo]   = fileparts(f_exo);
[~, base_noexo] = fileparts(f_noexo);

% ── Plot ───────────────────────────────────────────────────────────────────
fig = figure('Name', 'EXO vs NOEXO — todos los canales', ...
             'Position', [50 50 1400 900]);

for ch = 1:n_muscles
    subplot(n_muscles, 1, ch);
    plot(t_exo,   env_exo(:,ch),   'b', 'LineWidth', 0.9); hold on;
    plot(t_noexo, env_noexo(:,ch), 'r', 'LineWidth', 0.9);
    ylabel('\muV', 'FontSize', 8);
    title(muscle_names{ch}, 'FontSize', 9, 'FontWeight', 'bold');
    grid on; hold off;
    if ch == 1
        legend({'EXO', 'NOEXO'}, 'Location', 'northeast', 'FontSize', 8);
    end
end

xlabel('Tiempo (s)');
sgtitle(sprintf('Envelope RMS — %s vs %s', base_exo, base_noexo), ...
        'FontSize', 12, 'FontWeight', 'bold');

% ── Guardar PNG ────────────────────────────────────────────────────────────
[~, subj] = fileparts(f_exo);
subj = regexprep(subj, '_EXO.*', '');   % extrae "p16" de "p16_EXO"

results_dir = fullfile(p_exo, 'results');
if ~exist(results_dir, 'dir'), mkdir(results_dir); end

png_path = fullfile(results_dir, [subj '_compare_exo_noexo.png']);
fig.Position = [50 50 1400 900];
exportgraphics(fig, png_path, 'Resolution', 150);
fprintf('PNG guardado: %s\n', png_path);
