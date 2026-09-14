% =========================================================================
% step2_mdf_v2.m — MDF per repetition (fatigue block)
% =========================================================================
% Computes Median Frequency (MDF) on each individual repetition of the
% fatigue block, then fits a linear regression across repetitions to
% estimate the rate of spectral compression (fatigue slope).
%
% Method: AR model order 3, Burg's method (pburg) — as per Cap. 5.
%
% Output: mdf2_p[N].mat
% =========================================================================

clear; clc;

% ── Parámetros de adquisición y filtrado ─────────────────────────────────
Fs       = 2148;    % Hz
low_cut  = 20;      % Hz — límite inferior bandpass
high_cut = 450;     % Hz — límite superior bandpass
bp_order = 4;       % orden efectivo (butter orden 2 + filtfilt = orden 4)
notch_f  = 50;      % Hz — frecuencia de red
env_cut  = 8;       % Hz — LP envolvente (solo para detección de onsets)

% Músculos y columnas en el CSV de Delsys (formato europeo: ; separador, , decimal)
% Col 1=tiempo, cols 2-8 = 7 músculos en orden (sensor 3 ausente, no ocupa columna)
%   col 2=AD, col 3=LD, col 4=PD, col 5=UT, col 6=BB, col 7=TB, col 8=ECR
muscle_names = {'AD','LD','PD','UT','BB','TB','ECR'};
n_muscles    = numel(muscle_names);

% Canal de referencia para detección de onsets:
% se elige visualmente tras ver las envolventes (ver más abajo)

% ── Cargar archivo ────────────────────────────────────────────────────────
[fname, fpath] = uigetfile('*.csv', 'Selecciona el archivo EMG (EXO o NOEXO)');
if isequal(fname, 0), error('No se seleccionó archivo.'); end

fprintf('Cargando %s ...\n', fname);
% Formato europeo Delsys: separador=';', decimal=',', 1 línea de cabecera
raw = read_eu_csv(fullfile(fpath, fname));

t_all   = raw(:, 1);       % vector de tiempo (s)
emg_all = raw(:, 2:end);   % [N x 7] señales en voltios — cols 2-8
fprintf('  %d muestras | %.1f s | %d canales\n', numel(t_all), t_all(end), size(emg_all,2));

% ── Filtrado ──────────────────────────────────────────────────────────────
% Notch 50 Hz
[b_n, a_n] = butter(2, [notch_f-2, notch_f+2]/(Fs/2), 'stop');

% Bandpass 20-450 Hz (orden 2 + filtfilt = orden 4 efectivo)
[b_bp, a_bp] = butter(bp_order/2, [low_cut, high_cut]/(Fs/2), 'bandpass');

% LP 8 Hz para envolvente (orden 4 + filtfilt = orden 8 efectivo)
[b_lp, a_lp] = butter(4, env_cut/(Fs/2), 'low');

emg_bp  = zeros(size(emg_all));   % señal bandpass — para MDF
emg_env = zeros(size(emg_all));   % envolvente — para detección onsets

for ch = 1:n_muscles
    sig = emg_all(:, ch);
    sig = filtfilt(b_n,  a_n,  sig);   % notch
    sig = filtfilt(b_bp, a_bp, sig);   % bandpass
    emg_bp(:, ch) = sig;               % guardar bandpass para MDF
    emg_env(:, ch) = filtfilt(b_lp, a_lp, abs(sig));  % rectif + LP = envolvente
end

% ── Selección del bloque de fatiga ───────────────────────────────────────
% Normalizar cada envolvente a [0,1] y promediar → señal global de actividad
env_norm = nan(size(emg_env));
for ch = 1:n_muscles
    e = emg_env(:, ch);
    rng = max(e) - min(e);
    if rng > 1e-10          % ignorar canales planos (ej. UT)
        env_norm(:, ch) = (e - min(e)) / rng;
    end
end
global_act = nanmean(env_norm, 2);  % media ignorando canales NaN

fprintf('t_all: %d muestras, rango %.2f – %.2f s\n', numel(t_all), t_all(1), t_all(end));
fprintf('global_act: %d valores, NaN: %d\n', numel(global_act), sum(isnan(global_act)));

figure('Name','Selección bloque de fatiga');
plot(t_all, global_act, 'k', 'LineWidth', 1.2);
xlabel('Tiempo (s)'); ylabel('Actividad global (norm.)');
title('Haz click en el INICIO y FIN del bloque de fatiga');
grid on;

fprintf('Haz click en INICIO del bloque de fatiga...\n');
[t_fat_start, ~] = ginput(1);
xline(t_fat_start, 'b--', 'LineWidth', 1.5);

fprintf('Haz click en FIN del bloque de fatiga...\n');
[t_fat_end, ~] = ginput(1);
xline(t_fat_end, 'r--', 'LineWidth', 1.5);

fprintf('Bloque de fatiga: %.1f s → %.1f s (%.1f s)\n', ...
    t_fat_start, t_fat_end, t_fat_end - t_fat_start);

% Recortar al bloque de fatiga
idx_fat = t_all >= t_fat_start & t_all <= t_fat_end;
t_fat   = t_all(idx_fat);
bp_fat  = emg_bp(idx_fat, :);
env_fat = emg_env(idx_fat, :);

% ── Subplot envolventes del bloque de fatiga → elegir canal de referencia ─
figure('Name','Envolventes bloque de fatiga — elige canal de referencia');
for ch = 1:n_muscles
    subplot(n_muscles, 1, ch);
    plot(t_fat, env_fat(:, ch) * 1e6, 'b', 'LineWidth', 0.8);
    ylabel('\muV');
    title(muscle_names{ch});
    grid on;
end
xlabel('Tiempo (s)');
sgtitle('Envolventes del bloque de fatiga — ¿qué canal usamos como referencia?');

ref_muscle = input('Introduce el número del canal de referencia (1=AD, 2=LD, 3=PD, 4=UT, 5=BB, 6=TB, 7=ECR): ');

% ── Sliding window RMS para detección de onsets ───────────────────────────
% Ventana 0.5s → suaviza la señal mostrando la energía local en cada instante
win_rms = round(0.5 * Fs);   % 0.5s en muestras
rms_ref = sqrt(movmean(bp_fat(:, ref_muscle).^2, win_rms)) * 1e6;  % en µV

% Mostrar para verificar
figure('Name', ['Sliding window RMS — ' muscle_names{ref_muscle}]);
plot(t_fat, rms_ref, 'b', 'LineWidth', 1);
xlabel('Tiempo (s)'); ylabel('\muV');
title(['Sliding window RMS (0.5s) — ' muscle_names{ref_muscle}]);
grid on;

% ── Detección de onsets por umbral ────────────────────────────────────────
% Umbral = prctile(rms_ref, pct_thr).
% Percentil 25 por defecto: en el protocolo 1 rep/4s (~2.5s activo, ~1.5s reposo)
% los valles ocupan ~35-40% del tiempo → el percentil 25 cae en la distribución
% de valles, no de picos. Robusto a spikes que inflarían un umbral basado en el máximo.
pct_thr = 25;
thr     = prctile(rms_ref, pct_thr);

while true
    % Binarizar y detectar flancos
    above   = rms_ref > thr;
    d_above = diff([0; above; 0]);
    onset_idx  = find(d_above ==  1);
    offset_idx = find(d_above == -1);

    % Eliminar segmentos demasiado cortos (< 1s = artefacto o ruido)
    min_dur    = round(1.0 * Fs);
    valid      = (offset_idx - onset_idx) >= min_dur;
    onset_idx  = onset_idx(valid);
    offset_idx = offset_idx(valid);
    n_reps     = numel(onset_idx);

    fprintf('\n--- Detección de repeticiones ---\n');
    fprintf('  Percentil umbral: %d  (equivale a %.1f µV) | Repeticiones detectadas: %d\n', ...
            pct_thr, thr, n_reps);
    for k = 1:n_reps
        dur_k = (offset_idx(k) - onset_idx(k)) / Fs;
        fprintf('  Rep %2d: %.2f s – %.2f s (%.2f s)\n', k, ...
            t_fat(onset_idx(k)), t_fat(offset_idx(k)), dur_k);
    end

    % Visualizar onsets detectados sobre el RMS
    fig_onsets = figure('Name', ['Onsets detectados — ' muscle_names{ref_muscle}]);
    plot(t_fat, rms_ref, 'b', 'LineWidth', 1); hold on;
    yline(thr, 'k--', sprintf('Pct %d → %.1f µV', pct_thr, thr), 'LineWidth', 1);
    for k = 1:n_reps
        xline(t_fat(onset_idx(k)),  'g', 'LineWidth', 1.2);
        xline(t_fat(offset_idx(k)), 'r', 'LineWidth', 1.0);
    end
    xlabel('Tiempo (s)'); ylabel('\muV');
    title(sprintf('Onsets (verde) y offsets (rojo) — %d reps | Pct %d (%.1f µV)', ...
          n_reps, pct_thr, thr));
    grid on; hold off;

    % Ajuste manual: el usuario introduce un percentil (0-100), no un µV
    fprintf('\n  Pista: baja el percentil → umbral más bajo → más reps detectadas.\n');
    fprintf('         Sube el percentil → umbral más alto → solo picos claros.\n');
    resp = input('  ¿Ajustar? (Enter para aceptar, o escribe nuevo percentil 0–100): ', 's');
    if isempty(resp)
        close(fig_onsets);
        break;
    else
        new_pct = str2double(resp);
        if isnan(new_pct) || new_pct < 0 || new_pct > 100
            fprintf('  Valor no válido — introduce un número entre 0 y 100.\n');
        else
            pct_thr = new_pct;
            thr     = prctile(rms_ref, pct_thr);
        end
        close(fig_onsets);
    end
end

% ── MDF por repetición (AR orden 3, método de Burg) ──────────────────────
% Según Cap. 5: pburg en cada repetición individual → 1 valor MDF por rep
mdf_reps = nan(n_reps, n_muscles);   % [n_reps x 7]
t_rep    = nan(n_reps, 1);           % tiempo central de cada rep (s, relativo al bloque)

fprintf('\nCalculando MDF por repetición...\n');
for k = 1:n_reps
    seg_idx  = onset_idx(k) : offset_idx(k);
    t_rep(k) = mean(t_fat(seg_idx)) - t_fat(1);   % tiempo central relativo al inicio del bloque

    for ch = 1:n_muscles
        seg = bp_fat(seg_idx, ch);

        % PSD por Burg (AR orden 3)
        % pburg(x, p, [], Fs) → usa nfft por defecto, devuelve frecuencias en Hz
        [pxx, f] = pburg(seg, 3, [], Fs);

        % Restringir al rango del bandpass (20–450 Hz)
        mask = f >= low_cut & f <= high_cut;
        f_m  = f(mask);
        p_m  = pxx(mask);

        % MDF: frecuencia donde la potencia acumulada cruza el 50% del total
        cum_p = cumsum(p_m);
        mdf_reps(k, ch) = interp1(cum_p, f_m, cum_p(end)/2, 'linear', 'extrap');
    end
end
fprintf('  Hecho.\n');

% ── Regresión lineal MDF vs tiempo ───────────────────────────────────────
slopes = nan(1, n_muscles);   % pendiente (Hz/s)
r2     = nan(1, n_muscles);   % coeficiente de determinación

for ch = 1:n_muscles
    y     = mdf_reps(:, ch);
    valid = ~isnan(y);
    if sum(valid) >= 3
        p      = polyfit(t_rep(valid), y(valid), 1);
        slopes(ch) = p(1);                          % Hz/s (negativo = fatiga)
        y_fit  = polyval(p, t_rep(valid));
        ss_res = sum((y(valid) - y_fit).^2);
        ss_tot = sum((y(valid) - mean(y(valid))).^2);
        r2(ch) = 1 - ss_res / ss_tot;
    end
end

% ── Resultados en consola ─────────────────────────────────────────────────
fprintf('\n--- MDF slopes (Hz/s) ---\n');
fprintf('  %-5s  %8s  %6s\n', 'Musc.', 'Slope', 'R²');
for ch = 1:n_muscles
    fprintf('  %-5s  %+8.4f  %6.3f\n', muscle_names{ch}, slopes(ch), r2(ch));
end

% ── Plot MDF rep a rep + recta de regresión ───────────────────────────────
t_line = linspace(t_rep(1), t_rep(end), 100);
fig_mdf = figure('Name', 'MDF por repetición — todos los canales');
for ch = 1:n_muscles
    subplot(n_muscles, 1, ch);
    plot(t_rep, mdf_reps(:, ch), 'bo', 'MarkerFaceColor', 'b', 'MarkerSize', 4);
    hold on;
    if ~isnan(slopes(ch))
        p_ch = polyfit(t_rep(~isnan(mdf_reps(:,ch))), ...
                       mdf_reps(~isnan(mdf_reps(:,ch)), ch), 1);
        plot(t_line, polyval(p_ch, t_line), 'r-', 'LineWidth', 1.2);
    end
    ylabel('Hz');
    title(sprintf('%s — slope: %+.3f Hz/s, R²=%.2f', ...
        muscle_names{ch}, slopes(ch), r2(ch)));
    grid on; hold off;
end
xlabel('Tiempo relativo (s)');
sgtitle('MDF por repetición (AR-3 Burg) — bloque de fatiga');

% ── Guardar resultados ────────────────────────────────────────────────────
% Nombre base = nombre del CSV sin extensión (ej. "p5_NOEXO")
[~, base_name, ~] = fileparts(fname);

% Carpeta results/ junto al CSV — crearla si no existe
results_dir = fullfile(fpath, 'results');
if ~exist(results_dir, 'dir')
    mkdir(results_dir);
    fprintf('\nCarpeta creada: %s\n', results_dir);
end

% .mat: mdf_reps [n_reps x 7], slopes [1x7], r2 [1x7] + metadatos
mat_path = fullfile(results_dir, ['mdf2_' base_name '.mat']);
save(mat_path, 'mdf_reps', 't_rep', 'slopes', 'r2', ...
     'muscle_names', 'n_reps', 't_fat_start', 't_fat_end', ...
     'onset_idx', 'offset_idx', 'Fs', 'fname', 'pct_thr', 'thr');
fprintf('Mat guardado: %s\n', mat_path);

% PNG: figura de regresión lineal (300 dpi)
png_path = fullfile(results_dir, ['mdf_' base_name '_regression.png']);
fig_mdf.Position = [100 100 1000 900];   % tamaño fijo para consistencia entre sujetos
exportgraphics(fig_mdf, png_path, 'Resolution', 300);
fprintf('PNG guardado: %s\n', png_path);

