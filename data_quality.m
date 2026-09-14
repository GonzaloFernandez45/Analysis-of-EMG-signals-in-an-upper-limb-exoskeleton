%% STEP 0 — Control de calidad de datos EMG

clear;
clc;
close all;

%% ===== SELECCIONAR PARTICIPANTE =====

data_dir = 'C:\Users\gzomo\TFG\gente';

participant_dir = uigetdir(data_dir, 'Selecciona la carpeta del participante');

if participant_dir == 0
    error('No se selecciono ningun participante.');
end

[~, participant_name] = fileparts(participant_dir);

fprintf('\n========================================\n');
fprintf('CONTROL DE CALIDAD EMG\n');
fprintf('Participante: %s\n', participant_name);
fprintf('========================================\n\n');

%% ===== CLASIFICAR ARCHIVOS =====

csv_files = dir(fullfile(participant_dir, '*.csv'));

file_exo   = '';
file_noexo = '';

muscle_names = {'AD','LD','PD','UT','BB','TB','ECR'};
mvc_files = struct();

for i = 1:numel(csv_files)

    file_name = csv_files(i).name;

    if contains(file_name, '_NOEXO.csv')
        file_noexo = fullfile(participant_dir, file_name);

    elseif contains(file_name, '_EXO.csv')
        file_exo = fullfile(participant_dir, file_name);

    else
        for m = 1:numel(muscle_names)

            muscle = muscle_names{m};

            if contains(file_name, [muscle '_MVC.csv'])
                mvc_files.(muscle) = fullfile(participant_dir, file_name);
            end

        end
    end
end

%% ===== MOSTRAR CLASIFICACION =====

fprintf('\nArchivos clasificados:\n\n');

fprintf('EXO:    %s\n', file_exo);
fprintf('NOEXO:  %s\n\n', file_noexo);

fprintf('MVC:\n');

for m = 1:numel(muscle_names)
    muscle = muscle_names{m};

    if isfield(mvc_files, muscle)
        fprintf('  %-4s -> %s\n', muscle, mvc_files.(muscle));
    else
        fprintf('  %-4s -> NO ENCONTRADO\n', muscle);
    end
end

%% ===== CARGAR EXO Y NOEXO =====

data_exo   = xlsread(file_exo);
data_noexo = xlsread(file_noexo);

fprintf('\nDatos cargados correctamente.\n');;

%% ===== CARGAR MVC =====

mvc_data = struct();

for m = 1:numel(muscle_names)

    muscle = muscle_names{m};

    file_mvc = mvc_files.(muscle);

    data_mvc = xlsread(file_mvc);

    mvc_data.(muscle).time = data_mvc(:, 1);
    mvc_data.(muscle).emg  = data_mvc(:, 2:end);

end

fprintf('MVC cargados correctamente.\n');

%% ===== SEPARAR TIEMPO Y EMG =====

time_exo  = data_exo(:, 1);
emg_exo   = data_exo(:, 2:end);

time_noexo = data_noexo(:, 1);
emg_noexo  = data_noexo(:, 2:end);

%% ===== CONTROL BASICO =====

n_ch_exo   = size(emg_exo, 2);
n_ch_noexo = size(emg_noexo, 2);

duration_exo   = time_exo(end) - time_exo(1);
duration_noexo = time_noexo(end) - time_noexo(1);

Fs_exo   = 1 / median(diff(time_exo));
Fs_noexo = 1 / median(diff(time_noexo));

nan_exo   = sum(isnan(emg_exo(:)));
nan_noexo = sum(isnan(emg_noexo(:)));

inf_exo   = sum(isinf(emg_exo(:)));
inf_noexo = sum(isinf(emg_noexo(:)));

%% ===== CONTROL BASICO MVC =====

fprintf('\n========================================\n');
fprintf('INTEGRIDAD DE LOS MVC\n');
fprintf('========================================\n');

for m = 1:numel(muscle_names)

    muscle = muscle_names{m};

    time_mvc = mvc_data.(muscle).time;
    emg_mvc  = mvc_data.(muscle).emg;

    n_ch_mvc = size(emg_mvc, 2);

    duration_mvc = time_mvc(end) - time_mvc(1);

    Fs_mvc = 1 / median(diff(time_mvc));

    nan_mvc = sum(isnan(emg_mvc(:)));
    inf_mvc = sum(isinf(emg_mvc(:)));

    fprintf('\n%s MVC\n', muscle);
    fprintf('  Duracion:   %.1f s\n', duration_mvc);
    fprintf('  Fs:         %.2f Hz\n', Fs_mvc);
    fprintf('  Canales:    %d\n', n_ch_mvc);
    fprintf('  NaN:        %d\n', nan_mvc);
    fprintf('  Inf:        %d\n', inf_mvc);

end

%% ===== INFORME BASICO =====

fprintf('\n========================================\n');
fprintf('INTEGRIDAD DE LOS DATOS\n');
fprintf('========================================\n');

fprintf('\nEXO\n');
fprintf('  Duracion:   %.1f s\n', duration_exo);
fprintf('  Fs:         %.2f Hz\n', Fs_exo);
fprintf('  Canales:    %d\n', n_ch_exo);
fprintf('  NaN:        %d\n', nan_exo);
fprintf('  Inf:        %d\n', inf_exo);

fprintf('\nNOEXO\n');
fprintf('  Duracion:   %.1f s\n', duration_noexo);
fprintf('  Fs:         %.2f Hz\n', Fs_noexo);
fprintf('  Canales:    %d\n', n_ch_noexo);
fprintf('  NaN:        %d\n', nan_noexo);
fprintf('  Inf:        %d\n', inf_noexo);