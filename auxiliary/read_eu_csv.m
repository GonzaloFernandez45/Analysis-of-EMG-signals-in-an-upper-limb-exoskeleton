function raw = read_eu_csv(file_path)
% READ_EU_CSV  Lee CSVs de Delsys en formato europeo. Compatible R2016b+.
%
% Detecta automáticamente dos formatos:
%
%   FORMATO A  (separador=';', decimal=',')
%     Cabecera: "X [s];Avanti sensor 1: EMG 1 [V];..."
%     8 columnas: [tiempo, ch1..ch7]
%
%   FORMATO B  (separador=',' y decimal=',' — exportación antigua)
%     Cabecera opcional de metadata (Label:...) + fila de columnas
%     14 columnas: [t,ch1, t,ch2, ..., t,ch7]
%     → salida normalizada a 8 columnas [tiempo, ch1..ch7]

    % ── Leer archivo ──────────────────────────────────────────────────────
    fid = fopen(file_path, 'r');
    if fid == -1, error('No se puede abrir: %s', file_path); end
    content = fread(fid, '*char')';
    fclose(fid);

    % ── Separar en líneas ─────────────────────────────────────────────────
    lines = regexp(content, '\r?\n', 'split');

    % ── Localizar fila de cabecera de columnas ────────────────────────────
    hdr_line = 0;
    for k = 1:numel(lines)
        ln = strtrim(lines{k});
        if length(ln) >= 3 && (strncmpi(ln, 'X[s]', 4) || strncmpi(ln, 'X [s]', 5))
            hdr_line = k;
            break;
        end
    end
    if hdr_line == 0
        error('Cabecera X[s] no encontrada en: %s', file_path);
    end

    header = lines{hdr_line};

    % ── Detectar formato por el separador de la cabecera ─────────────────
    if ~isempty(strfind(header, ';'))
        fmt = 'A';
    else
        fmt = 'B';
    end

    % ── Extraer líneas de datos (sin vacías) ──────────────────────────────
    data_lines = lines(hdr_line+1 : end);
    data_lines = data_lines(~cellfun(@(x) isempty(strtrim(x)), data_lines));

    % ── Parsear según formato ─────────────────────────────────────────────
    if strcmp(fmt, 'A')
        % ── FORMATO A: simple ─────────────────────────────────────────────
        ncols       = numel(regexp(header, ';', 'split'));
        data_str    = strjoin(data_lines, ' ');
        data_str    = strrep(data_str, ',', '.');
        data_str    = strrep(data_str, ';', ' ');
        C           = textscan(data_str, '%f', 'Delimiter', ' ', ...
                               'MultipleDelimsAsOne', true, 'EmptyValue', NaN);
        data        = C{1};
        nrows       = floor(numel(data) / ncols);
        raw         = reshape(data(1:nrows*ncols), ncols, nrows)';

    else
        % ── FORMATO B: coma como separador Y decimal ──────────────────────
        % Estrategia: usar regexp para extraer tokens numéricos completos
        % incluyendo los que usan coma como decimal (ej. "23,27493" o "-1,493E-05").
        % Patrón: número con parte decimal opcional (punto O coma) y exponente opcional.
        num_pat = '[-+]?\d+(?:[.,]\d+)?(?:[Ee][-+]?\d+)?';

        data_str = strjoin(data_lines, ' ');
        tokens   = regexp(data_str, num_pat, 'match');

        % Convertir cada token: reemplazar coma decimal por punto
        data = zeros(numel(tokens), 1);
        for i = 1:numel(tokens)
            data(i) = str2double(strrep(tokens{i}, ',', '.'));
        end

        % Número de columnas desde la cabecera
        ncols = numel(regexp(header, ',', 'split'));
        nrows = floor(numel(data) / ncols);
        if nrows == 0
            error('No se pudieron parsear datos de %s', file_path);
        end
        raw_full = reshape(data(1:nrows*ncols), ncols, nrows)';

        % Estructura B: [t ch1 t ch2 ... t ch7] — 14 cols
        % Queremos:     [t ch1 ch2 ... ch7]      —  8 cols
        % col 1 = tiempo, cols pares = datos (2,4,6,8,10,12,14)
        raw = raw_full(:, [1, 2:2:ncols]);
    end

end
