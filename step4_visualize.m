%% step4_visualize.m
% Publication-ready figures for the TFG.
%
% Figures produced:
%   1. RMS (%MVC) — grouped bars + individual subject dots (baseline & post-fat)
%   2. MDF slope  — grouped bars + individual subject dots
%   3. Radar (spider) chart — mean EXO vs NOEXO activation per muscle
%   4. CSV summary table — all statistical results for copy-paste into Word
%
% Requires: statistics_results.mat  (from step3_statistics.m)
%           rms_pX.mat + mvc_reference.mat  (same files used in step3)
%           mdf_pX.mat

clear; clc; close all;

%% ── 0. STYLE ─────────────────────────────────────────────────────────────

col_exo   = [0.18 0.38 0.75];   % blue
col_noexo = [0.78 0.18 0.18];   % red
col_dot   = [0.15 0.15 0.15];   % dark gray for individual dots

set(groot, 'defaultAxesFontName',  'Helvetica', ...
           'defaultAxesFontSize',   10, ...
           'defaultAxesLineWidth',  0.8, ...
           'defaultAxesTickDir',    'out', ...
           'defaultAxesBox',        'off', ...
           'defaultFigureColor',    'w');

results_dir = 'C:\Users\gzomo\TFG\results';
mvc_dir     = 'C:\Users\gzomo\TFG\mvc';
output_dir  = results_dir;

muscle_names  = {'AD'; 'LD'; 'PD'; 'UT'; 'BB'; 'TB'; 'ECR'};
muscle_long   = {'Anterior Deltoid'; 'Lateral Deltoid'; 'Posterior Deltoid'; ...
                 'Upper Trapezius';  'Biceps Brachii';  'Triceps Brachii'; 'ECR'};
n_muscles = 7;

%% ── 1. LOAD DATA ─────────────────────────────────────────────────────────

% Statistics results
stat = load(fullfile(results_dir, 'statistics_results.mat'));
res  = stat.results;

% RMS + MVC files
[rms_files, rms_path] = uigetfile('rms_*.mat', ...
    'Select RMS .mat files (all participants)', results_dir, 'MultiSelect', 'on');
if ischar(rms_files), rms_files = {rms_files}; end
n_subjects = numel(rms_files);

[mvc_files, mvc_path] = uigetfile('*mvc_reference.mat', ...
    'Select MVC reference .mat files (same order)', mvc_dir, 'MultiSelect', 'on');
use_mvc = ~isequal(mvc_files, 0);
if use_mvc && ischar(mvc_files), mvc_files = {mvc_files}; end

[mdf_files, mdf_path] = uigetfile('mdf_*.mat', ...
    'Select MDF .mat files (all participants)', results_dir, 'MultiSelect', 'on');
if ischar(mdf_files), mdf_files = {mdf_files}; end

% Rebuild matrices
rms_base_exo   = zeros(n_subjects, n_muscles);
rms_base_noexo = zeros(n_subjects, n_muscles);
rms_post_exo   = zeros(n_subjects, n_muscles);
rms_post_noexo = zeros(n_subjects, n_muscles);
mdf_slope_exo   = zeros(n_subjects, n_muscles);
mdf_slope_noexo = zeros(n_subjects, n_muscles);

for s = 1:n_subjects
    r  = load(fullfile(rms_path, rms_files{s}));
    rr = r.rms_results;
    eb = rr.EXO.baseline;  ep = rr.EXO.postfat;
    nb = rr.NOEXO.baseline; np = rr.NOEXO.postfat;

    if use_mvc
        mvc_raw  = load(fullfile(mvc_path, mvc_files{s}));
        f_mvc    = fieldnames(mvc_raw);
        mvc      = mvc_raw.(f_mvc{1});
        mvc_norm = zeros(n_muscles,1);
        if isfield(mvc,'channel_idx') && isfield(mvc,'RMS_MVC')
            for m = 1:numel(mvc.channel_idx)
                ch = mvc.channel_idx(m);
                if ch >= 1 && ch <= n_muscles
                    mvc_norm(ch) = mvc.RMS_MVC(m) * 1e6;
                end
            end
        end
        mvc_norm(mvc_norm==0) = NaN;
        eb = eb ./ mvc_norm * 100;  ep = ep ./ mvc_norm * 100;
        nb = nb ./ mvc_norm * 100;  np = np ./ mvc_norm * 100;
    end

    rms_base_exo(s,:)   = eb';
    rms_base_noexo(s,:) = nb';
    rms_post_exo(s,:)   = ep';
    rms_post_noexo(s,:) = np';

    md = load(fullfile(mdf_path, mdf_files{s}));
    md = md.mdf_results;
    mdf_slope_exo(s,:)   = md.EXO.slope';
    mdf_slope_noexo(s,:) = md.NOEXO.slope';
end

rms_unit = iif(use_mvc, '%MVC', 'µV');

%% ── FIGURE 1 — RMS with individual dots ─────────────────────────────────

fig1 = figure('Name', 'Fig1 RMS individual', 'NumberTitle', 'off', ...
              'Position', [40 40 1400 580]);

titles_bl = {'Baseline', 'Post-fatigue'};
pairs = {rms_base_exo, rms_base_noexo; rms_post_exo, rms_post_noexo};

bw  = 0.28;
x   = 1:n_muscles;
jit = 0.06;   % dot jitter

for sub = 1:2
    ax = subplot(1,2,sub);
    hold(ax,'on');

    Xe = pairs{sub,1};
    Xn = pairs{sub,2};
    me = mean(Xe,1,'omitnan');
    mn = mean(Xn,1,'omitnan');
    se_e = std(Xe,0,1,'omitnan') / sqrt(n_subjects);
    se_n = std(Xn,0,1,'omitnan') / sqrt(n_subjects);

    % Bars
    b1 = bar(ax, x-bw/2, me, bw, 'FaceColor', col_exo,   'FaceAlpha', 0.75, 'EdgeColor','none');
    b2 = bar(ax, x+bw/2, mn, bw, 'FaceColor', col_noexo, 'FaceAlpha', 0.75, 'EdgeColor','none');

    % SEM error bars
    errorbar(ax, x-bw/2, me, se_e, 'k.', 'LineWidth',1.2, 'HandleVisibility','off');
    errorbar(ax, x+bw/2, mn, se_n, 'k.', 'LineWidth',1.2, 'HandleVisibility','off');

    % Individual dots + lines connecting EXO-NOEXO per subject
    for s = 1:n_subjects
        jx = (rand(1,n_muscles)-0.5)*jit;
        % connecting line
        for ch = 1:n_muscles
            plot(ax, [ch-bw/2+jx(ch), ch+bw/2+jx(ch)], ...
                 [Xe(s,ch), Xn(s,ch)], '-', 'Color', [0.5 0.5 0.5 0.4], ...
                 'LineWidth', 0.8, 'HandleVisibility','off');
        end
        plot(ax, x-bw/2+jx, Xe(s,:), 'o', 'MarkerFaceColor', col_dot, ...
             'MarkerEdgeColor','w', 'MarkerSize',5, 'HandleVisibility','off');
        plot(ax, x+bw/2+jx, Xn(s,:), 'o', 'MarkerFaceColor', col_dot, ...
             'MarkerEdgeColor','w', 'MarkerSize',5, 'HandleVisibility','off');
    end

    set(ax, 'XTick', x, 'XTickLabel', muscle_long, 'XTickLabelRotation', 35);
    ylabel(ax, sprintf('RMS (%s)', rms_unit));
    title(ax, titles_bl{sub}, 'FontSize',12, 'FontWeight','bold');
    if sub==1
        legend(ax, [b1 b2], {'EXO','NOEXO'}, 'Location','northwest', 'Box','off');
    end
    grid(ax,'on'); ax.GridAlpha = 0.15;
end

sgtitle(fig1, sprintf('Muscle Activation — EXO vs NOEXO (%s, mean ± SEM + individual)', rms_unit), ...
        'FontSize',13, 'FontWeight','bold');

exportgraphics(fig1, fullfile(output_dir,'fig1_RMS_individual.png'), 'Resolution',300);
fprintf('Fig 1 saved.\n');

%% ── FIGURE 2 — MDF slope with individual dots ────────────────────────────

fig2 = figure('Name','Fig2 MDF slope individual','NumberTitle','off', ...
              'Position',[60 60 1050 500]);
hold on;

me_e  = mean(mdf_slope_exo,   1,'omitnan');
me_n  = mean(mdf_slope_noexo, 1,'omitnan');
se_e  = std(mdf_slope_exo,   0,1,'omitnan') / sqrt(n_subjects);
se_n  = std(mdf_slope_noexo, 0,1,'omitnan') / sqrt(n_subjects);

b1 = bar(x-bw/2, me_e, bw, 'FaceColor',col_exo,   'FaceAlpha',0.75,'EdgeColor','none');
b2 = bar(x+bw/2, me_n, bw, 'FaceColor',col_noexo, 'FaceAlpha',0.75,'EdgeColor','none');
errorbar(x-bw/2, me_e, se_e, 'k.','LineWidth',1.2,'HandleVisibility','off');
errorbar(x+bw/2, me_n, se_n, 'k.','LineWidth',1.2,'HandleVisibility','off');
yline(0,'k--','LineWidth',1,'HandleVisibility','off');

for s = 1:n_subjects
    jx = (rand(1,n_muscles)-0.5)*jit;
    for ch = 1:n_muscles
        plot([ch-bw/2+jx(ch), ch+bw/2+jx(ch)], ...
             [mdf_slope_exo(s,ch), mdf_slope_noexo(s,ch)], ...
             '-','Color',[0.5 0.5 0.5 0.4],'LineWidth',0.8,'HandleVisibility','off');
    end
    plot(x-bw/2+jx, mdf_slope_exo(s,:),   'o','MarkerFaceColor',col_dot, ...
         'MarkerEdgeColor','w','MarkerSize',5,'HandleVisibility','off');
    plot(x+bw/2+jx, mdf_slope_noexo(s,:), 'o','MarkerFaceColor',col_dot, ...
         'MarkerEdgeColor','w','MarkerSize',5,'HandleVisibility','off');
end

set(gca,'XTick',x,'XTickLabel',muscle_long,'XTickLabelRotation',35);
ylabel('MDF regression slope (Hz/s)');
title('MDF Fatigue Rate — EXO vs NOEXO (mean ± SEM + individual)','FontSize',12,'FontWeight','bold');
legend([b1 b2],{'EXO','NOEXO'},'Location','southeast','Box','off');
grid on; ax2 = gca; ax2.GridAlpha = 0.15;

exportgraphics(fig2, fullfile(output_dir,'fig2_MDF_slope_individual.png'),'Resolution',300);
fprintf('Fig 2 saved.\n');

%% ── FIGURE 3 — Radar chart ───────────────────────────────────────────────
% Mean RMS post-fatigue EXO vs NOEXO, normalised to max across both
% conditions so all muscles fit on [0,1].

fig3 = figure('Name','Fig3 Radar','NumberTitle','off','Position',[80 80 700 650]);
ax3  = axes('Position',[0.05 0.05 0.9 0.85]);
hold(ax3,'on'); axis(ax3,'off'); axis(ax3,'equal');

N    = n_muscles;
angles = linspace(pi/2, pi/2 + 2*pi, N+1);   % start at top
angles = angles(1:end);

me_exo   = mean(rms_post_exo,   1,'omitnan');
me_noexo = mean(rms_post_noexo, 1,'omitnan');
ref      = max([me_exo; me_noexo], [], 1);
ref(ref==0) = 1;
ve = me_exo   ./ ref;
vn = me_noexo ./ ref;

% Grid rings
n_rings = 5;
for r = 1:n_rings
    rr = r/n_rings;
    th = linspace(0, 2*pi, 200);
    plot(ax3, rr*cos(th), rr*sin(th), '-', 'Color',[0.85 0.85 0.85],'LineWidth',0.5);
    text(ax3, 0, rr+0.03, sprintf('%.0f%%', rr*100), ...
         'HorizontalAlignment','center','FontSize',7,'Color',[0.6 0.6 0.6]);
end

% Spokes
for k = 1:N
    plot(ax3, [0 cos(angles(k))], [0 sin(angles(k))], '-','Color',[0.8 0.8 0.8],'LineWidth',0.5);
    lbl_r = 1.15;
    text(ax3, lbl_r*cos(angles(k)), lbl_r*sin(angles(k)), muscle_names{k}, ...
         'HorizontalAlignment','center','FontSize',9,'FontWeight','bold');
end

% Polygons
xe = [ve.*cos(angles(1:N)), ve(1)*cos(angles(1))];
ye = [ve.*sin(angles(1:N)), ve(1)*sin(angles(1))];
xn = [vn.*cos(angles(1:N)), vn(1)*cos(angles(1))];
yn = [vn.*sin(angles(1:N)), vn(1)*sin(angles(1))];

fill(ax3, xn, yn, col_noexo, 'FaceAlpha',0.20,'EdgeColor',col_noexo,'LineWidth',1.5);
fill(ax3, xe, ye, col_exo,   'FaceAlpha',0.20,'EdgeColor',col_exo,   'LineWidth',1.5);
plot(ax3, xn, yn, '-','Color',col_noexo,'LineWidth',2);
plot(ax3, xe, ye, '-','Color',col_exo,  'LineWidth',2);
plot(ax3, xe(1:N).*1, ye(1:N).*1, 'o','MarkerFaceColor',col_exo,   'MarkerEdgeColor','w','MarkerSize',7);
plot(ax3, xn(1:N).*1, yn(1:N).*1, 's','MarkerFaceColor',col_noexo, 'MarkerEdgeColor','w','MarkerSize',7);

legend(ax3, {'NOEXO','EXO'}, 'Location','southoutside','Orientation','horizontal','Box','off','FontSize',10);
title(ax3, sprintf('Post-fatigue RMS (%s) — EXO vs NOEXO\n(normalised per muscle to max)', rms_unit), ...
      'FontSize',11,'FontWeight','bold');

exportgraphics(fig3, fullfile(output_dir,'fig3_radar.png'),'Resolution',300);
fprintf('Fig 3 saved.\n');

%% ── FIGURE 4 — CSV summary table ─────────────────────────────────────────

analyses_labels = {'RMS Baseline', 'RMS Post-fat', 'MDF Slope', 'MDF Total Shift'};
fields_list     = {'RMS_Baseline', 'RMS_Post_fat', 'MDF_Slope', 'MDF_Total_Shift'};

csv_path = fullfile(output_dir, 'summary_statistics.csv');
fid = fopen(csv_path, 'w');
fprintf(fid, 'Analysis,Muscle,Mean_Diff,p_raw,p_adj_Bonferroni,Cohens_d,Test,Significant\n');

for a = 1:numel(analyses_labels)
    field = matlab.lang.makeValidName(analyses_labels{a});
    if ~isfield(res, field), continue; end
    r = res.(field);
    for ch = 1:n_muscles
        fprintf(fid, '%s,%s,%.4f,%.4f,%.4f,%.4f,%s,%s\n', ...
                analyses_labels{a}, muscle_long{ch}, ...
                r.mean_diff(ch), r.p_raw(ch), r.p_adj(ch), ...
                r.cohens_d(ch), r.test_used{ch}, ...
                iif(r.sig(ch),'YES','no'));
    end
end
fclose(fid);
fprintf('CSV saved → %s\n', csv_path);

fprintf('\nAll figures and CSV saved to: %s\n', output_dir);

%% ══ LOCAL FUNCTIONS ══════════════════════════════════════════════════════

function out = iif(cond, a, b)
    if cond, out = a; else, out = b; end
end
