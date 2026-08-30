% =========================================================================
% step3_statistics.m
% =========================================================================
% Paired statistical analysis: EXO vs NOEXO  (N subjects)
%
% Folder structure expected:
%   [root]/[subject]/results/pre_*.mat         ← from step1
%   [root]/[subject]/results/results/mdf5_*.mat ← from step2
%
% THREE ANALYSIS BLOCKS (Chapter 5):
%   BLOCK 1 — Baseline RMS %MVC       EXO vs NOEXO
%   BLOCK 2 — Post-fatigue RMS %MVC   EXO vs NOEXO
%   BLOCK 3 — Delta MDF + endurance   EXO vs NOEXO
%
% Normality: Shapiro-Wilk (Royston 1992) on paired differences
%   Normal     → paired t-test  + Cohen's d + 95% CI
%   Non-normal → Wilcoxon signed-rank + rank-biserial r
%
% Multiple comparisons: Bonferroni alpha = 0.05/7 for all 7-channel tests
%   n_reps and block_duration_s: uncorrected alpha = 0.05
%
% Output: console tables + step3_results.mat (saved in [root])
%
% Compatibility: MATLAB R2016b+  |  Requires Statistics Toolbox
% =========================================================================

clear; clc;

ALPHA      = 0.05;
N_CH       = 7;
ALPHA_BONF = ALPHA / N_CH;   % 0.00714

muscle_names = {'AD','LD','PD','UT','BB','TB','ECR'};

% -------------------------------------------------------------------------
% 1. Select root folder (e.g. TFG/gente/)
% -------------------------------------------------------------------------
root = uigetdir('C:\Users\gzomo\TFG\gente', ...
    'Select root subjects folder (e.g. gente/)');
if isequal(root, 0), error('No folder selected.'); end

% Get subject subdirectories
entries   = dir(root);
is_subdir = [entries.isdir] & ~ismember({entries.name}, {'.', '..'});
subj_dirs = entries(is_subdir);
fprintf('Subject folders found: %d\n', numel(subj_dirs));

% -------------------------------------------------------------------------
% 2. Collect all files, parse and store in database
% -------------------------------------------------------------------------
%  db.(subjID).pre_EXO   = struct from pre_*_EXO.mat
%  db.(subjID).pre_NOEXO = struct from pre_*_NOEXO.mat
%  db.(subjID).mdf5_EXO  = struct from mdf5_*_EXO.mat
%  db.(subjID).mdf5_NOEXO= struct from mdf5_*_NOEXO.mat

db = struct();

for di = 1:numel(subj_dirs)
    subj_folder = fullfile(root, subj_dirs(di).name);

    % --- pre_*.mat  (in [subj]/results/) ---
    pre_d = dir(fullfile(subj_folder, 'results', 'pre_*.mat'));
    for fi = 1:numel(pre_d)
        [sid, cond] = parse_subject(pre_d(fi).name, 'pre_');
        if isempty(sid) || isempty(cond), continue; end
        S = load(fullfile(pre_d(fi).folder, pre_d(fi).name));
        % Verify condition from file if ambiguous
        if isfield(S, 'condition') && isempty(cond)
            cond = upper(strtrim(S.condition));
        end
        db.(sid).(['pre_' cond]) = S;
    end

    % --- mdf5_*.mat  (in [subj]/results/results/) ---
    mdf5_d = dir(fullfile(subj_folder, 'results', 'results', 'mdf5_*.mat'));
    for fi = 1:numel(mdf5_d)
        [sid, cond] = parse_subject(mdf5_d(fi).name, 'mdf5_');
        if isempty(sid) || isempty(cond), continue; end
        S = load(fullfile(mdf5_d(fi).folder, mdf5_d(fi).name));
        if isfield(S, 'condition') && isempty(cond)
            cond = upper(strtrim(S.condition));
        end
        db.(sid).(['mdf5_' cond]) = S;
    end
end

subj_list = fieldnames(db);
n_subj    = numel(subj_list);
fprintf('Subjects with data: %d\n\n', n_subj);

% -------------------------------------------------------------------------
% 3. Build paired matrices  [n_subj × 7] or [n_subj × 1]
% -------------------------------------------------------------------------
rms_base_exo   = nan(n_subj, N_CH);
rms_base_noexo = nan(n_subj, N_CH);
rms_post_exo   = nan(n_subj, N_CH);
rms_post_noexo = nan(n_subj, N_CH);
dmdf_exo       = nan(n_subj, N_CH);
dmdf_noexo     = nan(n_subj, N_CH);
nreps_exo      = nan(n_subj, 1);
nreps_noexo    = nan(n_subj, 1);
dur_exo        = nan(n_subj, 1);
dur_noexo      = nan(n_subj, 1);

for si = 1:n_subj
    id = subj_list{si};
    d  = db.(id);

    if isfield(d, 'pre_EXO')
        rms_base_exo(si,:) = d.pre_EXO.rms_base_mvc;
        rms_post_exo(si,:) = d.pre_EXO.rms_post_mvc;
    end
    if isfield(d, 'pre_NOEXO')
        rms_base_noexo(si,:) = d.pre_NOEXO.rms_base_mvc;
        rms_post_noexo(si,:) = d.pre_NOEXO.rms_post_mvc;
    end
    if isfield(d, 'mdf5_EXO')
        dmdf_exo(si,:) = d.mdf5_EXO.delta_mdf_mean;
        nreps_exo(si)  = d.mdf5_EXO.n_reps;
        if isfield(d.mdf5_EXO, 'block_duration_s')
            dur_exo(si) = d.mdf5_EXO.block_duration_s;
        end
    end
    if isfield(d, 'mdf5_NOEXO')
        dmdf_noexo(si,:) = d.mdf5_NOEXO.delta_mdf_mean;
        nreps_noexo(si)  = d.mdf5_NOEXO.n_reps;
        if isfield(d.mdf5_NOEXO, 'block_duration_s')
            dur_noexo(si) = d.mdf5_NOEXO.block_duration_s;
        end
    end
end

% Report missing pairs
fprintf('Pairing summary:\n');
for si = 1:n_subj
    id = subj_list{si};
    has = fieldnames(db.(id));
    fprintf('  %-20s  %s\n', id, strjoin(has, '  '));
end
fprintf('\n');

% =========================================================================
% BLOCK 1 — Baseline RMS %MVC: EXO vs NOEXO
% =========================================================================
fprintf('=================================================================\n');
fprintf('BLOCK 1: Baseline RMS %%MVC  EXO vs NOEXO\n');
fprintf('         Bonferroni alpha = %.5f\n', ALPHA_BONF);
fprintf('=================================================================\n');
res_base = run_block_7ch(rms_base_exo, rms_base_noexo, ...
                          muscle_names, ALPHA_BONF, ALPHA);

% =========================================================================
% BLOCK 2 — Post-fatigue RMS %MVC: EXO vs NOEXO
% =========================================================================
fprintf('\n=================================================================\n');
fprintf('BLOCK 2: Post-fatigue RMS %%MVC  EXO vs NOEXO\n');
fprintf('         Bonferroni alpha = %.5f\n', ALPHA_BONF);
fprintf('=================================================================\n');
res_post = run_block_7ch(rms_post_exo, rms_post_noexo, ...
                          muscle_names, ALPHA_BONF, ALPHA);

% =========================================================================
% BLOCK 3A — Delta MDF: EXO vs NOEXO
% =========================================================================
fprintf('\n=================================================================\n');
fprintf('BLOCK 3A: Delta MDF (Hz)  EXO vs NOEXO\n');
fprintf('          Bonferroni alpha = %.5f\n', ALPHA_BONF);
fprintf('=================================================================\n');
res_mdf = run_block_7ch(dmdf_exo, dmdf_noexo, ...
                         muscle_names, ALPHA_BONF, ALPHA);

% =========================================================================
% BLOCK 3B — Endurance outcomes
% =========================================================================
fprintf('\n=================================================================\n');
fprintf('BLOCK 3B: Endurance  EXO vs NOEXO  (alpha = %.2f)\n', ALPHA);
fprintf('=================================================================\n');
res_end = struct();
res_end.nreps    = run_scalar_test(nreps_exo, nreps_noexo, ...
                                    'n_reps', ALPHA);
res_end.duration = run_scalar_test(dur_exo, dur_noexo, ...
                                    'block_duration_s', ALPHA);

% =========================================================================
% Save results
% =========================================================================
save_path = fullfile(root, 'step3_results.mat');
save(save_path, ...
    'res_base', 'res_post', 'res_mdf', 'res_end', ...
    'muscle_names', 'ALPHA', 'ALPHA_BONF', 'N_CH', ...
    'rms_base_exo', 'rms_base_noexo', ...
    'rms_post_exo', 'rms_post_noexo', ...
    'dmdf_exo', 'dmdf_noexo', ...
    'nreps_exo', 'nreps_noexo', ...
    'dur_exo', 'dur_noexo', ...
    'subj_list');

fprintf('\nResults saved: %s\n', save_path);
fprintf('Done.\n');

% =========================================================================
% LOCAL FUNCTIONS  (must be at end of script — MATLAB R2016b+)
% =========================================================================

function [sid, cond] = parse_subject(fname, prefix)
% Extract subject ID and condition from filename.
% e.g. 'pre_p1_EXO.mat' with prefix 'pre_' → sid='p1', cond='EXO'
    base = strrep(fname, '.mat', '');
    base = strrep(base, prefix, '');
    if ~isempty(regexpi(base, '_NOEXO$'))
        cond = 'NOEXO';
        sid  = regexprep(base, '_NOEXO$', '', 'ignorecase');
    elseif ~isempty(regexpi(base, '_EXO$'))
        cond = 'EXO';
        sid  = regexprep(base, '_EXO$', '', 'ignorecase');
    else
        cond = '';
        sid  = base;
    end
    sid = matlab.lang.makeValidName(sid);
end


function res = run_block_7ch(mat_exo, mat_noexo, mnames, alpha_bonf, alpha_nom)
% Paired tests for a [N×7] pair of matrices (one column per muscle).

    N_CH = size(mat_exo, 2);

    res.n           = nan(1, N_CH);
    res.mean_exo    = nan(1, N_CH);
    res.mean_noexo  = nan(1, N_CH);
    res.mean_diff   = nan(1, N_CH);
    res.sd_diff     = nan(1, N_CH);
    res.sw_p        = nan(1, N_CH);
    res.normal_diff = false(1, N_CH);
    res.p_test      = nan(1, N_CH);
    res.test_used   = cell(1, N_CH);
    res.effect_size = nan(1, N_CH);
    res.effect_name = cell(1, N_CH);
    res.ci_lo       = nan(1, N_CH);
    res.ci_hi       = nan(1, N_CH);
    res.sig_bonf    = false(1, N_CH);
    res.sig_nom     = false(1, N_CH);

    fprintf('%-5s  %7s  %7s  %7s  %6s  %6s  %8s  %-8s  %-5s %5s  %s\n', ...
        'Musc', 'M_EXO', 'M_NOEXO', 'Mdiff', 'SDdiff', 'SW_p', ...
        'p_test', 'Test', 'ES_name', 'ES', 'Sig');
    fprintf('%s\n', repmat('-', 1, 88));

    for ch = 1:N_CH
        e  = mat_exo(:, ch);
        no = mat_noexo(:, ch);
        ok = ~isnan(e) & ~isnan(no);
        e_ok  = e(ok);
        no_ok = no(ok);
        d_ok  = e_ok - no_ok;
        n     = sum(ok);

        if n < 5
            fprintf('%-5s  insufficient data (n=%d)\n', mnames{ch}, n);
            continue;
        end

        res.n(ch)          = n;
        res.mean_exo(ch)   = mean(e_ok);
        res.mean_noexo(ch) = mean(no_ok);
        res.mean_diff(ch)  = mean(d_ok);
        res.sd_diff(ch)    = std(d_ok);

        % Shapiro-Wilk on paired differences
        [sw_h, sw_p] = sw_test(d_ok, alpha_nom);
        res.sw_p(ch)        = sw_p;
        res.normal_diff(ch) = ~sw_h;

        if ~sw_h
            % Normal → paired t-test + Cohen's d
            [~, p_t, ci_t]      = ttest(e_ok, no_ok, 'Alpha', alpha_nom);
            res.p_test(ch)      = p_t;
            res.test_used{ch}   = 't-test';
            res.effect_size(ch) = mean(d_ok) / std(d_ok);   % Cohen's d
            res.effect_name{ch} = 'd';
            res.ci_lo(ch)       = ci_t(1);
            res.ci_hi(ch)       = ci_t(2);
        else
            % Non-normal → Wilcoxon + rank-biserial r
            [p_w, ~, stats_w]   = signrank(e_ok, no_ok);
            T_plus              = stats_w.signedrank;
            r_rb                = (2*T_plus / (n*(n+1)/2)) - 1;
            res.p_test(ch)      = p_w;
            res.test_used{ch}   = 'Wilcoxon';
            res.effect_size(ch) = r_rb;
            res.effect_name{ch} = 'r_rb';
            res.ci_lo(ch)       = NaN;
            res.ci_hi(ch)       = NaN;
        end

        res.sig_bonf(ch) = res.p_test(ch) < alpha_bonf;
        res.sig_nom(ch)  = res.p_test(ch) < alpha_nom;

        if res.sig_bonf(ch),     sig_str = '*** BONF';
        elseif res.sig_nom(ch),  sig_str = '*   nom';
        else,                    sig_str = '';
        end

        norm_flag = '';
        if sw_p < alpha_nom, norm_flag = '!'; end

        fprintf('%-5s  %7.2f  %7.2f  %7.2f  %6.2f  %5.3f%s  %8.4f  %-8s  %-5s %5.3f  %s\n', ...
            mnames{ch}, res.mean_exo(ch), res.mean_noexo(ch), ...
            res.mean_diff(ch), res.sd_diff(ch), sw_p, norm_flag, ...
            res.p_test(ch), res.test_used{ch}, res.effect_name{ch}, ...
            res.effect_size(ch), sig_str);
    end

    fprintf('\n  *** BONF = p < %.5f (Bonferroni-corrected)\n', alpha_bonf);
    fprintf('  *   nom  = p < %.2f (nominal, not corrected)\n', alpha_nom);
    fprintf('  !        = Shapiro-Wilk rejected normality → Wilcoxon used\n');
end


function res = run_scalar_test(vec_exo, vec_noexo, label, alpha)
% Paired test for a single scalar outcome.

    ok    = ~isnan(vec_exo) & ~isnan(vec_noexo);
    e_ok  = vec_exo(ok);
    no_ok = vec_noexo(ok);
    d_ok  = e_ok - no_ok;
    n     = sum(ok);

    res = struct();
    if n < 5
        fprintf('  %-24s  insufficient data (n=%d)\n', label, n);
        return;
    end

    res.n          = n;
    res.mean_exo   = mean(e_ok);
    res.mean_noexo = mean(no_ok);
    res.mean_diff  = mean(d_ok);
    res.sd_diff    = std(d_ok);

    [sw_h, sw_p] = sw_test(d_ok, alpha);
    res.sw_p   = sw_p;
    res.normal = ~sw_h;

    if ~sw_h
        [~, p_t, ci_t]  = ttest(e_ok, no_ok, 'Alpha', alpha);
        res.p           = p_t;
        res.test        = 't-test';
        res.es          = mean(d_ok) / std(d_ok);
        res.es_name     = 'd';
        res.ci          = ci_t;
    else
        [p_w, ~, stats_w] = signrank(e_ok, no_ok);
        T_plus  = stats_w.signedrank;
        r_rb    = (2*T_plus / (n*(n+1)/2)) - 1;
        res.p       = p_w;
        res.test    = 'Wilcoxon';
        res.es      = r_rb;
        res.es_name = 'r_rb';
        res.ci      = [NaN NaN];
    end

    res.sig = res.p < alpha;
    sig_str = ''; if res.sig, sig_str = ' *'; end
    norm_str = ''; if sw_p < alpha, norm_str = ' [SW!]'; end

    if strcmp(res.test, 't-test')
        fprintf('  %-24s  EXO=%7.2f  NOEXO=%7.2f  diff=%6.2f±%5.2f  SW_p=%.3f%s  p=%.4f (%s)  %s=%.3f  CI=[%.2f,%.2f]%s\n', ...
            label, res.mean_exo, res.mean_noexo, res.mean_diff, res.sd_diff, ...
            sw_p, norm_str, res.p, res.test, res.es_name, res.es, ...
            res.ci(1), res.ci(2), sig_str);
    else
        fprintf('  %-24s  EXO=%7.2f  NOEXO=%7.2f  diff=%6.2f±%5.2f  SW_p=%.3f%s  p=%.4f (%s)  %s=%.3f%s\n', ...
            label, res.mean_exo, res.mean_noexo, res.mean_diff, res.sd_diff, ...
            sw_p, norm_str, res.p, res.test, res.es_name, res.es, sig_str);
    end
end


function [H, p_val] = sw_test(x, alpha)
% -------------------------------------------------------------------------
% SW_TEST  Wrapper for Shapiro-Wilk normality test.
%
% Calls swtest() (Ben Saïda 2014, MATLAB File Exchange ID 13964) if
% available — place swtest.m in the same folder as this script.
%
% Falls back to lillietest() (Lilliefors) if swtest.m is not found,
% with a one-time console warning.
%
% [H, p] = sw_test(x, alpha)
%   x     : data vector (NaN values ignored)
%   alpha : significance level
%   H     : 1 = reject normality; 0 = cannot reject
%   p_val : p-value
%
% To obtain swtest.m:
%   https://www.mathworks.com/matlabcentral/fileexchange/13964
% -------------------------------------------------------------------------

    x     = x(:);
    x     = x(~isnan(x));
    H     = false;
    p_val = 1;

    if numel(x) < 3, return; end

    if ~exist('swtest', 'file')
        error(['swtest.m not found. Download it from MATLAB File Exchange ID 13964:\n' ...
               'https://www.mathworks.com/matlabcentral/fileexchange/13964\n' ...
               'Place swtest.m in the same folder as step3_statistics.m and re-run.']);
    end
    [H, p_val] = swtest(x, alpha);
end
