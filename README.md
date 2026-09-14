# EMG Signal Analysis in an Upper-Limb Exoskeleton

## Overview

This project processes and analyzes surface EMG (sEMG) signals recorded from
7 muscles (UT, AD, LD, PD, BB, TB, ECR) during a reach-to-grasp fatigue
protocol, performed with and without a passive upper-limb exoskeleton
(POWERUP). The MATLAB scripts cover the full pipeline: data quality control,
MVC extraction, filtering and block segmentation, fatigue analysis via
median frequency (MDF), paired group-level statistics (EXO vs NOEXO), and
publication-ready visualization.

Aquí está el README completo actualizado, listo para pegar:

# EMG Signal Analysis in an Upper-Limb Exoskeleton

## Overview

This project processes and analyzes surface EMG (sEMG) signals recorded from
7 muscles (UT, AD, LD, PD, BB, TB, ECR) during a reach-to-grasp fatigue
protocol, performed with and without a passive upper-limb exoskeleton
(POWERUP). The MATLAB scripts cover the full pipeline: data quality control,
MVC extraction, filtering and block segmentation, fatigue analysis via
median frequency (MDF), paired group-level statistics (EXO vs NOEXO), and
publication-ready visualization.

Developed for the Bachelor's Thesis (TFG) *"Adquisición y análisis de
señales biomédicas en un exoesqueleto de miembro superior para aplicaciones
de rehabilitación."*

## Project Structure

### 1. `data_quality.m`

**Purpose:**
- Runs a quality-control pass on raw EMG recordings before they enter the
  main pipeline (signal integrity, saturation, missing channels).

### 2. `step0_mvc.m`

**Purpose:**
- Loads all MVC (Maximum Voluntary Contraction) recordings of a subject.
- Filters each file and computes the peak RMS activation of every channel.
- Searches for the peak activation of each muscle across all MVC files,
  not just its own recording, to obtain the most representative reference.
- Displays a summary table (muscle × file) and lets the user manually
  adjust individual muscles if needed.
- Saves `mvc_reference.mat`, used by `step1_preprocess.m`.

**Key Features:**
- Peak located on a smoothed signal (robust to artefacts), measured on the
  raw sliding RMS to recover the true peak value.
- Grouped bar-chart summary across files and muscles.
- Manual override per muscle if the automatic selection looks wrong.

### 3. `step1_preprocess.m`

**Purpose:**
- Loads the raw CSV for each recording (EXO or NOEXO), auto-detecting the
  Delsys export format.
- Applies standard filtering (notch 50 Hz + bandpass 20-450 Hz).
- Lets the user select the three protocol blocks (baseline, fatigue,
  post-fatigue) via 6 clicks on the signal.
- Computes RMS %MVC for baseline and post-fatigue using `mvc_reference.mat`.
- Saves `pre_[name].mat`, consumed by `step2_mdf_v5.m` and
  `step3_statistics.m`.

**Key Features:**
- Single step fusing loading, filtering, block selection and normalization
  (replaces the older two-step workflow).
- Auto-detects CSV format (European semicolon/comma or standard).

### 4. `step2_mdf_v5.m`

**Purpose:**
- Reads `pre_[name].mat` and computes Median Frequency (MDF) per repetition
  of the fatigue block.
- Fits a linear regression across repetitions to estimate the fatigue slope
  (spectral compression rate).
- Reports delta MDF in three forms (mean, median, relative %).

**Key Features:**
- AR(3) model via Burg's method (`pburg`), NFFT 1024.
- Edge trim 10% each side (middle 80% of each repetition used).
- Adaptive percentile onset/offset detection referenced to the ECR channel.
- Theil-Sen robust slope and Spearman rho as monotonic-trend descriptors.
- Manual correction interface to add missed repetitions, plus a
  verification figure overlaying onset detection across all 7 channels.

### 5. `step3_statistics.m`

**Purpose:**
- Aggregates `pre_*.mat` and `mdf5_*.mat` across all subject folders.
- Runs paired EXO vs NOEXO statistical tests across three blocks: baseline
  RMS %MVC, post-fatigue RMS %MVC, and delta MDF + endurance outcomes.
- Saves `step3_results.mat`, consumed by `step4_visualize.m`.

**Key Features:**
- Shapiro-Wilk normality test (Royston, 1992) on paired differences.
- Paired t-test + Cohen's d (normal) or Wilcoxon signed-rank +
  rank-biserial r (non-normal).
- Bonferroni correction (alpha = 0.05/7) for the 7-channel tests.

### 6. `step4_visualize.m`

**Purpose:**
- Generates the publication-ready figures for the thesis from
  `step3_results.mat`.

**Key Features:**
- Figures: baseline RMS %MVC, post-fatigue RMS %MVC, delta MDF, MDF
  trajectory across the time-normalized fatigue block (0-100%), and
  endurance outcomes, all as boxplots with mean dot.
- Significance markers: `***` = Bonferroni-corrected, dagger = nominal
  (p < 0.05, uncorrected).

### `auxiliary/` — shared utilities and support analyses

- `read_eu_csv.m`: Delsys CSV loader, auto-detects the European (`;`/`,`)
  and standard export formats. Used by the active pipeline.
- `swtest.m`: Shapiro-Wilk test implementation, used by
  `step3_statistics.m`. Third-party code, see Dependencies.
- `compare_exo_noexo.m`: channel-by-channel visual comparison of RMS
  envelopes between EXO and NOEXO recordings of the same subject.
- `step2_mdf_window_comparison.m`: compares three spectral window
  strategies for MDF estimation; documents why the strategy used in
  `step2_mdf_v5.m` was chosen.

### `archive/` — superseded development versions

Earlier iterations of `step0_mvc`, `step1_rms`, and `step2_mdf` (v2 through
v4c), plus a discarded visualization trial. Kept for traceability, not part
of the active pipeline.

## How to Use

1. **Quality control:** run `step0_data_quality.m` on the raw recordings.
2. **Extract MVC reference:** run `step0_mvc_vFinal.m` on the subject's MVC
   recordings.
3. **Preprocess:** run `step1_preprocess.m` on each EXO/NOEXO CSV to filter,
   segment, and compute RMS %MVC.
4. **Fatigue analysis:** run `step2_mdf_v5.m` to compute MDF and fatigue
   slope per repetition.
5. **Group statistics:** run `step3_statistics.m` across all subject
   folders to get paired EXO vs NOEXO results.
6. **Visualize:** run `step4_visualize.m` to generate the final figures.

## Dependencies

- MATLAB R2016b+
- Statistics and Machine Learning Toolbox (required by
  `step3_statistics.m`)
- EMG data recorded with a Delsys system, exported as `.csv`
- `swtest.m` (Shapiro-Wilk test): third-party implementation by Ahmed Ben
  Saida (2014), included in `auxiliary/`

## Author

Gonzalo Fernández Eizaguirre

Supervisor: Rodrigo Rodríguez Merino
