%% ========================================================================
%  MAKASSAR THERMAL ANALYSIS - v12 REVISED
%  Thermal Contrast and Rainfall Extremes
%
%  STRUCTURAL REPAIR:
%  One end-to-end script. No embedded standalone scripts or mid-file local
%  function blocks. Main multivariable model is Rain + RH; RH-only is not
%  duplicated as a second main-text Table 5 model. Figures and tables are
%  generated from the same computed results.
%
%  IMPORTANT:
%  This file is structurally reorganized from the supplied combined script.
%  It has not been executed in MATLAB in this environment. Check data files,
%  MATLAB release compatibility, and required toolboxes before publication.
% ========================================================================
%% =========================================================================
%  MKS vs HND: Quality Control, Homogeneity, Thermal Contrast, ETCCDI
%  Stations : MKS (Coastal, Paotere-Makassar), HND (Airport, Hasanuddin)
%  Purpose  : Full revised analysis pipeline with tie-corrected MK and conventional rank-ACF diagnostics.
%             Supersedes uhi_hnd_mks_etccdi.m by adding the quality-control,
%             homogeneity, descriptive and inferential steps that the
%             manuscript's Methods section describes.
%
%  VERSION  : v12 revised - common complete-year ETCCDI sensitivity window.
%  NOTE     : HR is sensitivity-only for the short annual series (n=14-17).
%             Section 14d now enforces identical common complete years at MKS/HND.
%  PREREQUISITES:
%    - MKS.xlsx and HND.xlsx, Sheet 1 = daily data with columns:
%        Date | Tmean | Tmin | Tmax | Rain | Hummidity
%    - Statistics and Machine Learning Toolbox
%
%  ANALYSIS SEQUENCE:
%    Section 1  : Read data (MKS, HND only)
%    Section 2  : Quality control (completeness, range, outliers,
%                 internal consistency, spatial consistency)
%    Section 3  : Lapse-rate correction + inner join
%    Section 4  : Descriptive statistics                        [Methods 2.3]
%    Section 5  : Normality tests (KS + Lilliefors)             [Methods 2.3]
%    Section 6  : Paired inter-station comparison (MKS - HND) [Methods 2.6]
%    Section 7  : Compute DT = T_MKS - T_HND
%    Section 8  : RELATIVE homogeneity test of DT (SNHT+Pettitt)[Methods 2.3]
%    Section 9  : Annual DT trend, FULL RECORD                  [Table S3]
%    Section 10 : DT ~ HND rainfall + RH, lag 0-3               [Section 3.3]
%                 Wet/dry season models are supplementary sensitivity analyses.
%    Section 12 : ETCCDI rainfall extreme indices               [Table 6]
%    Section 13 : Visualization
%    Section 14d-iii: Hamed-Rao sensitivity diagnostic (short-series caution)
%
%  METHODOLOGICAL NOTE ON HOMOGENEITY (Section 8):
%    SNHT and Pettitt are applied to the MONTHLY, DESEASONALIZED DT series,
%    not to raw daily station series. This is deliberate:
%      (a) Applying these tests to raw daily data violates their assumptions
%          (independence, no seasonal cycle) and inflates false positives.
%      (b) Testing the inter-station DIFFERENCE series evaluates relative
%          homogeneity: a break indicates a change in the relationship between
%          the two station records. It does NOT identify which station or
%          observational factor caused the change, and shared regional signals
%          do not necessarily cancel perfectly.
% =========================================================================
clear; close all; clc;

%% =========================================================================
%  CONFIGURATION
% =========================================================================
ELEV_HND    = 14;     % m, BMKG station metadata
ELEV_MKS    = 1.75;   % m, BMKG station metadata
LAPSE_RATE  = 6.5;    % degC per 1000 m

WET_DAY_THRESH = 1.0; % mm, ETCCDI standard
BASE_YEARS     = [];  % [] = derive station-specific percentile thresholds from the full common record
% The full-record p95/p99 thresholds are reused unchanged for the
% homogeneous-period sensitivity analysis (Section 14d).

% QC thresholds
PHYS_LIM = struct( ...
    'Tmean',    [18 38], ...
    'Tmin',     [15 35], ...
    'Tmax',     [22 42], ...
    'Rain',     [0 500], ...
    'Humidity', [20 100]);
IQR_FACTOR = 3.0;     % outlier flagged if beyond Q1/Q3 -/+ IQR_FACTOR*IQR
Z_THRESH   = 4.0;     % ...AND |z| > Z_THRESH
REMOVE_FLATLINE_RAINFALL = false; % delete flagged runs only after metadata review

% Variables subject to STATISTICAL outlier screening.
% RAINFALL IS DELIBERATELY EXCLUDED. Daily rainfall is strongly right-skewed
% (median = 0 mm, skewness > 2.5), so IQR- and z-score-based rules flag the
% heaviest events as "outliers" - precisely the events the ETCCDI indices are
% designed to measure. Applying such screening to rainfall censors the upper
% tail and destroys the extreme-rainfall signal (verified: it cut RX1day from
% ~160 mm to ~90 mm and eliminated all significant trends). Rainfall errors
% are screened by the physical range check (Section 2b) and the flat-line /
% repeated-value check (Section 2c-ii) instead, which is standard practice.
OUTLIER_VARS = {'Tmean','Tmin','Tmax','Humidity'};

% Temperature outlier screening is applied to DESEASONALIZED anomalies rather
% than raw values, so that legitimate seasonal extremes are not flagged.
DESEASON_OUTLIER = true;

SNHT_MC    = 1000;    % Monte Carlo replicates for SNHT p-value
rng(42,'twister');    % reproducible Monte Carlo significance assessment

% MKS observing-system timeline: a new instrument was installed at a new
% site in July 2023 for trial operation. The previous instrument at the
% original site continued to provide the observations analysed here through
% December 2023. The new system officially commenced operation in January
% 2024. Thus, statistical breaks detected in July-October 2023 cannot be
% directly attributed to the new instrument or site relocation.
%
% BREAK_YM is the first month after the selected pre-break segment. It is
% updated from Section 8 using the earliest statistically supported break.
% Complete calendar years 2009-2022 form the primary temperature-trend
% period; full-record estimates are retained as sensitivity results.
BREAK_YM = 202308;    % configured fallback; Section 8 overwrites if a break is detected
BREAK_DETECTED = false;

%% =========================================================================
%  SECTION 1: Read Data
% =========================================================================
fprintf('=== SECTION 1: READ DATA (MKS, HND) ===\n');
numCols = {'Tmean','Tmin','Tmax','Rain','Hummidity'};

mks_d = readtable('MKS.xlsx', 'Sheet', 1, 'VariableNamingRule','preserve');
hnd_d = readtable('HND.xlsx', 'Sheet', 1, 'VariableNamingRule','preserve');

% Standardize Date/Tanggal and humidity headers before accessing the columns.
% This prevents the common "Unrecognized table variable name 'Tanggal'"
% error when the Excel files use different date-header conventions.
mks_d = standardizeDailyHeaders(mks_d, 'MKS');
hnd_d = standardizeDailyHeaders(hnd_d, 'HND');

mks_d = fixNumericCols(mks_d, numCols); hnd_d = fixNumericCols(hnd_d, numCols);
mks_d = renameDaily(mks_d);             hnd_d = renameDaily(hnd_d);
if numel(unique(mks_d.Date)) ~= height(mks_d)
    error('MKS contains duplicate daily dates. Resolve duplicates before joining.');
end
if numel(unique(hnd_d.Date)) ~= height(hnd_d)
    error('HND contains duplicate daily dates. Resolve duplicates before joining.');
end

fprintf('Data loaded. MKS=%d days, HND=%d days\n', height(mks_d), height(hnd_d));

%% =========================================================================
%  SECTION 2: QUALITY CONTROL
% =========================================================================
fprintf('\n=== SECTION 2: QUALITY CONTROL ===\n');

stations = {'MKS','HND'};
staTbl   = {mks_d, hnd_d};
varsQC   = {'Tmean','Tmin','Tmax','Rain','Humidity'};

% --- 2a. Completeness -----------------------------------------------------
fprintf('\n--- 2a. Data completeness ---\n');
fprintf('%-8s %-10s %8s %10s %9s\n','Station','Variable','N','N_missing','Percent');
fprintf('%s\n', repmat('-',1,48));
for s = 1:2
    tbl = staTbl{s};
    for v = 1:numel(varsQC)
        nm = varsQC{v};
        if ~ismember(nm, tbl.Properties.VariableNames); continue; end
        nmiss = sum(isnan(tbl.(nm)));
        pct   = 100*nmiss/height(tbl);
        flag  = ''; if pct > 5, flag = '  <-- >5% WARNING'; end
        fprintf('%-8s %-10s %8d %10d %8.2f%%%s\n', stations{s}, nm, height(tbl), nmiss, pct, flag);
    end
end

% --- 2b. Range check (physical limits) -----------------------------------
fprintf('\n--- 2b. Range check (values outside physical limits set to NaN) ---\n');
fprintf('%-8s %-10s %10s %9s\n','Station','Variable','N_flagged','Percent');
fprintf('%s\n', repmat('-',1,42));
for s = 1:2
    tbl = staTbl{s};
    for v = 1:numel(varsQC)
        nm = varsQC{v};
        if ~ismember(nm, tbl.Properties.VariableNames); continue; end
        lim  = PHYS_LIM.(nm);
        x    = tbl.(nm);
        bad  = (x < lim(1) | x > lim(2)) & ~isnan(x);
        x(bad) = NaN;
        staTbl{s}.(nm) = x;
        fprintf('%-8s %-10s %10d %8.3f%%\n', stations{s}, nm, sum(bad), 100*sum(bad)/sum(~isnan(x)));
    end
end

% --- 2c-i. Outlier detection, TEMPERATURE & HUMIDITY ONLY ----------------
fprintf('\n--- 2c-i. Outliers (IQR x %.1f AND |Z| > %.1f -> NaN) ---\n', IQR_FACTOR, Z_THRESH);
fprintf('Applied to %s only; rainfall excluded by design (see config).\n', strjoin(OUTLIER_VARS,', '));
if DESEASON_OUTLIER
    fprintf('Screening performed on deseasonalized (calendar-month) anomalies.\n');
end
fprintf('%-8s %-10s %8s %9s %10s\n','Station','Variable','N_IQR','N_Zscore','N_removed');
fprintf('%s\n', repmat('-',1,48));
for s = 1:2
    tbl = staTbl{s};
    moS = month(tbl.Date);
    for v = 1:numel(OUTLIER_VARS)
        nm = OUTLIER_VARS{v};
        if ~ismember(nm, tbl.Properties.VariableNames); continue; end
        x  = tbl.(nm);
        ok = ~isnan(x);

        if DESEASON_OUTLIER
            % Remove calendar-month climatology before screening
            xs = x;
            for m = 1:12
                sel = (moS == m) & ok;
                if any(sel); xs(sel) = x(sel) - mean(x(sel)); end
            end
        else
            xs = x;
        end

        xok = xs(ok);
        Q1 = prctile(xok,25); Q3 = prctile(xok,75); IQRv = Q3 - Q1;
        iqrOut = (xs < Q1 - IQR_FACTOR*IQRv) | (xs > Q3 + IQR_FACTOR*IQRv);
        zOut   = abs((xs - mean(xok)) / std(xok)) > Z_THRESH;
        combined = iqrOut & zOut & ok;
        x(combined) = NaN;
        staTbl{s}.(nm) = x;
        fprintf('%-8s %-10s %8d %9d %10d\n', stations{s}, nm, ...
            sum(iqrOut & ok), sum(zOut & ok), sum(combined));
    end
end

% --- 2c-ii. Rainfall-specific error screening ----------------------------
% Instead of statistical outlier rules, rainfall is screened for the error
% signatures that actually occur in rain-gauge records: negative values
% (already caught by the range check) and long runs of an identical non-zero
% value, which indicate a stuck gauge or a transcription fault.
fprintf('\n--- 2c-ii. Rainfall error screening (flat-line / stuck-gauge) ---\n');
FLAT_RUN = 5;
rainFlatFlags = cell(1,2);
rainFlatValues = cell(1,2);  % original rainfall values at flagged dates
for s = 1:2
    r = staTbl{s}.Rain;
    flagged = false(size(r));
    i = 1;
    while i <= numel(r)
        if isfinite(r(i)) && r(i) > 0
            j = i;
            while j < numel(r) && isfinite(r(j+1)) && r(j+1) == r(i); j = j + 1; end
            if (j - i + 1) >= FLAT_RUN; flagged(i:j) = true; end
            i = j + 1;
        else
            i = i + 1;
        end
    end
    rainFlatFlags{s} = flagged;
    rainFlatValues{s} = r(flagged);
    if REMOVE_FLATLINE_RAINFALL
        rr = staTbl{s}.Rain; rr(flagged) = NaN; staTbl{s}.Rain = rr;
    end
    fprintf('%-8s : %d days flagged as repeated non-zero runs (>= %d days); removed=%d\n', ...
        stations{s}, sum(flagged), FLAT_RUN, REMOVE_FLATLINE_RAINFALL*sum(flagged));
end
fprintf(['Note: repeated non-zero rainfall runs are flagged but retained by default.\n' ...
         'Set REMOVE_FLATLINE_RAINFALL=true only after confirming gauge/transcription error.\n' ...
         'No IQR/Z-score rule is applied to heavy rainfall observations.\n']);

% --- 2d. Internal consistency: Tmin <= Tmean <= Tmax ---------------------
fprintf('\n--- 2d. Internal consistency (Tmin <= Tmean <= Tmax) ---\n');
for s = 1:2
    tbl = staTbl{s};
    bad = (tbl.Tmin > tbl.Tmean) | (tbl.Tmean > tbl.Tmax) | (tbl.Tmin > tbl.Tmax);
    bad = bad & ~isnan(tbl.Tmin) & ~isnan(tbl.Tmean) & ~isnan(tbl.Tmax);
    fprintf('%-8s : %d days violate ordering -> all three set to NaN\n', stations{s}, sum(bad));
    staTbl{s}.Tmin(bad)  = NaN;
    staTbl{s}.Tmean(bad) = NaN;
    staTbl{s}.Tmax(bad)  = NaN;
end

mks_d = staTbl{1}; hnd_d = staTbl{2};

% --- 2e. Spatial consistency ---------------------------------------------
fprintf('\n--- 2e. Spatial consistency (HND vs MKS, matched dates) ---\n');
[~, im, ih] = intersect(mks_d.Date, hnd_d.Date);
for v = 1:numel(varsQC)
    nm = varsQC{v};
    r = corr(mks_d.(nm)(im), hnd_d.(nm)(ih), 'rows','complete');
    flag = ''; if r < 0.7, flag = '  <-- below 0.7'; end
    fprintf('  r(%s) = %.4f%s\n', nm, r, flag);
end
fprintf(['  Note: rainfall correlation between two stations is expected to be\n' ...
         '  lower than temperature (convective rainfall is spatially patchy);\n' ...
         '  a low value here is a physical result, not a data fault.\n']);

fprintf('\n=== QC COMPLETE ===\n');

%% =========================================================================
%  SECTION 3: Lapse-Rate Correction + Inner Join
% =========================================================================
fprintf('\n=== SECTION 3: LAPSE-RATE CORRECTION + JOIN ===\n');

dz_km = (ELEV_HND - ELEV_MKS) / 1000;
corr_degC = LAPSE_RATE * dz_km;
fprintf('Elevation difference (HND - MKS): %.2f m -> correction %+.4f degC applied to HND\n', ...
    ELEV_HND - ELEV_MKS, corr_degC);
hnd_d.Tmean = hnd_d.Tmean + corr_degC;
hnd_d.Tmin  = hnd_d.Tmin  + corr_degC;
hnd_d.Tmax  = hnd_d.Tmax  + corr_degC;

T = innerjoin(mks_d, hnd_d, 'Keys','Date', ...
    'LeftVariables',  {'Date','Tmean','Tmin','Tmax','Rain','Humidity'}, ...
    'RightVariables', {'Tmean','Tmin','Tmax','Rain','Humidity'});
T.Properties.VariableNames = {'Date', ...
    'Tmean_MKS','Tmin_MKS','Tmax_MKS','Rain_MKS','Hum_MKS', ...
    'Tmean_HND','Tmin_HND','Tmax_HND','Rain_HND','Hum_HND'};

date = T.Date;
fprintf('After inner join: %d days\n', height(T));

%% =========================================================================
%  SECTION 4: Descriptive Statistics                          [Methods 2.3]
% =========================================================================
fprintf('\n=== SECTION 4: DESCRIPTIVE STATISTICS ===\n');
fprintf('%-12s %8s %8s %8s %9s %9s %8s\n', ...
    'Variable','Mean','Median','SD','Skewness','Kurtosis','N');
fprintf('%s\n', repmat('-',1,66));
descVars = {
    'Tmean_MKS', T.Tmean_MKS; 'Tmin_MKS', T.Tmin_MKS; 'Tmax_MKS', T.Tmax_MKS;
    'Rain_MKS',  T.Rain_MKS;  'Hum_MKS',  T.Hum_MKS;
    'Tmean_HND', T.Tmean_HND; 'Tmin_HND', T.Tmin_HND; 'Tmax_HND', T.Tmax_HND;
    'Rain_HND',  T.Rain_HND;  'Hum_HND',  T.Hum_HND };
for i = 1:size(descVars,1)
    x = descVars{i,2}; x = x(~isnan(x));
    fprintf('%-12s %8.3f %8.3f %8.3f %9.3f %9.3f %8d\n', descVars{i,1}, ...
        mean(x), median(x), std(x), skewness(x), kurtosis(x), numel(x));
end

%% =========================================================================
%  SECTION 5: Normality Tests                                 [Methods 2.3]
% =========================================================================
fprintf('\n=== SECTION 5: NORMALITY TESTS ===\n');
fprintf(['Marginal-distribution diagnostics only; these tests do not determine ' ...
    'the normality\nof paired station differences. Lilliefors is the ' ...
    'parameter-estimated normality test.\n']);
fprintf('%-12s %10s %12s %12s %12s   %s\n', ...
    'Variable','KS-stat','KS-p','Lillie-p','Decision','Implication');
fprintf('%s\n', repmat('-',1,86));
for i = 1:size(descVars,1)
    x = descVars{i,2}; x = x(~isnan(x));
    z = (x - mean(x)) / std(x);
    [~, pKS, statKS] = kstest(z);
    try
        [~, pLL] = lillietest(x);
    catch
        pLL = NaN;
    end
    if isfinite(pLL) && pLL < 0.05
        dec = 'non-normal'; imp = 'marginal distribution';
    elseif isfinite(pLL)
        dec = 'not rejected'; imp = 'marginal distribution';
    else
        dec = 'not tested'; imp = 'inspect paired differences';
    end
    fprintf('%-12s %10.4f %12s %12s %12s   %s\n', ...
        descVars{i,1}, statKS, fmtP(pKS), fmtP(pLL), dec, imp);
end

%% =========================================================================
%  SECTION 6: Paired Inter-Station Comparison                 [Methods 2.6]
% =========================================================================
fprintf('\n=== SECTION 6: INTER-STATION COMPARISON (MKS - HND) ===\n');
fprintf('Observations are matched by date, so PAIRED tests are used:\n');
fprintf('paired t-test (parametric) and Wilcoxon signed-rank (nonparametric); reported mean difference is MKS-HND.\n');
fprintf(['Variance equality is reported via vartest2 for completeness. Daily paired-test p-values do not adjust for serial dependence and should be treated as supporting/descriptive evidence, not autocorrelation-robust inference. ' ...
    'For inference,\ninterpret the paired tests on matched dates; the ' ...
    'marginal normality tests above\nare not a substitute for checking ' ...
    'the distribution of paired differences.\n\n']);
fprintf('%-10s %12s %14s %14s %14s\n', ...
    'Variable','mean MKS-HND','paired t p','Wilcoxon p','vartest2 p');
fprintf('%s\n', repmat('-',1,62));
cmpPairs = {
    'Tmean', T.Tmean_HND, T.Tmean_MKS;
    'Tmin',  T.Tmin_HND,  T.Tmin_MKS;
    'Tmax',  T.Tmax_HND,  T.Tmax_MKS;
    'Rain',  T.Rain_HND,  T.Rain_MKS;
    'RH',    T.Hum_HND,   T.Hum_MKS };
for i = 1:size(cmpPairs,1)
    a = cmpPairs{i,2}; b = cmpPairs{i,3};
    ok = ~isnan(a) & ~isnan(b);
    a = a(ok); b = b(ok);
    [~, pT]  = ttest(a, b);              % paired
    pW       = signrank(a, b);           % Wilcoxon signed-rank
    [~, pV]  = vartest2(a, b);
    fprintf('%-10s %12.4f %14s %14s %14s\n', ...
        cmpPairs{i,1}, mean(b-a), fmtP(pT), fmtP(pW), fmtP(pV));
end
fprintf(['NOTE: These paired daily tests do not adjust for serial dependence. ' ...
    'Treat their p-values as supporting/descriptive evidence, not as ' ...
    'autocorrelation-robust inference.\\n\\n']);

%% =========================================================================
%  SECTION 7: Compute DT = T_MKS - T_HND
% =========================================================================
fprintf('\n=== SECTION 7: THERMAL CONTRAST (DT = MKS - HND) ===\n');
DTmean = T.Tmean_MKS - T.Tmean_HND;
DTmin  = T.Tmin_MKS  - T.Tmin_HND;
DTmax  = T.Tmax_MKS  - T.Tmax_HND;

fprintf('DTmean       : mean=%+.3f degC, std=%.3f\n', mean(DTmean,'omitnan'), std(DTmean,'omitnan'));
fprintf('DTmin (night): mean=%+.3f degC, std=%.3f\n', mean(DTmin,'omitnan'),  std(DTmin,'omitnan'));
fprintf('DTmax (day)  : mean=%+.3f degC, std=%.3f\n', mean(DTmax,'omitnan'),  std(DTmax,'omitnan'));
fprintf('SIGN CONVENTION: positive DT = MKS warmer; negative DT = MKS cooler.\n');

%% =========================================================================
%  SECTION 8: RELATIVE Homogeneity Test of DT               [Methods 2.3]
% =========================================================================
fprintf('\n=== SECTION 8: RELATIVE HOMOGENEITY OF DT (SNHT + Pettitt) ===\n');
fprintf('Applied to MONTHLY, DESEASONALIZED DT series.\n');

yr = year(date); mo = month(date);
ym = yr*100 + mo;
uym = unique(ym);

dtSeries = {DTmean, DTmin, DTmax};
dtNames  = {'DTmean','DTmin','DTmax'};
shiftInfo = cell(1,3);
homogResults = repmat(struct('name','','monthly',[],'deseasonalized',[], ...
    'uym',[],'snhtStatistic',NaN,'snhtP',NaN,'snhtIndex',NaN, ...
    'pettittStatistic',NaN,'pettittP',NaN,'pettittIndex',NaN, ...
    'selectedIndex',NaN,'lastPreYM',NaN,'firstPostYM',NaN, ...
    'rawPreMean',NaN,'rawPostMean',NaN,'rawShift',NaN,'significant',false),1,3);

fprintf('\n%-8s %6s %10s %8s %10s %10s %10s %10s   %s\n', ...
    'Series','n_mon','SNHT_T','SNHT_p','Pettitt_U','Pettitt_p','last pre','first post','Conclusion');
fprintf('%s\n', repmat('-',1,105));

for k = 1:3
    x = dtSeries{k};
    monthlyAll = NaN(numel(uym),1);
    monLblAll  = NaN(numel(uym),1);
    for i = 1:numel(uym)
        idx = (ym == uym(i));
        monthlyAll(i) = mean(x(idx), 'omitnan');
        monLblAll(i)  = mod(uym(i),100);
    end
    good = isfinite(monthlyAll);
    monthly = monthlyAll(good);
    monLbl = monLblAll(good);
    uymGood = uym(good);

    % Remove calendar-month climatology.
    deseas = monthly;
    for m = 1:12
        sel = (monLbl == m);
        if any(sel)
            deseas(sel) = monthly(sel) - mean(monthly(sel));
        end
    end

    [Tsnht, kSnht, pSnht] = snht_test(deseas, SNHT_MC);
    [Upett, kPett, pPett] = pettitt_test(deseas);
    sigSNHT = pSnht < 0.05;
    sigPett = pPett < 0.05;

    % Select only a statistically supported location. If both tests are
    % significant, use the earlier of their maxima as a conservative split.
    if sigSNHT && sigPett
        kb = min(kSnht, kPett);
    elseif sigSNHT
        kb = kSnht;
    elseif sigPett
        kb = kPett;
    else
        kb = NaN;
    end

    if isfinite(kb) && kb >= 1 && kb < numel(uymGood)
        lastPreYM = uymGood(kb);
        firstPostYM = uymGood(kb+1);
        rawPreMean = mean(monthly(1:kb), 'omitnan');
        rawPostMean = mean(monthly(kb+1:end), 'omitnan');
        rawShift = rawPostMean - rawPreMean;
        concl = 'INHOMOGENEOUS';
        shiftInfo{k} = sprintf(['  %s: raw monthly pre-mean = %+.3f degC (n=%d mo), ' ...
            'post-mean = %+.3f degC (n=%d mo), shift = %+.3f degC\n'], ...
            dtNames{k}, rawPreMean, kb, rawPostMean, numel(monthly)-kb, rawShift);
        brkLastPre = lastPreYM;
        brkFirstPost = firstPostYM;
        isSig = true;
    else
        lastPreYM = NaN; firstPostYM = NaN;
        rawPreMean = NaN; rawPostMean = NaN; rawShift = NaN;
        concl = 'no significant break';
        shiftInfo{k} = '';
        brkLastPre = NaN; brkFirstPost = NaN; isSig = false;
    end

    homogResults(k).name = dtNames{k};
    homogResults(k).monthly = monthly;
    homogResults(k).deseasonalized = deseas;
    homogResults(k).uym = uymGood;
    homogResults(k).snhtStatistic = Tsnht;
    homogResults(k).snhtP = pSnht;
    homogResults(k).snhtIndex = kSnht;
    homogResults(k).pettittStatistic = Upett;
    homogResults(k).pettittP = pPett;
    homogResults(k).pettittIndex = kPett;
    homogResults(k).selectedIndex = kb;
    homogResults(k).lastPreYM = lastPreYM;
    homogResults(k).firstPostYM = firstPostYM;
    homogResults(k).rawPreMean = rawPreMean;
    homogResults(k).rawPostMean = rawPostMean;
    homogResults(k).rawShift = rawShift;
    homogResults(k).significant = isSig;

    if isSig
        fprintf('%-8s %6d %10.3f %8.4f %10.0f %10.4f %10d %10d   %s\n', ...
            dtNames{k}, numel(deseas), Tsnht, pSnht, Upett, pPett, ...
            brkLastPre, brkFirstPost, concl);
    else
        fprintf('%-8s %6d %10.3f %8.4f %10.0f %10.4f %10s %10s   %s\n', ...
            dtNames{k}, numel(deseas), Tsnht, pSnht, Upett, pPett, '-', '-', concl);
    end
end

% Conservative common cutoff: earliest FIRST-post month among significant series.
validFirstPost = [homogResults.firstPostYM];
validFirstPost = validFirstPost(isfinite(validFirstPost));
if ~isempty(validFirstPost)
    BREAK_YM = min(validFirstPost);
    BREAK_DETECTED = true;
    fprintf('\nConservative cutoff (first post-break month) = %d\n', BREAK_YM);
    fprintf('The pre-break segment ends in the preceding month.\n');
else
    BREAK_DETECTED = false;
    warning('No significant change point detected; the configured cutoff is used for sensitivity only.');
end

fprintf('\n--- Magnitude of detected shifts (raw monthly means) ---\n');
for k = 1:3
    if ~isempty(shiftInfo{k}); fprintf('%s', shiftInfo{k}); end
end
fprintf(['\nThe detected breaks indicate a change in the relative MKS-HND series;\n' ...
    'they do not identify which station or observational factor caused it.\n' ...
    'Check station history as contextual evidence, not as proof of attribution.\n']);

%% =========================================================================
%  SECTION 9: Annual DT Trend, FULL RECORD                     [Table S3]
%
%  WARNING: this is the FULL-RECORD trend. It includes the post-2023
%  discontinuity and reverses sign relative to the pre-discontinuity period
%  estimate, which is what the manuscript reports as Table 3. Table 3 is
%  filled from Section 14b-ii, not from here. This table belongs in
%  Table S3, where both windows are shown side by side.
% =========================================================================
fprintf('\n=== SECTION 9: ANNUAL DT TREND, FULL RECORD (Table S3) ===\n');
uniqYr = unique(yr);
DTmean_a = NaN(numel(uniqYr),1);
DTmin_a  = NaN(numel(uniqYr),1);
DTmax_a  = NaN(numel(uniqYr),1);
for i = 1:numel(uniqYr)
    idx = (yr == uniqYr(i));
    DTmean_a(i) = mean(DTmean(idx),'omitnan');
    DTmin_a(i)  = mean(DTmin(idx), 'omitnan');
    DTmax_a(i)  = mean(DTmax(idx), 'omitnan');
end

S_mean = trendtest(DTmean_a, uniqYr); MK_mean = mann_kendall(DTmean_a);
S_min  = trendtest(DTmin_a,  uniqYr); MK_min  = mann_kendall(DTmin_a);
S_max  = trendtest(DTmax_a,  uniqYr); MK_max  = mann_kendall(DTmax_a);

fprintf('\n%-16s %8s %9s %20s %9s %20s %7s %7s %7s %6s\n', ...
    'Definition','Mean','OLS/yr','OLS 95% CI','Sen/yr','Sen 95% CI','R2','MK-Z','MK p','n');
fprintf('%s\n', repmat('-',1,124));
printTrendRow('DT Tmean',        mean(DTmean_a,'omitnan'), S_mean, MK_mean);
printTrendRow('DT Tmax (day)',   mean(DTmax_a,'omitnan'),  S_max,  MK_max);
printTrendRow('DT Tmin (night)', mean(DTmin_a,'omitnan'),  S_min,  MK_min);

fprintf(['\nNOTE: full-record trends include the post-break segment. Compare their\n' ...
         'direction and magnitude with the complete pre-break trends in Section 14.\n' ...
         'A sign reversal should be reported only if it is present in the computed values.\n']);

%% =========================================================================
% SECTION 10: DAILY MULTIVARIABLE MODEL (RAIN + RH)
% Main model for Section 3.3. DT(t) is related to HND rainfall and RH at
% exact calendar lag 0-3 days, controlling calendar-month fixed effects.
% HAC covariance uses Bartlett weights and exact calendar-day separation.
% Sign convention: DT = MKS - HND; HND temperature corrected once above.
% =========================================================================
fprintf('\n=== SECTION 10: DAILY MULTIVARIABLE MODELS (RAIN + RH) ===\n');
LAGS = 0:3;
HAC_LAG = 5;
dn = datenum(date(:));
mo = month(date(:));
rainH = T.Rain_HND(:);
rhH = T.Hum_HND(:);
responses = {DTmax(:), DTmin(:)};
responseNames = {'DTmax_day','DTmin_night'};
modelRows = cell(0,12);
modelFit = struct();
for ir = 1:2
    for lagDay = LAGS
        [b,ci,p,R2,n] = fitLaggedModel(responses{ir},rainH,rhH,mo,dn,lagDay,HAC_LAG);
        modelRows(end+1,:) = {responseNames{ir},lagDay,n,b(2),ci(2,1),ci(2,2),p(2), ...
            b(3),ci(3,1),ci(3,2),p(3),R2}; %#ok<SAGROW>
        key = sprintf('r%d_l%d',ir,lagDay);
        modelFit.(key) = struct('beta',b,'CI',ci,'p',p,'R2',R2,'n',n);
        fprintf('%s | lag %d | Rain=%+.5f [%+.5f,%+.5f], %s | RH=%+.5f [%+.5f,%+.5f], %s | R2=%.3f, n=%d\n', ...
            responseNames{ir},lagDay,b(2),ci(2,1),ci(2,2),pstr(p(2)), ...
            b(3),ci(3,1),ci(3,2),pstr(p(3)),R2,n);
    end
end
T_model_daily = cell2table(modelRows,'VariableNames', ...
    {'Response','Lag_calendar_days','N_days','Beta_Rain_degC_per_mm', ...
     'Rain_CI95_low','Rain_CI95_high','Rain_HAC_p','Beta_RH_degC_per_RHpoint', ...
     'RH_CI95_low','RH_CI95_high','RH_HAC_p','R2'});

% Supplementary season-stratified lag-0 models: wet Nov-Apr, dry May-Oct.
wetMonths = [11 12 1 2 3 4];
wetMask = ismember(mo,wetMonths);
seasonCases = {'DTmax_day','Wet',DTmax(:),wetMask; ...
               'DTmax_day','Dry',DTmax(:),~wetMask; ...
               'DTmin_night','Wet',DTmin(:),wetMask; ...
               'DTmin_night','Dry',DTmin(:),~wetMask};
seasonRows = cell(0,11);
for q = 1:size(seasonCases,1)
    [bS,ciS,pS,R2S,nS] = fitSeason(seasonCases{q,3},rainH,rhH,mo,dn,seasonCases{q,4},HAC_LAG);
    seasonRows(end+1,:) = {seasonCases{q,1},seasonCases{q,2},nS, ...
        bS(2),ciS(2,1),ciS(2,2),pS(2),bS(3),ciS(3,1),ciS(3,2),pS(3)}; %#ok<SAGROW>
end
T_model_seasonal = cell2table(seasonRows,'VariableNames', ...
    {'Response','Season','N_days','Beta_Rain','Rain_CI95_low','Rain_CI95_high','Rain_HAC_p', ...
     'Beta_RH','RH_CI95_low','RH_CI95_high','RH_HAC_p'});

% RH-only sensitivity is evaluated in Section 14c using the same exact-date HAC routine.
%% =========================================================================
%  SECTION 12: ETCCDI Rainfall Extreme Indices                 [Table 6]
% =========================================================================
fprintf('\n=== SECTION 12: ETCCDI RAINFALL EXTREME INDICES ===\n');
% ETCCDI is calculated from each station's own daily calendar, NOT the
% inter-station inner join. This prevents a missing date at one station
% from unnecessarily removing a valid rainfall observation at the other.
etccdi_MKS = computeETCCDI(mks_d.Date, mks_d.Rain, WET_DAY_THRESH, BASE_YEARS);
etccdi_HND = computeETCCDI(hnd_d.Date, hnd_d.Rain, WET_DAY_THRESH, BASE_YEARS);

fprintf('\n--- Annual index means over record ---\n');
fprintf('%-8s %9s %9s %8s %9s %9s %8s %10s\n', ...
    'Station','RX1day','RX5day','SDII','R95p','R99p','CDD','PRCPTOT');
fprintf('%s\n', repmat('-',1,74));
printETCCDIrow('MKS', etccdi_MKS);
printETCCDIrow('HND', etccdi_HND);

% Percentile thresholds and how often they are actually exceeded. R99p is
% zero in a year when no day passed the threshold, which at one to two
% qualifying days per year is an ordinary outcome rather than missing data.
% Fig. 6c shows several such years at MKS, and a reader who is not told this
% will read them as a gap in the record.
fprintf('\n--- Percentile thresholds and exceedance frequency ---\n');
fprintf('%-8s %10s %10s %14s %14s %12s\n', ...
    'Station','p95 (mm)','p99 (mm)','R95p days/yr','R99p days/yr','yrs R99p=0');
fprintf('%s\n', repmat('-',1,74));
for s = 1:2
    pairsP = {'MKS', etccdi_MKS; 'HND', etccdi_HND};
    E = pairsP{s,2};
    fprintf('%-8s %10.1f %10.1f %14.2f %14.2f %12d\n', pairsP{s,1}, ...
        E.p95, E.p99, mean(E.nR95day), mean(E.nR99day), sum(E.R99p == 0));
end
fprintf(['\nQuote the "R99p days/yr" figure in the Fig. 6 caption, so that the\n' ...
         'zero years are explained by a computed number rather than asserted.\n']);

% Trends are reported with BOTH ordinary least squares and the nonparametric
% Mann-Kendall test. Extreme-value indices are strongly right-skewed and a
% small number of exceptional years carries high leverage in OLS; Mann-Kendall
% depends only on the ordering of values, so agreement between the two is
% evidence that a trend is not driven by a few extreme years.
fprintf('\n--- Trends (OLS and Mann-Kendall) ---\n');
fprintf('%-8s %-9s %10s %20s %7s %8s %8s %5s\n', ...
    'Station','Index','Slope/yr','95% CI','R2','OLS p','MK p','agree');
fprintf('%s\n', repmat('-',1,82));
idxNames = {'RX1day','RX5day','SDII','R95p','R99p','CDD','PRCPTOT'};
for i = 1:numel(idxNames)
    nm = idxNames{i};
    pairs = {'MKS', etccdi_MKS; 'HND', etccdi_HND};
    for s = 1:2
        S  = trendtest(pairs{s,2}.(nm), pairs{s,2}.Year);
        MK = mann_kendall(pairs{s,2}.(nm));
        if (S.p < 0.05) == (MK.p < 0.05); tag = 'yes'; else; tag = 'NO'; end
        fprintf('%-8s %-9s %10.3f  [%+.3f, %+.3f] %7.3f %8.4f %8.4f %5s\n', ...
            pairs{s,1}, nm, S.slope, S.ci(1), S.ci(2), S.R2, S.p, MK.p, tag);
    end
end
fprintf(['\nA row marked NO means OLS and Mann-Kendall disagree on significance.\n' ...
         'For skewed extreme indices this usually indicates that the OLS result\n' ...
         'is driven by a few exceptional years; report the Mann-Kendall result\n' ...
         'and say so explicitly.\n']);

%% =========================================================================
%  SECTION 13: Visualization
% =========================================================================
% Figures are produced in SECTION 16, after the robustness and attribution
% results, so that trend lines can be drawn over the correct period and the
% post-relocation segment can be marked.

%% =========================================================================
%  SECTION 14: ROBUSTNESS CHECK - pre-discontinuity analysis window
%  Repeats key analyses before the conservative cutoff. The cutoff is
%  based on the earliest statistically supported change point when detected;
%  it is not interpreted as the date of instrument replacement or relocation.
% =========================================================================
fprintf('\n=== SECTION 14: ROBUSTNESS - PRE-DISCONTINUITY WINDOW ===\n');
if BREAK_DETECTED
    fprintf('Cutoff is based on a statistically supported breakpoint.\n');
else
    fprintf('WARNING: no breakpoint detected; configured cutoff %d is sensitivity-only.\n',BREAK_YM);
end
ymAll = year(date)*100 + month(date);
sub   = ymAll < BREAK_YM;
fprintf('Full period      : %s to %s (n = %d days)\n', ...
    datestr(min(date),'yyyy-mm'), datestr(max(date),'yyyy-mm'), height(T));
fprintf('Pre-break window : %s to %s (n = %d days, %.1f%% of record)\n', ...
    datestr(min(date),'yyyy-mm'), datestr(max(date(sub)),'yyyy-mm'), ...
    sum(sub), 100*sum(sub)/height(T));

% --- 14a. Mean DT ---------------------------------------------------------
fprintf('\n--- 14a. Mean thermal contrast (degC) ---\n');
fprintf('%-16s %14s %14s %12s\n','Metric','Full period','Pre-break','Difference');
fprintf('%s\n', repmat('-',1,60));
dtAll = {DTmean, DTmax, DTmin};
dtLab = {'DT Tmean','DT Tmax (day)','DT Tmin (night)'};
for k = 1:3
    a = mean(dtAll{k},'omitnan');
    b = mean(dtAll{k}(sub),'omitnan');
    fprintf('%-16s %14.3f %14.3f %12.3f\n', dtLab{k}, a, b, b-a);
end

% --- 14b. Annual trends ---------------------------------------------------
% NOTE ON PARTIAL YEARS: the pre-discontinuity period ends in July 2023, so a
% 2023 annual mean would cover only January-July - a seasonally unbalanced
% partial year at the regression endpoint, where leverage is highest. Trends
% are therefore computed over COMPLETE calendar years only (2009-2022), and
% the partial-year version is printed alongside so the sensitivity is visible.
fprintf('\n--- 14b. Annual DT trend, full vs pre-discontinuity period ---\n');
yrSub  = year(date(sub));
uYrSub = unique(yrSub);
% complete calendar years within the pre-discontinuity period
cntYr  = arrayfun(@(y) numel(unique(month(date(sub & yr == y)))), uYrSub);
uYrSubC = uYrSub(cntYr == 12);
fprintf('Complete years in pre-break window: %d-%d (n=%d); ', ...
    min(uYrSubC), max(uYrSubC), numel(uYrSubC));
fprintf('partial year(s) excluded: %s\n', mat2str(uYrSub(cntYr < 12)'));

fprintf('\n%-16s %10s %8s %10s %8s %10s %8s   %s\n', ...
    'Metric','Full slp','Full OLS p','Pre-break slp','Pre-break OLS p','CompleteYr','OLS p','Same?');
fprintf('%s\n', repmat('-',1,92));
for k = 1:3
    x = dtAll{k};
    aFull = NaN(numel(uniqYr),1);
    for i = 1:numel(uniqYr)
        aFull(i) = mean(x(yr == uniqYr(i)),'omitnan');
    end
    Sf = trendtest(aFull, uniqYr);

    xs = x(sub);
    aSub = NaN(numel(uYrSub),1);
    for i = 1:numel(uYrSub)
        aSub(i) = mean(xs(yrSub == uYrSub(i)),'omitnan');
    end
    Ss = trendtest(aSub, uYrSub);

    aSubC = NaN(numel(uYrSubC),1);
    for i = 1:numel(uYrSubC)
        aSubC(i) = mean(xs(yrSub == uYrSubC(i)),'omitnan');
    end
    Sc = trendtest(aSubC, uYrSubC);

    sameConcl = (Sf.p < 0.05) == (Sc.p < 0.05);
    if sameConcl; tag = 'YES'; else; tag = 'NO  <-- CHECK'; end
    fprintf('%-16s %10.4f %8.4f %10.4f %8.4f %10.4f %8.4f   %s\n', ...
        dtLab{k}, Sf.slope, Sf.p, Ss.slope, Ss.p, Sc.slope, Sc.p, tag);
end
fprintf(['\nThe "CompleteYr" column is the one to report: it excludes the\n' ...
         'seasonally unbalanced partial year. If it differs materially from\n' ...
         'the partial-year estimate, the partial year materially affects the trend estimate.\n']);

fprintf('\n--- 14b-ii. Full detail for the complete-year trends (for Table 3) ---\n');
fprintf('%-16s %9s %20s %9s %20s %7s %7s %8s %6s\n', ...
    'Metric','OLS/yr','OLS 95% CI','Sen/yr','Sen 95% CI','R2','OLS p','MK p','n');
fprintf('%s\n', repmat('-',1,125));
for k = 1:3
    x = dtAll{k}; xs = x(sub);
    aSubC = NaN(numel(uYrSubC),1);
    for i = 1:numel(uYrSubC)
        aSubC(i) = mean(xs(yrSub == uYrSubC(i)),'omitnan');
    end
    S  = trendtest(aSubC, uYrSubC);
    % Mann-Kendall on the SAME series Table 3 reports. Section 9 also computes
    % it, but on the full record - a different quantity, and not the one the
    % table prints. Without this the seasonal rows of Table 3 carry MK and the
    % annual rows do not, which a reviewer will ask about.
    MK = mann_kendall(aSubC);
    fprintf('%-16s %9.4f [%+.4f,%+.4f] %9.4f [%+.4f,%+.4f] %7.3f %7.4f %8.4f %6d\n', ...
        dtLab{k}, S.slope, S.ci(1), S.ci(2), S.senSlope, ...
        S.senCI(1), S.senCI(2), S.R2, S.p, MK.p, S.n);
end

% --- 14c. Lag-0 RH regression ------------------------------------------------
fprintf('\n--- 14c. Lag-0 regression, DT ~ RH (HND predictor) ---\n');
fprintf('%-16s %10s %10s %10s %10s   %s\n', ...
    'Model / term','Full beta','Full HACp','Pre-break beta','Pre-break HACp','Same sign & sig?');
fprintf('%s\n', repmat('-',1,78));
rhPairs = {
    'DTmax', DTmax, T.Hum_HND;
    'DTmin', DTmin, T.Hum_HND };
for k = 1:2
    y = rhPairs{k,2}; H = rhPairs{k,3};
    [bF,~,pF] = fitRH(y,H,mo,dn,0,HAC_LAG);
    [bS,~,pS] = fitRH(y(sub),H(sub),mo(sub),dn(sub),0,HAC_LAG);
    okSign = sign(bF) == sign(bS);
    okSig  = (pF < 0.05) == (pS < 0.05);
    if okSign && okSig; tag = 'YES'; else; tag = 'NO  <-- CHECK'; end
    fprintf('%-16s %10.4f %10s %10.4f %10s   %s\n', ...
        sprintf('%s / RH', rhPairs{k,1}), bF, fmtP(pF), bS, fmtP(pS), tag);
end

% --- 14d. ETCCDI trends ---------------------------------------------------
fprintf('\n--- 14d. ETCCDI trends: OLS + Mann-Kendall, full vs pre-discontinuity ---\n');
% Keep the percentile reference thresholds fixed at the full-record values.
% This isolates the effect of changing the analysis window; it does not
% simultaneously change the R95p/R99p reference distribution.
fixedThresholdsMKS = struct('p95',etccdi_MKS.p95,'p99',etccdi_MKS.p99);
fixedThresholdsHND = struct('p95',etccdi_HND.p95,'p99',etccdi_HND.p99);
mksPre = year(mks_d.Date)*100 + month(mks_d.Date) < BREAK_YM;
hndPre = year(hnd_d.Date)*100 + month(hnd_d.Date) < BREAK_YM;
etSubMKS = computeETCCDI(mks_d.Date(mksPre), mks_d.Rain(mksPre), ...
                         WET_DAY_THRESH, BASE_YEARS, fixedThresholdsMKS);
etSubHND = computeETCCDI(hnd_d.Date(hndPre), hnd_d.Rain(hndPre), ...
                         WET_DAY_THRESH, BASE_YEARS, fixedThresholdsHND);

% Use the SAME complete calendar years at both stations. This prevents a
% partial endpoint year (e.g. Jan-Jul 2023) from entering any pre-break trend,
% permutation, ACF, or Hamed-Rao calculation. It also prevents the two
% stations from being compared over different sets of years if data gaps occur.
completeMKS = etSubMKS.Year(etSubMKS.isCompleteYear);
completeHND = etSubHND.Year(etSubHND.isCompleteYear);
commonPreYears = intersect(completeMKS, completeHND, 'stable');
if numel(commonPreYears) < 3
    error('ETCCDI:InsufficientCommonYears', ...
        'Fewer than three common complete pre-discontinuity years are available.');
end
etSubMKS = subsetETCCDIYears(etSubMKS, commonPreYears);
etSubHND = subsetETCCDIYears(etSubHND, commonPreYears);
fprintf('Pre-discontinuity ETCCDI years used at BOTH stations: %d-%d (n=%d)\n', ...
    commonPreYears(1), commonPreYears(end), numel(commonPreYears));
if any(diff(commonPreYears) ~= 1)
    fprintf('WARNING: common complete years are not contiguous: %s\n', mat2str(commonPreYears(:)'));
end

fprintf(['OLS slope/CI and Theil-Sen slope/CI are reported separately.\n' ...
         'Mann-Kendall tests monotonic trend with tie-corrected variance.\n' ...
         'Pre-discontinuity annual analysis uses the same common complete calendar years at both stations.\n']);
fprintf(['\n%-8s %-9s %9s %9s %9s %9s %9s %9s %9s %9s %5s %5s\n'], ...
    'Station','Index','FullSlp','FullOLS','FullMK','PreSlp','PreOLS','PreMK', ...
    'FullZ','PreZ','OLS?','MK?');
fprintf('%s\n', repmat('-',1,112));
for i = 1:numel(idxNames)
    nm = idxNames{i};
    pairs = {'MKS', etccdi_MKS, etSubMKS; 'HND', etccdi_HND, etSubHND};
    for s = 1:2
        % Full-period estimates and tests
        Sf  = trendtest(pairs{s,2}.(nm), pairs{s,2}.Year);
        MKf = mann_kendall(pairs{s,2}.(nm));
        % Pre-discontinuity estimates and tests
        Sh  = trendtest(pairs{s,3}.(nm), pairs{s,3}.Year);
        MKh = mann_kendall(pairs{s,3}.(nm));

        sameOLS = (Sf.p < 0.05) == (Sh.p < 0.05);
        sameMK  = (MKf.p < 0.05) == (MKh.p < 0.05);
        if sameOLS; tagOLS = 'YES'; else; tagOLS = 'NO'; end
        if sameMK;  tagMK  = 'YES'; else; tagMK  = 'NO'; end
        fprintf('%-8s %-9s %9.3f %9.4g %9.4g %9.3f %9.4g %9.4g %9.3f %9.3f %5s %5s\n', ...
            pairs{s,1}, nm, Sf.slope, Sf.p, MKf.p, Sh.slope, Sh.p, MKh.p, ...
            MKf.Z, MKh.Z, tagOLS, tagMK);
    end
end
fprintf(['\nFullSlp and PreSlp are OLS slopes (mm/year); Sen slopes and rank-based CIs\n' ...
         'are provided in the separate pre-discontinuity validation table.\n' ...
         'FullOLS/PreOLS are OLS p-values; FullMK/PreMK are Mann-Kendall p-values;\n' ...
         'FullZ/PreZ are MK standardized Z.\n' ...
         'OLS? and MK? indicate whether each test gives the same significance\n' ...
         'classification (p < 0.05) in the two analysis periods.\n']);

fprintf(['\nINTERPRETATION: a row marked YES indicates that the significance\n' ...
         'classification is unchanged when the analysis is restricted to the\n' ...
         'pre-discontinuity period. This sensitivity result does not establish\n' ...
         'that station relocation caused the detected relative discontinuity.\n' ...
         'Any difference in OLS or MK significance classification between periods\n' ...
         'should be reported and discussed explicitly.\n']);


% --- 14d-i. Small-sample permutation MK and FDR-BH -----------------------
% Permutation p-values are an IID/exchangeability sensitivity only. They do
% not preserve serial dependence and must not replace the primary MK test
% when autocorrelation is present.
B_MK_PERM = 9999;
rng(20261110,'twister'); % fixed seed for reproducible Monte Carlo permutation p-values
fprintf('\n--- 14d-i. Pre-discontinuity MK permutation and FDR-BH ---\n');
fprintf(['Permutation test: %d random permutations, two-sided, ties preserved; ', ...
         'IID/exchangeability assumption.\n'], B_MK_PERM);
fprintf('BH adjustment is applied separately to the 7 ETCCDI indices within each station.\n');
fprintf('%-8s %-9s %9s %9s %9s %9s %9s\n', ...
    'Station','Index','MK_p_asym','MK_p_perm','BH_q','SenSlope','SenCI excludes 0?');
fprintf('%s\n', repmat('-',1,82));
for sIdx = 1:2
    if sIdx==1
        st='MKS'; E=etSubMKS;
    else
        st='HND'; E=etSubHND;
    end
    pRaw=NaN(numel(idxNames),1);
    pPerm=NaN(numel(idxNames),1);
    senVals=NaN(numel(idxNames),1);
    senCIs=NaN(numel(idxNames),2);
    for ii=1:numel(idxNames)
        nm=idxNames{ii};
        y=E.(nm);
        MK=mann_kendall(y);
        pRaw(ii)=MK.p;
        pPerm(ii)=mk_permutation_p(y,B_MK_PERM);
        S=trendtest(y,E.Year);
        senVals(ii)=S.senSlope;
        senCIs(ii,:)=S.senCI;
    end
    qBH=bh_adjust(pRaw);
    for ii=1:numel(idxNames)
        if all(isfinite(senCIs(ii,:)))
            ciFlag = (senCIs(ii,1)>0 || senCIs(ii,2)<0);
        else
            ciFlag = false;
        end
        fprintf('%-8s %-9s %9.4g %9.4g %9.4g %9.3f %9s\n', ...
            st,idxNames{ii},pRaw(ii),pPerm(ii),qBH(ii),senVals(ii), ...
            string(ciFlag));
    end
end
fprintf(['Permutation p-values assume exchangeable annual values; they are not ', ...
         'autocorrelation-adjusted. Interpret alongside the lag-1 diagnostics.\n' ...
         'The fixed RNG seed makes the Monte Carlo results reproducible.\n']);
fprintf(['BH q-values control the expected false-discovery rate within each ', ...
         'station-specific family of seven indices (pre-discontinuity period).\n']);

% --- 14d-ii. Lag-1 serial-correlation diagnostic -------------------------
% Diagnose serial dependence in linearly detrended annual-index series.
% This is a screening diagnostic, not a substitute for a modified MK test.
fprintf(['\n--- 14d-ii. Lag-1 autocorrelation of detrended annual indices ---\n' ...
         'Residuals are obtained after removing a linear time trend.\n' ...
         'The reported p-value is an approximate two-sided test of lag-1 correlation.\n' ...
         'A significant result (p < 0.05) flags the standard MK result for review.\n']);
fprintf('%-8s %-9s %9s %9s %9s %9s %8s %8s\n', ...
    'Station','Index','Full_r1','Full_p','Pre_r1','Pre_p','FullFlag','PreFlag');
fprintf('%s\n', repmat('-',1,82));
for i = 1:numel(idxNames)
    nm = idxNames{i};
    pairs = {'MKS', etccdi_MKS, etSubMKS; 'HND', etccdi_HND, etSubHND};
    for sIdx = 1:2
        ACFf = detrended_lag1(pairs{sIdx,2}.(nm), pairs{sIdx,2}.Year);
        ACFh = detrended_lag1(pairs{sIdx,3}.(nm), pairs{sIdx,3}.Year);
        if ACFf.p < 0.05; flagF = 'CHECK'; else; flagF = 'no'; end
        if ACFh.p < 0.05; flagH = 'CHECK'; else; flagH = 'no'; end
        fprintf('%-8s %-9s %9.3f %9.4f %9.3f %9.4f %8s %8s\n', ...
            pairs{sIdx,1}, nm, ACFf.r1, ACFf.p, ACFh.r1, ACFh.p, flagF, flagH);
    end
end
fprintf(['\nInterpretation: non-significant lag-1 correlation does not prove independence,\n' ...
         'especially with short annual records. Significant correlation indicates that\n' ...
         'standard MK p-values may require a serial-correlation-adjusted sensitivity\n' ...
         'test (e.g., a modified Mann-Kendall procedure). Do not automatically apply\n' ...
         'pre-whitening without checking its effect on the trend magnitude.\n' ...
         'In this dataset, HND RX1day and RX5day, and SDII at both stations,\n' ...
         'show significant lag-1 dependence; interpret their standard MK\n' ...
         'p-values cautiously and report the sensitivity diagnostics.\n']);

% --- 14d-iii. Hamed-Rao modified MK sensitivity diagnostic ---------------
% v8 uses tie-corrected Var(S) in both standard MK and Hamed-Rao. The
% detrended-rank ACF is calculated with a conventional full-series
% autocovariance denominator (not correlation between separately centered
% lagged segments). Approximate +/-1.96/sqrt(n) bounds are used only as a
% screening rule; because n is short, HR remains a sensitivity diagnostic.
fprintf('\n--- 14d-iii. Hamed-Rao modified MK sensitivity diagnostic ---\n');
fprintf(['Method: Sen-slope detrending; standard ACF of detrended ranks; max lag=3.\n' ...
         'Approximate screening bounds +/-1.96/sqrt(n); tie-corrected Var(S).\n' ...
         'HR p-values are withheld if the variance factor is non-positive.\n' ...
         'Extreme correction factors are reported but flagged for caution.\n']);
flaggedSeries = {'HND','RX5day'; 'MKS','SDII'; 'HND','SDII'};
fprintf('%-8s %-9s %-7s %3s %8s %8s %8s %8s %8s %8s %8s %s\n', ...
    'Station','Index','Period','N','MK_p','HR_p','VarRatio','rho1','rho2','rho3','ACFbound','Status');
fprintf('%s\n', repmat('-',1,116));
for q = 1:size(flaggedSeries,1)
    st = flaggedSeries{q,1}; nm = flaggedSeries{q,2};
    if strcmp(st,'MKS')
        Efull = etccdi_MKS; Ehom = etSubMKS;
    else
        Efull = etccdi_HND; Ehom = etSubHND;
    end
    seriesSet = {Efull,'Full'; Ehom,'Pre-break'};
    for q2 = 1:2
        E = seriesSet{q2,1}; per = seriesSet{q2,2};
        yv = E.(nm);
        MKstd = mann_kendall(yv);
        MKhr = mann_kendall_hamed_rao_v8(yv);
        rho = MKhr.rankACF;
        r1=NaN; r2=NaN; r3=NaN;
        if numel(rho)>=1; r1=rho(1); end
        if numel(rho)>=2; r2=rho(2); end
        if numel(rho)>=3; r3=rho(3); end
        if isfinite(MKhr.p)
            if MKhr.extremeFactor
                status='CAUTION: extreme factor';
            else
                status='Computed; sensitivity only';
            end
        else
            status=['NOT COMPUTABLE: ' MKhr.reason];
        end
        fprintf('%-8s %-9s %-7s %3d %8.4g %8.4g %8.3f %8.3f %8.3f %8.3f %8.3f %s\n', ...
            st,nm,per,MKhr.n,MKstd.p,MKhr.p,MKhr.varRatio,r1,r2,r3,MKhr.acfBound,status);
    end
end
fprintf(['\nInterpretation: the Hamed-Rao result is a sensitivity analysis, not a\n' ...
         'replacement for the pre-specified standard MK. A correction factor\n' ...
         'far below or above 1 is especially unstable for n=14-17 and must not\n' ...
         'be used alone to claim or reject a trend. Report the standard MK,\n' ...
         'Sen slope and period sensitivity; discuss HR only as a diagnostic.\n']);

%% =========================================================================
%  SECTION 15: STATION-SPECIFIC TREND COMPARISON
%  The DT trend is computed on the pre-discontinuity window only. Comparing the
%  station-specific trends describes which series changes more rapidly, but
%  does not establish a physical driver or causal attribution.
% =========================================================================
fprintf('\n=== SECTION 15: STATION-SPECIFIC TRENDS (pre-break period) ===\n');
fprintf('Period: %s to %s\n\n', datestr(min(date),'yyyy-mm'), datestr(max(date(sub)),'yyyy-mm'));

attrVars = {
    'Tmean_MKS', T.Tmean_MKS; 'Tmean_HND', T.Tmean_HND;
    'Tmax_MKS',  T.Tmax_MKS;  'Tmax_HND',  T.Tmax_HND;
    'Tmin_MKS',  T.Tmin_MKS;  'Tmin_HND',  T.Tmin_HND };

fprintf('Complete pre-break calendar years only (%d-%d, n=%d)\n\n', ...
    min(uYrSubC), max(uYrSubC), numel(uYrSubC));
fprintf('%-12s %9s %20s %7s %7s %8s %4s\n', ...
    'Series','Slope/yr','95% CI','R2','OLS p','MK p','sig');
fprintf('%s\n', repmat('-',1,74));
attrSlope = NaN(size(attrVars,1),1);
for i = 1:size(attrVars,1)
    x  = attrVars{i,2};
    xs = x(sub);
    a  = NaN(numel(uYrSubC),1);
    for j = 1:numel(uYrSubC)
        a(j) = mean(xs(yrSub == uYrSubC(j)),'omitnan');
    end
    S  = trendtest(a, uYrSubC);
    MK = mann_kendall(a);
    attrSlope(i) = S.slope;
    star = '';
    if S.p < 0.001; star = '***'; elseif S.p < 0.01; star = '**'; elseif S.p < 0.05; star = '*'; end
    fprintf('%-12s %9.4f  [%+.4f, %+.4f] %7.3f %7.4f %8.4f %4s\n', ...
        attrVars{i,1}, S.slope, S.ci(1), S.ci(2), S.R2, S.p, MK.p, star);
end

fprintf('\n--- Decomposition of the DT trend ---\n');
fprintf('%-16s %12s %12s %14s\n','Metric','MKS slope','HND slope','DT slope (M-H)');
fprintf('%s\n', repmat('-',1,58));
fprintf('%-16s %12.4f %12.4f %14.4f\n','Tmean', attrSlope(1), attrSlope(2), attrSlope(1)-attrSlope(2));
fprintf('%-16s %12.4f %12.4f %14.4f\n','Tmax',  attrSlope(3), attrSlope(4), attrSlope(3)-attrSlope(4));
fprintf('%-16s %12.4f %12.4f %14.4f\n','Tmin',  attrSlope(5), attrSlope(6), attrSlope(5)-attrSlope(6));

fprintf(['\nREADING THIS TABLE (DESCRIPTIVE, NOT CAUSAL ATTRIBUTION):\n' ...
    '  - DT slope is defined as MKS slope minus HND slope.\n' ...
    '  - A positive DT slope means the MKS temperature series increases\n' ...
    '    relative to HND over the selected pre-discontinuity period.\n' ...
    '  - Compare station-specific slopes to describe which component changes\n' ...
    '    more rapidly; do not infer a physical cause from these trends alone.\n' ...
    '  - Absolute station trends include regional climate variability and\n' ...
    '    possible residual station-specific artefacts. The difference series\n' ...
    '    removes only the component shared by both stations, not all confounding.\n' ...
    '  - In particular, a significant DTmean trend does not make the separate\n' ...
    '    DTmin or DTmax trends significant. Report each test independently.\n']);

%% =========================================================================
%  SECTION 16: FIGURES
%
%  Fig. A - annual DT series, showing full record and pre-break trend fits.
%  Integrated Fig. 3-5 below show breakpoints, monthly climatology and coupling.
%  Fig. B - Figure 6: selected indices representing daily maximum,
%           multi-day maximum, and upper-tail rainfall accumulation.
%
%  COLOURS
%  -------
%  The station pair was green/orange, which is the hardest combination for
%  the most common colour-vision deficiencies and the one that collapses
%  first in greyscale printing. Blue and vermillion separate under every
%  common CVD type, differ in lightness so they survive greyscale, and match
%  the palette used in Figs. 3, 4 and 5 - so one colour means one thing
%  across the whole manuscript.
%
%  P-VALUES
%  --------
%  Printed through pstr(), which renders anything below 0.001 as "p < 0.001"
%  rather than "p = 0.000". A p-value is never exactly zero, and a reviewer
%  who sees one written that way has been handed a reason to doubt the rest
%  of the arithmetic.
% =========================================================================
fprintf('\n=== SECTION 16: FIGURES ===\n');

SAVE_FIGS = true;    % export publication figures at 300 dpi
FIG_DIR = 'Output_Figures';

% --- Output housekeeping: archive previous figures, then start clean ---
% Existing outputs are NOT deleted. They are moved to a timestamped archive
% beside Output_Figures so the current folder contains only this run's files.
if exist(FIG_DIR,'dir')
    oldItems = dir(FIG_DIR);
    oldItems = oldItems(~ismember({oldItems.name},{'.','..'}));
    if ~isempty(oldItems)
        if ~exist('Output_Figures_Archive','dir')
            mkdir('Output_Figures_Archive');
        end
        archiveStamp = char(datetime('now','Format','yyyyMMdd_HHmmss'));
        archiveDir = fullfile('Output_Figures_Archive', ...
            ['Output_Figures_' archiveStamp]);
        [mkOK,mkMsg] = mkdir(archiveDir);
        if ~mkOK
            error('Cannot create figure archive folder: %s',mkMsg);
        end
        for iOld = 1:numel(oldItems)
            oldPath = fullfile(FIG_DIR,oldItems(iOld).name);
            [mvOK,mvMsg] = movefile(oldPath,archiveDir);
            if ~mvOK
                error('Could not archive previous figure output "%s": %s', ...
                    oldItems(iOld).name,mvMsg);
            end
        end
        fprintf('Previous figure outputs archived to: %s\n',archiveDir);
    end
else
    mkdir(FIG_DIR);
end

% Main publication output contains the six current JPGs only.
% No .fig files are generated unless explicitly added in a future revision.

cMKS = [0.00 0.45 0.70];    % blue       - coastal harbour
cHND = [0.84 0.37 0.00];    % vermillion - inland airport

% ---------------- Figure A: inter-station thermal contrast ---------------
figA = figure('Name','Fig A - Annual thermal contrast','NumberTitle','off', ...
              'Position',[80 80 1100 470], 'Color','w');
hold on;

% Fixed vertical range requested for a more readable comparison.
% All three annual thermal-contrast series are shown on the same -2 to +2 °C scale.
yLo = -2;
yHi =  2;
% Shade the calendar year containing the detected discontinuity/transition.
% Annual shading is a coarse visual cue, not an exact monthly breakpoint marker.
brkYear = floor(BREAK_YM/100);
if BREAK_DETECTED
    patch([brkYear - 0.5, max(uniqYr)+0.5, max(uniqYr)+0.5, brkYear - 0.5], ...
          [yLo yLo yHi yHi], [0.93 0.93 0.93], 'EdgeColor','none', ...
          'HandleVisibility','off');
    text(brkYear + 0.1, yHi - 0.12, ' detected discontinuity / transition', ...
         'FontSize', 8, 'Color', [0.35 0.35 0.35], 'FontAngle','italic');
else
    text(0.02,0.96,'No significant breakpoint; cutoff shown as configured sensitivity window', ...
        'Units','normalized','FontSize',8,'Color',[0.35 0.35 0.35], ...
        'VerticalAlignment','top');
end

dtSeriesA = {DTmean, DTmax, DTmin};
% Consistent series colours: neutral grey for the mean, vermillion for
% daytime contrast, blue for nighttime contrast.
dtColors  = {[0.25 0.25 0.25], [0.84 0.37 0.00], [0.00 0.45 0.70]};
dtMarks   = {'-o','-s','-^'};
dtLabels  = {'\DeltaT_{mean}','\DeltaT_{max} (day)','\DeltaT_{min} (night)'};

for k = 1:3
    x = dtSeriesA{k};

    aAll = NaN(numel(uniqYr),1);
    for i = 1:numel(uniqYr)
        aAll(i) = mean(x(yr == uniqYr(i)), 'omitnan');
    end
    plot(uniqYr, aAll, dtMarks{k}, 'Color', dtColors{k}, 'LineWidth', 1.4, ...
         'MarkerFaceColor', dtColors{k}, 'MarkerSize', 5, ...
         'DisplayName', dtLabels{k});

    xs = x(sub);
    aC = NaN(numel(uYrSubC),1);
    for i = 1:numel(uYrSubC)
        aC(i) = mean(xs(yrSub == uYrSubC(i)), 'omitnan');
    end
    S = trendtest(aC, uYrSubC);
    yfit = mean(aC,'omitnan') + S.slope*(uYrSubC - mean(uYrSubC));
    if S.p < 0.05
        lw = 2.2; ls = '--';
        lbl = sprintf('  trend %+.3f °C yr^{-1} (%s)', S.slope, pstr(S.p));
    else
        lw = 1.1; ls = ':';
        lbl = sprintf('  trend %+.3f °C yr^{-1} (n.s.)', S.slope);
    end
    plot(uYrSubC, yfit, ls, 'Color', dtColors{k}, 'LineWidth', lw, ...
         'DisplayName', lbl);
end

yline(0, 'k:', 'LineWidth', 1, 'HandleVisibility','off');
grid on; box on;
ylim([yLo yHi]); xlim([min(uniqYr)-0.5, max(uniqYr)+0.5]);
set(gca,'YTick',-2:0.5:2,'FontSize',10,'LineWidth',0.9,'GridAlpha',0.18);
xlabel('Year'); ylabel('\DeltaT = T_{MKS} - T_{HND} (\circC)');
title(sprintf('Annual thermal contrast (%d-%d); pre-break trend fitted to %d-%d', ...
    min(uniqYr), max(uniqYr), min(uYrSubC), max(uYrSubC)));
legend('Location','eastoutside','FontSize',8);

% ---------------- Figure B: extreme rainfall indices (Fig. 6) ------------
% Journal-style layout: three wide plotting panels and colour-coded
% statistical cards INSIDE the figure canvas. The cards sit in a dedicated
% right-hand column and therefore never cover the observations or trends.
figB = figure('Name','Fig 6 - Extreme Rainfall Indices','NumberTitle','off', ...
              'Position',[50 40 1600 1120], 'Color','w');

idx3  = {'RX1day','RX5day','R99p'};
ylab3 = {'RX1day (mm)','RX5day (mm)','R99p (mm)'};
ttl3  = {'Annual maximum 1-day rainfall', ...
         'Annual maximum 5-day rainfall', ...
         'Rainfall from extremely wet days (>99th percentile)'};

rowY = [0.695 0.385 0.090];
rowH = 0.235;
plotX = 0.075; plotW = 0.665;
statX = 0.765; statW = 0.215;
blueTint = [0.91 0.96 1.00];
orangeTint = [1.00 0.95 0.89];
blueDark = [0.00 0.25 0.48];
orangeDark = [0.55 0.20 0.00];

for pnl = 1:3
    ax = axes('Parent',figB,'Units','normalized', ...
              'Position',[plotX rowY(pnl) plotW rowH]);
    hold(ax,'on'); box(ax,'on');

    nm = idx3{pnl};
    sets = {etccdi_MKS, cMKS, 'MKS (coastal)'; ...
            etccdi_HND, cHND, 'HND (airport)'};

    % Handles are retained for a compact legend in panel (a).
    hAnnual = gobjects(2,1);
    hTrend = gobjects(2,1);
    stats = cell(2,1);

    for ss = 1:2
        E = sets{ss,1};
        c = sets{ss,2};
        yy = E.(nm);
        xx = E.Year(:);
        valid = isfinite(xx) & isfinite(yy);
        xxv = xx(valid); yyv = yy(valid);

        hAnnual(ss) = plot(ax,xxv,yyv,'-o','Color',c,'LineWidth',1.55, ...
            'MarkerFaceColor','w','MarkerEdgeColor',c,'MarkerSize',5.5);

        S = trendtest(yy,xx);
        MK = mann_kendall(yy);
        yfit = mean(yyv,'omitnan') + S.slope*(xxv-mean(xxv,'omitnan'));

        if isfinite(S.p) && S.p < 0.05
            trendStyle = '--'; trendWidth = 2.5;
        else
            trendStyle = ':'; trendWidth = 1.8;
        end
        hTrend(ss) = plot(ax,xxv,yfit,trendStyle,'Color',c,'LineWidth',trendWidth);

        stats{ss} = struct('slope',S.slope,'ci',S.ci, ...
            'pOLS',S.p,'pMK',MK.p);
    end

    grid(ax,'on'); box(ax,'on');
    set(ax,'Layer','top','FontSize',10.5,'TickDir','out', ...
        'LineWidth',0.8,'GridAlpha',0.18);
    ylabel(ax,ylab3{pnl});
    allETyears=[etccdi_MKS.Year(:);etccdi_HND.Year(:)];
    xlim(ax,[min(allETyears)-0.5,max(allETyears)+0.5]);
    title(ax,sprintf('(%s) %s',char('a'+pnl-1),ttl3{pnl}), ...
        'FontWeight','bold','FontSize',12);

    if pnl == 3
        xlabel(ax,'Year');
    else
        set(ax,'XTickLabel',[]);
    end

    if pnl == 1
        lg = legend(ax,[hAnnual(1),hAnnual(2),hTrend(1),hTrend(2)], ...
            {'MKS (coastal)','HND (airport)','MKS trend','HND trend'}, ...
            'Location','northwest','FontSize',8.5,'Box','on');
        lg.Color = 'w';
    end

    % Statistics card axes: all cards are contained within the figure,
    % aligned to the corresponding rainfall panel.
    axStat = axes('Parent',figB,'Units','normalized', ...
        'Position',[statX rowY(pnl) statW rowH], ...
        'XLim',[0 1],'YLim',[0 1],'Visible','off');
    hold(axStat,'on');

    cardPos = [0.015 0.525 0.97 0.455; ...
               0.015 0.020 0.97 0.455];
    cardFill = {blueTint,orangeTint};
    cardEdge = {cMKS,cHND};
    cardDark = {blueDark,orangeDark};

    for ss = 1:2
        st = stats{ss};
        rectangle(axStat,'Position',cardPos(ss,:), ...
            'Curvature',[0.08 0.12], ...
            'FaceColor',cardFill{ss},'EdgeColor',cardEdge{ss}, ...
            'LineWidth',1.25);

        xText = 0.075;
        yTop = cardPos(ss,2)+cardPos(ss,4)-0.060;
        text(axStat,xText,yTop,sets{ss,3}, ...
            'FontSize',9.8,'FontWeight','bold','Color',cardDark{ss}, ...
            'VerticalAlignment','top','Interpreter','none');
        text(axStat,xText,yTop-0.105, ...
            sprintf('Slope: %+.2f mm yr^{-1}',st.slope), ...
            'FontSize',8.8,'Color',[0.08 0.08 0.08], ...
            'VerticalAlignment','top','Interpreter','tex');
        text(axStat,xText,yTop-0.205, ...
            sprintf('95%% CI: [%+.2f, %+.2f]',st.ci(1),st.ci(2)), ...
            'FontSize',8.8,'Color',[0.08 0.08 0.08], ...
            'VerticalAlignment','top','Interpreter','tex');
        text(axStat,xText,yTop-0.300, ...
            sprintf('OLS %s | MK %s',pstr(st.pOLS),pstr(st.pMK)), ...
            'FontSize',8.5,'Color',[0.08 0.08 0.08], ...
            'VerticalAlignment','top','Interpreter','tex');
    end
end

% Caption is part of the figure and explains symbols, period, and inference.
trendPeriodLabel = sprintf('%d-%d',min(allETyears),max(allETyears));
annotation(figB,'textbox',[0.055 0.004 0.91 0.032], ...
    'String',['NOTE: Annual trend period = ' trendPeriodLabel '. ' ...
    'Blue = MKS (coastal); orange = HND (airport). ' ...
    'Solid line + open circles = annual index. Dashed trend = OLS p < 0.05; ' ...
    'dotted trend = OLS p >= 0.05. Slope is in mm yr^{-1}; CI = 95% confidence interval; ' ...
    'OLS = ordinary least squares; MK = Mann-Kendall.'], ...
    'Interpreter','tex','EdgeColor',[0.45 0.45 0.45], ...
    'BackgroundColor',[0.98 0.98 0.98],'LineWidth',0.8, ...
    'HorizontalAlignment','center','VerticalAlignment','middle', ...
    'FontSize',8.0,'Margin',4);

% ---------------- Optional export at journal resolution ------------------
if SAVE_FIGS
    % Save publication figures as JPG only, at 300 dpi.
    exportgraphics(figA,fullfile(FIG_DIR,'Fig1_Annual_Thermal_Contrast.jpg'), ...
        'Resolution',300,'ContentType','image');
    exportgraphics(figB,fullfile(FIG_DIR,'Fig6_Rainfall_Extremes.jpg'), ...
        'Resolution',300,'ContentType','image');

    fprintf('Main figures saved as JPG (300 dpi) in %s.\n',FIG_DIR);
else
    fprintf('Figures displayed on screen only. Set SAVE_FIGS = true to export.\n');
end

%% =========================================================================
% SECTION 16A: PUBLICATION FIGURE - MONTHLY DESEASONALIZED BREAKS (FIG. 3)
% =========================================================================
fig3 = figure('Color','w','Position',[80 50 980 820]);
tl3 = tiledlayout(3,1,'Padding','compact','TileSpacing','compact');
plotOrder = [1 3 2];
plotNames = {'\DeltaT_{mean}','\DeltaT_{max} (day)','\DeltaT_{min} (night)'};
plotCols = {[0.25 0.25 0.25],[0.84 0.37 0.00],[0.00 0.45 0.70]};
for ip = 1:3
    R = homogResults(plotOrder(ip));
    ax = nexttile(tl3); hold(ax,'on'); box(ax,'on');
    tt = datetime(floor(R.uym/100),mod(R.uym,100),15);
    plot(ax,tt,R.deseasonalized,'-','Color',[0.45 0.45 0.45],'LineWidth',0.8);
    if R.significant
        pre = R.uym <= R.lastPreYM; post = R.uym >= R.firstPostYM;
        plot(ax,tt(pre),R.deseasonalized(pre),'o','Color',plotCols{ip}, ...
            'MarkerFaceColor',plotCols{ip},'MarkerSize',3.5);
        plot(ax,tt(post),R.deseasonalized(post),'o','Color',[0.75 0.25 0.15], ...
            'MarkerFaceColor',[0.75 0.25 0.15],'MarkerSize',3.5);
        xline(ax,datetime(floor(R.firstPostYM/100),mod(R.firstPostYM,100),1),'--k');
        msg=sprintf('last pre %d | raw shift %+.3f ^oC | SNHT %s | Pettitt %s', ...
            R.lastPreYM,R.rawShift,pstr(R.snhtP),pstr(R.pettittP));
        text(ax,0.01,0.94,msg,'Units','normalized','VerticalAlignment','top', ...
            'FontSize',8,'BackgroundColor','w','Margin',2);
    else
        text(ax,0.01,0.94,'No statistically supported break','Units','normalized', ...
            'VerticalAlignment','top','FontSize',8,'BackgroundColor','w');
    end
    yline(ax,0,':k'); grid(ax,'on'); ylabel(ax,'Anomaly (^oC)');
    title(ax,plotNames{ip},'FontWeight','bold');
end
xlabel(tl3,'Year'); title(tl3,'Monthly deseasonalized inter-station thermal contrast and detected breaks');
exportgraphics(fig3,fullfile(FIG_DIR,'Fig3_DT_Breakpoints.jpg'),'Resolution',300);

%% =========================================================================
% SECTION 16B: TWO-PANEL CLIMATOLOGY + ANNUAL TREND COMPARISON (FIG. 4)
% Panel (a): monthly climatology across complete pre-discontinuity years.
% Panel (b): annual mean contrast, comparing pre-break trends (solid) with
%            full-record trends (dashed). Post-break annual observations are
%            shown as open markers and are excluded from the pre-break fit.
% =========================================================================
fig4 = figure('Color','w','Position',[70 45 1080 790]);
tl4 = tiledlayout(fig4,2,1,'Padding','compact','TileSpacing','compact');

% ---------------- Panel (a): monthly climatology --------------------------
ax4a = nexttile(tl4,1); hold(ax4a,'on'); box(ax4a,'on');
climSet = {DTmax,DTmean,DTmin};
climNames = {'\DeltaT_{max} (day)','\DeltaT_{mean}','\DeltaT_{min} (night)'};
climCols = {[0.84 0.37 0.00],[0.25 0.25 0.25],[0.00 0.45 0.70]};
climMean = NaN(12,3); climCI = NaN(12,3); climN = zeros(12,3);
for kk=1:3
    for mm=1:12
        vals=NaN(numel(uYrSubC),1);
        for iy=1:numel(uYrSubC)
            selY=(yr==uYrSubC(iy)) & (mo==mm);
            vals(iy)=mean(climSet{kk}(selY),'omitnan');
        end
        vals=vals(isfinite(vals));
        climN(mm,kk)=numel(vals);
        if ~isempty(vals), climMean(mm,kk)=mean(vals); end
        if numel(vals)>=2
            climCI(mm,kk)=tinv(0.975,numel(vals)-1)*std(vals)/sqrt(numel(vals));
        end
    end
    errorbar(ax4a,1:12,climMean(:,kk),climCI(:,kk),'-o','Color',climCols{kk}, ...
        'MarkerFaceColor',climCols{kk},'LineWidth',1.5,'CapSize',4, ...
        'DisplayName',climNames{kk});
end
% Shade the dry season (May-October), matching the manuscript's seasonal definition.
yl=ylim(ax4a);
patch(ax4a,[4.5 10.5 10.5 4.5],[yl(1) yl(1) yl(2) yl(2)], ...
    [0.94 0.94 0.90],'EdgeColor','none','FaceAlpha',0.55, ...
    'HandleVisibility','off');
uistack(findobj(ax4a,'Type','ErrorBar'),'top');
yline(ax4a,0,':k','HandleVisibility','off');
grid(ax4a,'on'); xlim(ax4a,[0.5 12.5]);
set(ax4a,'XTick',1:12,'XTickLabel',{'Jan','Feb','Mar','Apr','May','Jun', ...
    'Jul','Aug','Sep','Oct','Nov','Dec'},'FontSize',9);
xlabel(ax4a,'Month'); ylabel(ax4a,'\DeltaT (°C)');
title(ax4a,'(a) Monthly climatology of the inter-station thermal contrast');
text(ax4a,7.5,yl(2)-0.04*range(yl),'dry season (May-Oct)', ...
    'HorizontalAlignment','center','VerticalAlignment','top', ...
    'FontSize',8,'FontAngle','italic','Color',[0.35 0.35 0.35]);
legend(ax4a,'Location','northwest','NumColumns',1,'FontSize',8);

% ---------------- Panel (b): annual values and two trend windows ----------
ax4b = nexttile(tl4,2); hold(ax4b,'on'); box(ax4b,'on');
dtAnnualSet = {DTmax_a(:),DTmean_a(:),DTmin_a(:)};
dtAnnualCols = climCols;
dtAnnualMarks = {'o','o','o'};
dtAnnualLabels = climNames;
postYear = uniqYr > max(uYrSubC);
% Mark post-break years as excluded from pre-break fitting.
if any(postYear)
    xPost0=min(uniqYr(postYear))-0.5;
    patch(ax4b,[xPost0 max(uniqYr)+0.5 max(uniqYr)+0.5 xPost0], ...
        [-2 -2 3.2 3.2],[0.94 0.94 0.94], ...
        'EdgeColor','none','FaceAlpha',0.65,'HandleVisibility','off');
end

trendText = strings(1,3);
for kk=1:3
    aa=dtAnnualSet{kk};
    % Annual observations: filled markers for pre-break years, open markers
    % for years excluded from the pre-break fit.
    hPre=plot(ax4b,uniqYr(~postYear),aa(~postYear),dtAnnualMarks{kk}, ...
        'Color',dtAnnualCols{kk},'MarkerFaceColor',dtAnnualCols{kk}, ...
        'MarkerSize',4.5,'LineStyle','none','HandleVisibility','off');
    if any(postYear)
        plot(ax4b,uniqYr(postYear),aa(postYear),dtAnnualMarks{kk}, ...
            'Color',dtAnnualCols{kk},'MarkerFaceColor','w', ...
            'MarkerSize',5,'LineStyle','none','HandleVisibility','off');
    end

    % Pre-break trend: complete calendar years only.
    preIdx=ismember(uniqYr,uYrSubC);
    xPre=uniqYr(preIdx); yPre=aa(preIdx);
    SPre=trendtest(yPre,xPre);
    yPreFit=mean(yPre,'omitnan')+SPre.slope*(xPre-mean(xPre));
    plot(ax4b,xPre,yPreFit,'-','Color',dtAnnualCols{kk}, ...
        'LineWidth',1.8,'HandleVisibility','off');

    % Full-record trend: all available annual values.
    validFull=isfinite(aa)&isfinite(uniqYr);
    xFull=uniqYr(validFull); yFull=aa(validFull);
    SFull=trendtest(yFull,xFull);
    yFullFit=mean(yFull,'omitnan')+SFull.slope*(xFull-mean(xFull));
    plot(ax4b,xFull,yFullFit,'--','Color',dtAnnualCols{kk}, ...
        'LineWidth',1.35,'HandleVisibility','off');

    trendText(kk)=sprintf('%s: pre %+.3f (%s); full %+.3f (%s)', ...
        dtAnnualLabels{kk},SPre.slope,pstr(SPre.p),SFull.slope,pstr(SFull.p));
end
yline(ax4b,0,':k','HandleVisibility','off');
grid(ax4b,'on'); xlim(ax4b,[min(uniqYr)-0.5 max(uniqYr)+0.5]);
ylim(ax4b,[-2 3.2]);
set(ax4b,'FontSize',9);
xlabel(ax4b,'Year'); ylabel(ax4b,'\DeltaT (°C)');
title(ax4b,'(b) Annual mean contrast: pre-break trends (solid) and full-record trends (dashed)');
% Compact slope summary at the top of the panel.
summaryStr=strjoin(cellstr(trendText),newline);
text(ax4b,0.02,0.98,summaryStr,'Units','normalized', ...
    'VerticalAlignment','top','HorizontalAlignment','left', ...
    'FontSize',7.3,'BackgroundColor','w','Margin',3,'Interpreter','tex');
% Explain filled/open annual markers and trend-line styles.
text(ax4b,0.98,0.04,'filled = pre-break; open = excluded post-break; -- full-record', ...
    'Units','normalized','HorizontalAlignment','right', ...
    'FontSize',7.2,'Color',[0.25 0.25 0.25]);

title(tl4,sprintf('Monthly climatology and annual thermal-contrast trends (%d-%d)', ...
    min(uniqYr),max(uniqYr)),'FontWeight','bold');
exportgraphics(fig4,fullfile(FIG_DIR,'Fig4_Monthly_Climatology_and_Annual_Trends.jpg'), ...
    'Resolution',300,'ContentType','image');

T_monthly_climatology=table((1:12)',climMean(:,1),climCI(:,1),climN(:,1), ...
    climMean(:,2),climCI(:,2),climN(:,2),climMean(:,3),climCI(:,3),climN(:,3), ...
    'VariableNames',{'Month','DTmax_mean_degC','DTmax_CI95_halfwidth','DTmax_N_years', ...
    'DTmean_mean_degC','DTmean_CI95_halfwidth','DTmean_N_years', ...
    'DTmin_mean_degC','DTmin_CI95_halfwidth','DTmin_N_years'});

%% =========================================================================
% SECTION 16C: DAILY-SCALE COUPLING (FIG. 5)
% Partial/binned relationships are adjusted for the other predictor and
% calendar-month fixed effects. Table coefficients are from the same lag-0
% multivariable model. Unadjusted rainfall means are shown faintly as context.
% =========================================================================
fig5 = figure('Color','w','Position',[70 50 1120 820]);
tl5 = tiledlayout(2,2,'Padding','compact','TileSpacing','compact');
resp5 = {DTmax(:),DTmin(:)}; respLab5 = {'\DeltaT_{max} (day)','\DeltaT_{min} (night)'};
respCol5 = {[0.84 0.37 0.00],[0.00 0.45 0.70]};
rainEdges=[0 1 5 20 50 Inf]; rainLabs={'<1','1-<5','5-<20','20-<50','>=50'};
binRows=cell(0,8); panelLabels={'(a)','(b)','(c)','(d)'}; pp=0;
for ir=1:2
    y=resp5{ir}; good=isfinite(y)&isfinite(rainH)&isfinite(rhH)&isfinite(mo);
    yy=y(good); rr=rainH(good); hh=rhH(good); mm=mo(good); dnGood=dn(good);
    DD=zeros(numel(yy),11); for mm0=2:12, DD(:,mm0-1)=double(mm==mm0); end
    XX=[ones(numel(yy),1),rr,hh,DD]; bb=XX\yy;
    monthEff=DD*bb(4:end); monthCtr=monthEff-mean(monthEff);
    yAdjRain=yy-bb(3)*(hh-mean(hh))-monthCtr;
    yAdjRH=yy-bb(2)*(rr-mean(rr))-monthCtr;
    % Rainfall categories: adjusted and unadjusted means use identical days.
    pp=pp+1; ax=nexttile(tl5); hold(ax,'on'); box(ax,'on');
    rawM=NaN(5,1); adjM=rawM; adjCI=rawM; nBin=zeros(5,1);
    for j=1:5
        sel=rr>=rainEdges(j)&rr<rainEdges(j+1);
        nBin(j)=sum(sel);
        if nBin(j)>=30
            rawM(j)=mean(yy(sel)); adjM(j)=mean(yAdjRain(sel));
            [~,adjCI(j)]=meanHACCI(yAdjRain(sel),dnGood(sel),HAC_LAG);
        end
    end
    % Make the unadjusted series visually distinct and sufficiently strong for print.
    hRaw = plot(ax,1:5,rawM,'--o','Color',[0.48 0.48 0.48], ...
        'MarkerFaceColor','w','MarkerEdgeColor',[0.35 0.35 0.35], ...
        'LineWidth',1.8,'MarkerSize',6,'DisplayName','Unadjusted');
    hAdj = errorbar(ax,1:5,adjM,adjCI,'-o','Color',respCol5{ir}, ...
        'MarkerFaceColor',respCol5{ir},'MarkerEdgeColor',respCol5{ir}, ...
        'LineWidth',1.8,'MarkerSize',6,'CapSize',6, ...
        'DisplayName','Adjusted for RH and month');
    yline(ax,0,':k','HandleVisibility','off');
    set(ax,'XTick',1:5,'XTickLabel',rainLabs,'FontSize',10, ...
        'LineWidth',0.9,'GridAlpha',0.18);
    xlabel(ax,'Daily rainfall (mm)');
    ylabel(ax,'Adjusted \DeltaT (°C)');
    f0=modelFit.(sprintf('r%d_l0',ir));
    % sprintf treats \b as backspace; escape the slash so MATLAB's TeX interpreter receives \beta.
    ttl = sprintf('%s %s vs rainfall | \\beta=%+.4f, %s', ...
        panelLabels{pp},respLab5{ir},f0.beta(2),pstr(f0.p(2)));
    title(ax,ttl,'Interpreter','tex','FontWeight','bold');
    legend(ax,[hRaw hAdj],'Location','best','Box','off');
    for j=1:5
        binRows(end+1,:)={respLab5{ir},'Rainfall',j,rainLabs{j},adjM(j),adjCI(j),rawM(j),nBin(j)}; %#ok<SAGROW>
    end
    % RH quantile bins, after removing rainfall and month contributions.
    pp=pp+1; ax=nexttile(tl5); hold(ax,'on'); box(ax,'on');
    edgesRH=unique(quantile(hh,linspace(0,1,11)));
    xb=NaN(max(numel(edgesRH)-1,0),1); yb=xb; eb=xb; nb=zeros(size(xb));
    if numel(edgesRH)>=2
        binRH=discretize(hh,edgesRH); binRH(hh==edgesRH(end))=numel(edgesRH)-1;
        for jb=1:numel(xb)
            selRH=(binRH==jb); nb(jb)=sum(selRH);
            if nb(jb)>=30
                xb(jb)=mean(hh(selRH)); yb(jb)=mean(yAdjRH(selRH));
                [~,eb(jb)]=meanHACCI(yAdjRH(selRH),dnGood(selRH),HAC_LAG);
            end
        end
    end
    keep=isfinite(xb)&isfinite(yb)&isfinite(eb);
    if any(keep)
        errorbar(ax,xb(keep),yb(keep),eb(keep),'-o','Color',respCol5{ir}, ...
            'MarkerFaceColor',respCol5{ir},'MarkerEdgeColor',respCol5{ir}, ...
            'LineWidth',1.8,'MarkerSize',6,'CapSize',5);
    end
    f0=modelFit.(sprintf('r%d_l0',ir));
    xm=mean(hh); ym=mean(yAdjRH);
    if any(keep)
        xl=[min(xb(keep)) max(xb(keep))];
        plot(ax,xl,ym+f0.beta(3)*(xl-xm),'--k','LineWidth',1.4,'HandleVisibility','off');
    end
    yline(ax,0,':k','HandleVisibility','off');
    set(ax,'FontSize',10,'LineWidth',0.9,'GridAlpha',0.18);
    xlabel(ax,'Relative humidity (%)'); ylabel(ax,'Adjusted \DeltaT (°C)');
    ttl = sprintf('%s %s vs RH | \\beta=%+.4f, %s', ...
        panelLabels{pp},respLab5{ir},f0.beta(3),pstr(f0.p(3)));
    title(ax,ttl,'Interpreter','tex','FontWeight','bold');
    for j=1:numel(xb)
        binRows(end+1,:)={respLab5{ir},'RH',j,xb(j),yb(j),eb(j),NaN,nb(j)}; %#ok<SAGROW>
    end
end
title(tl5,'Daily-scale modulation of inter-station thermal contrast by rainfall and humidity');
exportgraphics(fig5,fullfile(FIG_DIR,'Fig5_Daily_Coupling.jpg'),'Resolution',300);
T_binned_coupling=cell2table(binRows,'VariableNames', ...
    {'Response','Predictor','Bin','Bin_label_or_midpoint','Adjusted_mean_degC','CI95_degC','Unadjusted_mean_degC','N_days'});

%% =========================================================================
% SECTION 16D: EXTREME-RAINFALL EVENT COMPOSITE (FIG. 2)
% Restored to the master pipeline. Uses station-specific p99 wet-day events,
% separated by at least 10 days, and calendar-month anomalies. Event day is
% lag 0; the antecedent window is visually shaded from -7 to -1 days.
% =========================================================================
LAG_LO_EVT = -10;
LAG_HI_EVT = 3;
MIN_SEP_EVT = 10;
WET_EVT = 1.0;
lagsEvt = LAG_LO_EVT:LAG_HI_EVT;
monthEvt = month(date);
if isdatetime(date)
    dayNumEvt = days(date-date(1));
else
    dayNumEvt = date-date(1);
end
dayNumEvt = round(dayNumEvt(:));
nEvt = height(T);

% Calendar-month anomalies, consistent with the original event-composite method.
anDTmaxEvt = local_anom(DTmax(:),monthEvt(:));
anTmaxH = local_anom(T.Tmax_HND(:),monthEvt(:));
anTmaxM = local_anom(T.Tmax_MKS(:),monthEvt(:));
anHumH = local_anom(T.Hum_HND(:),monthEvt(:));
anHumM = local_anom(T.Hum_MKS(:),monthEvt(:));

stationEvt = struct([]);
stationEvt(1).name = 'HND'; stationEvt(1).rain = T.Rain_HND(:);
stationEvt(1).tmaxAnom = anTmaxH; stationEvt(1).humAnom = anHumH;
stationEvt(1).col = cHND;
stationEvt(2).name = 'MKS'; stationEvt(2).rain = T.Rain_MKS(:);
stationEvt(2).tmaxAnom = anTmaxM; stationEvt(2).humAnom = anHumM;
stationEvt(2).col = cMKS;

fig2 = figure('Name','Fig 2 - Extreme-rainfall event composite', ...
    'NumberTitle','off','Position',[70 50 1050 850],'Color','w');
tl2 = tiledlayout(fig2,3,1,'Padding','compact','TileSpacing','compact');
comp2 = struct([]);
for ss=1:2
    rrEvt=stationEvt(ss).rain;
    wetEvt=rrEvt(rrEvt>=WET_EVT & isfinite(rrEvt));
    if isempty(wetEvt)
        comp2(ss).n=0; comp2(ss).P=struct();
        continue
    end
    p99Evt=prctile(wetEvt,99);
    evEvt=local_pickevents(rrEvt,p99Evt,MIN_SEP_EVT);
    P2.Hum=NaN(numel(evEvt),numel(lagsEvt));
    P2.Tmax=NaN(numel(evEvt),numel(lagsEvt));
    P2.DTmax=NaN(numel(evEvt),numel(lagsEvt));
    for ee=1:numel(evEvt)
        for qq=1:numel(lagsEvt)
            ii=evEvt(ee)+lagsEvt(qq);
            if ii<1 || ii>nEvt, continue; end
            if dayNumEvt(ii)-dayNumEvt(evEvt(ee))~=lagsEvt(qq), continue; end
            if ss==1
                P2.Hum(ee,qq)=anHumH(ii);
                P2.Tmax(ee,qq)=anTmaxH(ii);
            else
                P2.Hum(ee,qq)=anHumM(ii);
                P2.Tmax(ee,qq)=anTmaxM(ii);
            end
            P2.DTmax(ee,qq)=anDTmaxEvt(ii);
        end
    end
    comp2(ss).n=numel(evEvt); comp2(ss).P=P2; comp2(ss).events=evEvt;
    comp2(ss).p99=p99Evt;
end

rowsEvt = {'Hum','RH anomaly (points)','Relative humidity at the affected station'; ...
           'Tmax','T_{max} anomaly (°C)','Daily maximum temperature at the affected station'; ...
           'DTmax','\DeltaT_{max} anomaly (°C)','Daytime inter-station thermal contrast'};
for rr=1:3
    ax2=nexttile(tl2); hold(ax2,'on'); box(ax2,'on');
    hEvt=gobjects(0); labEvt={};
    for ss=1:2
        if comp2(ss).n<3, continue; end
        matEvt=comp2(ss).P.(rowsEvt{rr,1});
        muEvt=mean(matEvt,1,'omitnan');
        nEvtLag=sum(isfinite(matEvt),1);
        seEvt=std(matEvt,0,1,'omitnan')./max(sqrt(nEvtLag),1);
        okEvt=nEvtLag>=3;
        if ~any(okEvt), continue; end
        % 95% CI around the event mean; use t critical values at each lag.
        ciEvt=NaN(size(seEvt));
        for jj=find(okEvt)
            ciEvt(jj)=tinv(0.975,nEvtLag(jj)-1)*seEvt(jj);
        end
        xfEvt=[lagsEvt(okEvt),fliplr(lagsEvt(okEvt))];
        yfEvt=[muEvt(okEvt)+ciEvt(okEvt),fliplr(muEvt(okEvt)-ciEvt(okEvt))];
        fill(ax2,xfEvt,yfEvt,stationEvt(ss).col,'FaceAlpha',0.15, ...
            'EdgeColor','none','HandleVisibility','off');
        hEvt(end+1)=plot(ax2,lagsEvt(okEvt),muEvt(okEvt),'-o', ...
            'Color',stationEvt(ss).col,'LineWidth',1.8, ...
            'MarkerFaceColor',stationEvt(ss).col,'MarkerSize',4.5); %#ok<SAGROW>
        labEvt{end+1}=sprintf('%s, p99 events (n = %d)', ...
            stationEvt(ss).name,comp2(ss).n); %#ok<SAGROW>
    end
    yline(ax2,0,':k','HandleVisibility','off');
    xline(ax2,0,'k-','LineWidth',1,'HandleVisibility','off');
    yl2=ylim(ax2);
    patch(ax2,[ -7.5 -0.5 -0.5 -7.5],[yl2(1) yl2(1) yl2(2) yl2(2)], ...
        [0.4 0.4 0.4],'FaceAlpha',0.07,'EdgeColor','none','HandleVisibility','off');
    % Reorder children so the antecedent-window shading remains behind data.
    ch2=get(ax2,'Children');
    if numel(ch2)>1, uistack(ch2(1),'bottom'); end
    ylim(ax2,yl2); xlim(ax2,[LAG_LO_EVT-0.5 LAG_HI_EVT+0.5]);
    grid(ax2,'on'); set(ax2,'Layer','top','FontSize',10,'GridAlpha',0.18);
    ylabel(ax2,rowsEvt{rr,2});
    title(ax2,sprintf('(%s) %s',char('a'+rr-1),rowsEvt{rr,3}), ...
        'FontWeight','normal','Interpreter','tex');
    if ~isempty(hEvt), legend(ax2,hEvt,labEvt,'Location','best','FontSize',8,'Box','off'); end
    if rr==3, xlabel(ax2,'Days relative to the extreme-rainfall day (0 = event)'); end
end
title(tl2,'Antecedent conditions and daily response around p99 rainfall events');
if SAVE_FIGS
    exportgraphics(fig2,fullfile(FIG_DIR,'Fig2_Extreme_Rainfall_Event_Composite.jpg'), ...
        'Resolution',300,'ContentType','image');
end
fprintf('Fig2 event composite: HND n=%d; MKS n=%d independent p99 events.\n', ...
    comp2(1).n,comp2(2).n);

fprintf('\n=== FIGURES EXPORTED ===\n');
fprintf('Fig1: annual thermal contrast; Fig2: p99 event composite; Fig3: breakpoints; Fig4: monthly climatology + annual trend comparison; Fig5: daily coupling; Fig6: rainfall extremes.\n');
%% =========================================================================
%  SECTION 17: EXPORT NUMERICAL TABLES
%  All key numerical tables are written to XLSX workbooks only.
%  Figures remain exported separately at 300 dpi in Section 16.
% =========================================================================
OUT_DIR = 'Output_Tables';
if ~exist(OUT_DIR,'dir'); mkdir(OUT_DIR); end

% 1) Monthly relative-homogeneity breakpoint results
nBreak = numel(homogResults);
Series = strings(nBreak,1); N_months=NaN(nBreak,1); SNHT_stat=NaN(nBreak,1);
SNHT_p=NaN(nBreak,1); Pettitt_stat=NaN(nBreak,1); Pettitt_p=NaN(nBreak,1);
Last_pre_YYYYMM=NaN(nBreak,1); First_post_YYYYMM=NaN(nBreak,1);
Pre_mean=NaN(nBreak,1); Post_mean=NaN(nBreak,1); Shift=NaN(nBreak,1);
for ii=1:nBreak
    q=homogResults(ii); Series(ii)=string(q.name); N_months(ii)=numel(q.monthly);
    SNHT_stat(ii)=q.snhtStatistic; SNHT_p(ii)=q.snhtP;
    Pettitt_stat(ii)=q.pettittStatistic; Pettitt_p(ii)=q.pettittP;
    Last_pre_YYYYMM(ii)=q.lastPreYM; First_post_YYYYMM(ii)=q.firstPostYM;
    Pre_mean(ii)=q.rawPreMean; Post_mean(ii)=q.rawPostMean; Shift(ii)=q.rawShift;
end
T_breakpoints=table(Series,N_months,SNHT_stat,SNHT_p,Pettitt_stat,Pettitt_p, ...
    Last_pre_YYYYMM,First_post_YYYYMM,Pre_mean,Post_mean,Shift);

% 2) Annual thermal-contrast time series
T_thermal_annual=table(uniqYr(:),'VariableNames',{'Year'});
for ii=1:3
    xx=dtSeries{ii}; av=NaN(numel(uniqYr),1);
    for jj=1:numel(uniqYr); av(jj)=mean(xx(yr==uniqYr(jj)),'omitnan'); end
    T_thermal_annual.(matlab.lang.makeValidName(dtNames{ii}))=av;
end

% 3) Annual ETCCDI values, long-format table
T_etccdi_annual = table();
for ss=1:2
    if ss==1; stName="MKS"; E=etccdi_MKS; else; stName="HND"; E=etccdi_HND; end
    nn=numel(E.Year); tmp=table(repmat(stName,nn,1),E.Year(:), ...
        E.RX1day(:),E.RX5day(:),E.SDII(:),E.R95p(:),E.R99p(:),E.CDD(:),E.PRCPTOT(:), ...
        E.nValidDays(:),E.nExpectedDays(:),E.isCompleteYear(:), ...
        'VariableNames',{'Station','Year','RX1day','RX5day','SDII','R95p','R99p','CDD','PRCPTOT', ...
        'N_valid_days','N_expected_days','Complete_year'});
    T_etccdi_annual=[T_etccdi_annual;tmp]; %#ok<AGROW>
end

% 4) Full-record and pre-discontinuity ETCCDI trend comparison
TrendRows={}; row=0;
for ss=1:2
    if ss==1; stName="MKS"; Ef=etccdi_MKS; Eh=etSubMKS; else; stName="HND"; Ef=etccdi_HND; Eh=etSubHND; end
    for jj=1:numel(idxNames)
        nm=idxNames{jj};
        Sf=trendtest(Ef.(nm),Ef.Year); MKf=mann_kendall(Ef.(nm));
        Sh=trendtest(Eh.(nm),Eh.Year); MKh=mann_kendall(Eh.(nm));
        row=row+1;
        TrendRows(row,:)={stName,string(nm),Sf.slope,Sf.ci(1),Sf.ci(2),Sf.senSlope,Sf.senCI(1),Sf.senCI(2),Sf.R2,Sf.p,MKf.p, ...
            Sh.slope,Sh.ci(1),Sh.ci(2),Sh.senSlope,Sh.senCI(1),Sh.senCI(2),Sh.R2,Sh.p,MKh.p};
    end
end
T_etccdi_trends=cell2table(TrendRows,'VariableNames', ...
    {'Station','Index','Full_OLS_slope','Full_OLS_CI_low','Full_OLS_CI_high','Full_Sen_slope','Full_Sen_CI_low','Full_Sen_CI_high','Full_R2','Full_OLS_p','Full_MK_p', ...
     'PreBreak_OLS_slope','PreBreak_OLS_CI_low','PreBreak_OLS_CI_high','PreBreak_Sen_slope','PreBreak_Sen_CI_low','PreBreak_Sen_CI_high','PreBreak_R2','PreBreak_OLS_p','PreBreak_MK_p'});

% 5) Annual thermal-contrast trends for complete pre-discontinuity years
Metric=strings(3,1); OLS_slope=NaN(3,1); OLS_CI_low=NaN(3,1); OLS_CI_high=NaN(3,1);
Sen_slope=NaN(3,1); Sen_CI_low=NaN(3,1); Sen_CI_high=NaN(3,1);
R2=NaN(3,1); OLS_p=NaN(3,1); MK_p=NaN(3,1); N_years=NaN(3,1);
for ii=1:3
    xx=dtAll{ii}; xs=xx(sub); aa=NaN(numel(uYrSubC),1);
    for jj=1:numel(uYrSubC); aa(jj)=mean(xs(yrSub==uYrSubC(jj)),'omitnan'); end
    S=trendtest(aa,uYrSubC); MK=mann_kendall(aa);
    Metric(ii)=string(dtLab{ii}); OLS_slope(ii)=S.slope; OLS_CI_low(ii)=S.ci(1); OLS_CI_high(ii)=S.ci(2);
    Sen_slope(ii)=S.senSlope; Sen_CI_low(ii)=S.senCI(1); Sen_CI_high(ii)=S.senCI(2);
    R2(ii)=S.R2; OLS_p(ii)=S.p; MK_p(ii)=MK.p; N_years(ii)=S.n;
end
T_thermal_trends=table(Metric,OLS_slope,OLS_CI_low,OLS_CI_high,Sen_slope,Sen_CI_low,Sen_CI_high,R2,OLS_p,MK_p,N_years);

% 6) Save numerical tables in one multi-sheet Excel workbook
xlsxFile=fullfile(OUT_DIR,'Main_Analysis_Tables.xlsx');
writetable(T_breakpoints,xlsxFile,'Sheet','Breakpoints');
writetable(T_thermal_annual,xlsxFile,'Sheet','Annual_Thermal');
writetable(T_thermal_trends,xlsxFile,'Sheet','Thermal_Trends');
writetable(T_etccdi_annual,xlsxFile,'Sheet','Annual_ETCCDI');
writetable(T_etccdi_trends,xlsxFile,'Sheet','ETCCDI_Trends');
fprintf('Numerical tables exported to %s (XLSX workbooks).\n',OUT_DIR);

% Additional core result tables: daily multivariable models and seasonal sensitivity.
writetable(T_model_daily,xlsxFile,'Sheet','Daily_Rain_RH_Models');
writetable(T_model_seasonal,xlsxFile,'Sheet','Seasonal_Rain_RH');
writetable(T_binned_coupling,xlsxFile,'Sheet','Figure5_Binned_Values');
writetable(T_monthly_climatology,xlsxFile,'Sheet','Monthly_Climatology');

% QC completeness summary from the retained, screened station records.
qRows=cell(0,4);
for ss=1:2
    if ss==1, qTbl=mks_d; st='MKS'; else, qTbl=hnd_d; st='HND'; end
    qVars={'Tmean','Tmin','Tmax','Rain','Humidity'};
    for vv=1:numel(qVars)
        xx=qTbl.(qVars{vv});
        qRows(end+1,:)={st,qVars{vv},numel(xx),sum(isfinite(xx))}; %#ok<SAGROW>
    end
end
T_QC=cell2table(qRows,'VariableNames',{'Station','Variable','N_records','N_valid_after_QC'});
T_QC.Completeness_percent=100*T_QC.N_valid_after_QC./T_QC.N_records;
writetable(T_QC,xlsxFile,'Sheet','QC_Completeness');
% Dates flagged for repeated non-zero rainfall runs; values are retained unless
% REMOVE_FLATLINE_RAINFALL was explicitly enabled above.
flatRows=cell(0,4);
for ss=1:2
    if ss==1, qTbl=mks_d; st='MKS'; else, qTbl=hnd_d; st='HND'; end
    ff=rainFlatFlags{ss};
    ix=find(ff);
    for jj=1:numel(ix)
        flatRows(end+1,:)={st,qTbl.Date(ix(jj)),rainFlatValues{ss}(jj),REMOVE_FLATLINE_RAINFALL}; %#ok<SAGROW>
    end
end
if isempty(flatRows)
    T_Flatline=table(string.empty(0,1),NaT(0,1),zeros(0,1),false(0,1), ...
        'VariableNames',{'Station','Date','Rain_mm','Removed'});
else
    T_Flatline=cell2table(flatRows,'VariableNames',{'Station','Date','Rain_mm','Removed'});
end
writetable(T_Flatline,xlsxFile,'Sheet','Rain_Flatline_Flags');

% Descriptive statistics for the joined analysis records.
descNames={'Tmean_MKS','Tmin_MKS','Tmax_MKS','Rain_MKS','Hum_MKS', ...
           'Tmean_HND','Tmin_HND','Tmax_HND','Rain_HND','Hum_HND'};
Drows=cell(numel(descNames),7);
for ii=1:numel(descNames)
    xx=T.(descNames{ii}); xx=xx(isfinite(xx));
    Drows(ii,:)={descNames{ii},numel(xx),mean(xx),median(xx),std(xx),skewness(xx),kurtosis(xx)};
end
T_Descriptive=cell2table(Drows,'VariableNames', ...
    {'Variable','N','Mean','Median','SD','Skewness','Kurtosis'});
writetable(T_Descriptive,xlsxFile,'Sheet','Descriptive');

% Paired station comparisons (MKS-HND) with matched daily observations.
PairRows=cell(0,5);
pairsOut={'Tmean',T.Tmean_HND,T.Tmean_MKS;'Tmin',T.Tmin_HND,T.Tmin_MKS; ...
          'Tmax',T.Tmax_HND,T.Tmax_MKS;'Rain',T.Rain_HND,T.Rain_MKS; ...
          'RH',T.Hum_HND,T.Hum_MKS};
for ii=1:size(pairsOut,1)
    aa=pairsOut{ii,2}; bb=pairsOut{ii,3}; ok=isfinite(aa)&isfinite(bb); aa=aa(ok); bb=bb(ok);
    [~,pt]=ttest(aa,bb); pw=signrank(aa,bb);
    PairRows(end+1,:)={pairsOut{ii,1},numel(aa),mean(bb-aa),pt,pw}; %#ok<SAGROW>
end
T_Paired=cell2table(PairRows,'VariableNames',{'Variable','N_pairs','Mean_MKS_minus_HND','Paired_t_p','Wilcoxon_p'});
writetable(T_Paired,xlsxFile,'Sheet','Paired_Comparison');

% Save exact binned data behind Figure 5 as a separate flat file as well.
writetable(T_binned_coupling,fullfile(OUT_DIR,'Figure5_Binned_Values.csv'));
writetable(T_model_daily,fullfile(OUT_DIR,'Table5_Daily_Rain_RH.csv'));
fprintf('All numerical tables saved in %s\\n',xlsxFile);
fprintf('Figures saved in %s at 300 dpi.\\n',FIG_DIR);
fprintf('=== ANALYSIS COMPLETE ===\\n');

%% ========================================================================
% LOCAL FUNCTIONS - ALL FUNCTIONS KEPT AT THE END
% ========================================================================
function s = fmtP(p)
% Display small/underflowed p-values without ever reporting p = 0.
    if ~isfinite(p)
        s = 'NA';
    elseif p <= 0.001
        s = '<=0.001';
    else
        s = sprintf('%.3g',p);
    end
end


function s = pstr(p)
    if p < 0.001
        s = 'p <= 0.001';
    else
        s = sprintf('p = %.3f', p);
    end
end


function printTrendRow(label, meanVal, S, MK)
    fprintf('%-16s %8.3f %9.4f [%+.4f,%+.4f] %9.4f [%+.4f,%+.4f] %7.3f %7.3f %7.4f %6d\n', ...
        label, meanVal, S.slope, S.ci(1), S.ci(2), S.senSlope, ...
        S.senCI(1), S.senCI(2), S.R2, MK.Z, MK.p, S.n);
end


function [betaRH, pRH] = lag0HAC_rh(DT, Hum, Dm)
    good = isfinite(DT) & isfinite(Hum) & all(isfinite(Dm),2);
    U = DT(good); H = Hum(good); D = Dm(good,:);
    X = [ones(size(U)) H D];
    [b, ~, ~, ~, resid] = ols(X, U);
    se = neweyWest(X, resid, 5);
    n = numel(U); k = size(X,2);
    tH = b(2)/se(2);
    betaRH = b(2);
    pRH = 2*(1 - tcdf(abs(tH), n-k));
end


function T = fixNumericCols(T, cols)
    for i = 1:numel(cols)
        c = cols{i};
        if ismember(c, T.Properties.VariableNames)
            v = T.(c);
            if iscell(v) || isstring(v)
                v = strrep(string(v), ',', '.');
                T.(c) = str2double(v);
            end
        end
    end
end


function T = renameDaily(T)
    T.Properties.VariableNames{'Date'} = 'Date';
    T.Properties.VariableNames{'Tmean'}   = 'Tmean';
    T.Properties.VariableNames{'Rain'}      = 'Rain';
    T.Properties.VariableNames{'Hummidity'}      = 'Humidity';
end


function S = trendtest(y, x)
    % OLS linear trend and Theil-Sen slope with an approximate rank-based CI.
    % OLS CI and Sen CI are deliberately reported separately.
    x = x(:); y = y(:);
    good = isfinite(x) & isfinite(y);
    x = x(good); y = y(good);
    n = length(y); df = n - 2;
    S.slope=NaN; S.se=NaN; S.ci=[NaN NaN]; S.p=NaN;
    S.R2=NaN; S.n=n; S.senSlope=NaN; S.senCI=[NaN NaN];
    if n < 3 || numel(unique(x)) < 2; return; end

    X = [ones(size(x)) x];
    [b,~,~,~,stats] = regress(y, X);
    se = sqrt((sum((y - X*b).^2) / df) / sum((x - mean(x)).^2));
    tval = b(2) / se;
    tcrit = tinv(0.975, df);
    S.slope = b(2);
    S.se    = se;
    S.ci    = [b(2) - tcrit*se, b(2) + tcrit*se];
    S.p     = 2*(1 - tcdf(abs(tval), df));
    S.R2    = stats(1);

    % Theil-Sen slope: median of all pairwise slopes using actual x spacing.
    N = n*(n-1)/2;
    slopes = NaN(N,1); k=0;
    for ii=1:n-1
        jj=(ii+1):n;
        m=numel(jj);
        slopes(k+1:k+m)=(y(jj)-y(ii))./(x(jj)-x(ii));
        k=k+m;
    end
    slopes=sort(slopes(isfinite(slopes)));
    if isempty(slopes); return; end
    S.senSlope=median(slopes);

    % Approximate 95% CI for Sen slope from the tie-corrected MK variance.
    % This is the conventional rank-based order-statistic interval.
    [~,~,g]=unique(y);
    counts=accumarray(g,1);
    tieTerm=sum(counts.*(counts-1).*(2*counts+5));
    varS=(n*(n-1)*(2*n+5)-tieTerm)/18;
    NN=numel(slopes);
    if isfinite(varS) && varS>0 && NN>=2
        C=norminv(0.975)*sqrt(varS);
        lo=max(1,floor((NN-C)/2));
        hi=min(NN,ceil((NN+C)/2)+1);
        if lo<=hi
            S.senCI=[slopes(lo), slopes(hi)];
        end
    end
end


function MK = mann_kendall(y)
    % Standard Mann-Kendall test with tie-corrected variance.
    % The S statistic is unchanged by ties; Var(S) is reduced according
    % to the sizes of tied groups, as required for repeated annual indices.
    y = y(:); y = y(isfinite(y));
    n = numel(y);
    MK.n=n; MK.S=NaN; MK.Z=NaN; MK.p=NaN; MK.varS=NaN;
    if n < 2; return; end

    S=0;
    for i=1:n-1
        S=S+sum(sign(y(i+1:end)-y(i)));
    end
    [~,~,g]=unique(y);
    counts=accumarray(g,1);
    tieTerm=sum(counts.*(counts-1).*(2*counts+5));
    varS=(n*(n-1)*(2*n+5)-tieTerm)/18;
    if ~isfinite(varS) || varS<=0; return; end

    if S>0; Z=(S-1)/sqrt(varS);
    elseif S<0; Z=(S+1)/sqrt(varS);
    else; Z=0;
    end
    MK.S=S; MK.Z=Z; MK.p=2*(1-normcdf(abs(Z))); MK.varS=varS;
end


function p = mk_permutation_p(y, B)
    % Two-sided Monte Carlo permutation MK p-value under exchangeability.
    % Ties are retained because the observed values are randomly reordered.
    y=y(:); y=y(isfinite(y)); n=numel(y);
    p=NaN;
    if n<3 || B<1; return; end
    Sobs=mk_S_stat(y);
    exceed=0;
    for b=1:B
        yp=y(randperm(n));
        if abs(mk_S_stat(yp)) >= abs(Sobs)
            exceed=exceed+1;
        end
    end
    p=(exceed+1)/(B+1);
end


function S = mk_S_stat(y)
    n=numel(y); S=0;
    for ii=1:n-1
        S=S+sum(sign(y(ii+1:end)-y(ii)));
    end
end


function q = bh_adjust(p)
    % Benjamini-Hochberg adjusted p-values, preserving original order.
    p=p(:); q=NaN(size(p));
    good=find(isfinite(p));
    m=numel(good);
    if m==0; return; end
    [ps,ord]=sort(p(good));
    adj=ps.*m./(1:m)';
    adj=min(1,flipud(cummin(flipud(adj))));
    q(good(ord))=adj;
end


function MK = mann_kendall_hamed_rao_v8(y)
    % Hamed & Rao (1998) variance correction, used here as a sensitivity
    % diagnostic for short annual series. Both independent and corrected
    % variances use the standard tie correction.
    %
    % Rank ACF is calculated as a conventional sample autocorrelation:
    % numerator uses the full-series rank mean and denominator uses the
    % full-series sum of squares. This avoids correlating two separately
    % centered lagged segments, which is not the usual ACF estimator.
    y = y(:); y = y(isfinite(y));
    n = numel(y);
    MK.n=n; MK.S=NaN; MK.Z=NaN; MK.p=NaN; MK.varS=NaN;
    MK.varRatio=NaN; MK.rankACF=NaN(1,3); MK.sigLags=[];
    MK.acfBound=NaN; MK.extremeFactor=false;
    MK.reason='insufficient sample'; MK.senSlope=NaN;
    if n < 8; return; end

    % S and tie-corrected independent variance.
    S=0;
    for i=1:n-1
        S=S+sum(sign(y(i+1:end)-y(i)));
    end
    [~,~,g]=unique(y);
    counts=accumarray(g,1);
    tieTerm=sum(counts.*(counts-1).*(2*counts+5));
    varS=(n*(n-1)*(2*n+5)-tieTerm)/18;
    MK.S=S; MK.varS=varS;
    if ~isfinite(varS) || varS<=0
        MK.reason='zero or invalid tie-corrected Var(S)'; return;
    end

    % Sen-slope detrending (intercept is immaterial for ranks).
    slopes=zeros(n*(n-1)/2,1); k=0;
    for i=1:n-1
        m=n-i;
        slopes(k+1:k+m)=(y(i+1:end)-y(i))./(1:m)';
        k=k+m;
    end
    senSlope=median(slopes);
    MK.senSlope=senSlope;
    detrended=y-senSlope*(1:n)';
    rr=tiedrank(detrended);
    rr=rr-mean(rr);
    denom=sum(rr.^2);
    if ~isfinite(denom) || denom<=0
        MK.reason='zero variance in detrended ranks'; return;
    end

    % Conventional rank autocorrelation using the full-series mean/variance.
    maxLag=min(3,n-1);
    bound=1.96/sqrt(n);
    MK.acfBound=bound;
    rhoAll=NaN(1,maxLag); sigLags=[];
    for lag=1:maxLag
        rhoAll(lag)=sum(rr(1:n-lag).*rr(1+lag:n))/denom;
        if isfinite(rhoAll(lag)) && abs(rhoAll(lag))>bound
            sigLags(end+1)=lag; %#ok<AGROW>
        end
    end
    MK.rankACF=rhoAll;
    MK.sigLags=sigLags;

    % Hamed-Rao n/n* variance multiplier using significant rank-ACF lags.
    corrFactor=1;
    for lag=sigLags
        corrFactor=corrFactor + ...
            2*(n-lag)*(n-lag-1)*(n-lag-2)*rhoAll(lag)/(n*(n-1)*(n-2));
    end
    MK.varRatio=corrFactor;
    if ~isfinite(corrFactor) || corrFactor<=0
        MK.reason='non-positive variance correction factor'; return;
    end
    varSstar=varS*corrFactor;
    if ~isfinite(varSstar) || varSstar<=0
        MK.reason='non-positive corrected variance'; return;
    end
    if S>0; Z=(S-1)/sqrt(varSstar);
    elseif S<0; Z=(S+1)/sqrt(varSstar);
    else; Z=0;
    end
    MK.Z=Z; MK.p=2*(1-normcdf(abs(Z)));
    MK.extremeFactor=(corrFactor<0.5 || corrFactor>2);
    MK.reason='computed';
end


function A = detrended_lag1(y, x)
    % Lag-1 autocorrelation of residuals after removing a linear trend.
    % The p-value is an approximate two-sided t test with df = n-3.
    y = y(:); x = x(:);
    good = isfinite(y) & isfinite(x);
    y = y(good); x = x(good);
    A.n = numel(y);
    A.r1 = NaN; A.p = NaN;
    if numel(y) < 5
        return;
    end
    X = [ones(size(x)) x];
    b = X \ y;
    e = y - X*b;
    e1 = e(1:end-1);
    e2 = e(2:end);
    if std(e1) == 0 || std(e2) == 0
        A.r1 = 0;
        A.p = 1;
        return;
    end
    R = corrcoef(e1, e2);
    r = max(min(R(1,2), 0.999999), -0.999999);
    A.r1 = r;
    df = numel(y) - 3;
    tval = r * sqrt(df / max(1 - r^2, eps));
    A.p = 2 * (1 - tcdf(abs(tval), df));
end


function [b, se, p, R2, resid] = ols(X, y)
    good = all(isfinite(X),2) & isfinite(y);
    X = X(good,:); y = y(good);
    b = X \ y; yhat = X*b; resid = y - yhat;
    n = size(X,1); k = size(X,2);
    s2 = sum(resid.^2) / (n - k);
    covb = s2 * inv(X.'*X);
    se = sqrt(diag(covb));
    t = b ./ se;
    p = 2*(1 - tcdf(abs(t), n-k));
    R2 = 1 - sum(resid.^2)/sum((y - mean(y)).^2);
end


function tbl = run_dt_rh_model(DT, Hum, Dm, label)
    fprintf('\n%s\n', label);
    tbl = [];
    for L = 0:3
        U = DT(L+1:end); H = Hum(1:end-L); D = Dm(L+1:end,:);
        good = isfinite(U) & isfinite(H) & all(isfinite(D),2);
        U = U(good); H = H(good); D = D(good,:);
        X = [ones(size(U)) H D];
        [b, ~, p, R2, resid] = ols(X, U);
        Xbp = [ones(size(H)) H]; bpb = Xbp \ (resid.^2); e2hat = Xbp * bpb;
        e2 = resid.^2;
        R2bp = sum((e2hat - mean(e2)).^2) / sum((e2 - mean(e2)).^2);
        pbp = 1 - chi2cdf(numel(H) * R2bp, 1);
        seHAC = neweyWest(X, resid, 5);
        tH = b(2)/seHAC(2);
        n = numel(U); k = size(X,2);
        pH_hac = 2*(1 - tcdf(abs(tH), n-k));
        tcrit = tinv(0.975, n-k);
        ciH = [b(2) - tcrit*seHAC(2), b(2) + tcrit*seHAC(2)];
        fprintf(['Lag %d | bRH=%+7.4f (p=%.4f, HAC p=%.4f, 95%%CI [%+.4f,%+.4f]) | ' ...
                 'R2=%.3f | BP p=%.4f\n'], ...
            L, b(2), p(2), pH_hac, ciH(1), ciH(2), R2, pbp);
        tbl = [tbl; L, b(2), p(2), pH_hac, ciH(1), ciH(2), R2, pbp]; %#ok<AGROW>
    end
end


function se = neweyWest(X, resid, L)
    n = size(X,1); k = size(X,2);
    XtX_inv = inv(X.'*X);
    S = (X .* resid).' * (X .* resid);
    for l = 1:L
        w = 1 - l/(L+1);
        Xl = X(1:end-l,:) .* resid(1:end-l);
        Xr = X(l+1:end,:) .* resid(l+1:end);
        G = Xl.' * Xr;
        S = S + w * (G + G.');
    end
    V = XtX_inv * S * XtX_inv * (n/(n-k));
    se = sqrt(diag(V));
end


function print_spearman_rh(DT, Hum, label)
    fprintf('- %s -\n', label);
    for L = 0:3
        U = DT(L+1:end); H = Hum(1:end-L);
        good = isfinite(U) & isfinite(H);
        U = U(good); H = H(good);
        [rhoH, pH] = corr(U, H, 'type','Spearman');
        fprintf('  Lag %d: rho(RH)=%+5.3f (p=%.3g)\n', L, rhoH, pH);
    end
end


function [Tmax, kBreak, p] = snht_test(x, nMC)
    x = x(:); x = x(isfinite(x));
    [Tmax, kBreak] = snhtStat(x);
    n = numel(x);
    count = 0;
    for b = 1:nMC
        Tb = snhtStat(randn(n,1));
        if Tb >= Tmax, count = count + 1; end
    end
    p = (count + 1) / (nMC + 1);
end


function [Tmax, kBreak] = snhtStat(x)
    n = numel(x);
    z = (x - mean(x)) / std(x);
    Tk = zeros(n-1,1);
    cs = cumsum(z);
    total = cs(end);
    for k = 1:n-1
        z1 = cs(k)/k;
        z2 = (total - cs(k))/(n-k);
        Tk(k) = k*z1^2 + (n-k)*z2^2;
    end
    [Tmax, kBreak] = max(Tk);
end


function [Umax, kBreak, p] = pettitt_test(x)
    x = x(:); x = x(isfinite(x));
    n = numel(x);
    r = tiedrank(x);
    U = zeros(n,1);
    for k = 1:n
        U(k) = 2*sum(r(1:k)) - k*(n+1);
    end
    [Umax, kBreak] = max(abs(U));
    p = 2 * exp(-6*Umax^2 / (n^3 + n^2));
    p = min(p, 1);
end


function out = computeETCCDI(dates, rain, wetThresh, baseYears, fixedThresholds)
    % ETCCDI indices from a daily series with explicit calendar-day handling.
    % A calendar year is retained only when every date is present exactly once
    % and every daily rainfall value is finite. Missing days are never treated
    % as dry days and are never bridged when constructing 5-day accumulations.
    %
    % fixedThresholds (optional) is a struct with fields p95 and p99. Supplying
    % it keeps the full-record percentile thresholds unchanged in sensitivity
    % analyses that use a shorter analysis window.

    dates = dates(:);
    rain  = rain(:);
    if numel(dates) ~= numel(rain)
        error('computeETCCDI:InputSizeMismatch', ...
              'dates and rain must have the same number of elements.');
    end
    if ~isdatetime(dates)
        error('computeETCCDI:InvalidDates','dates must be a datetime vector.');
    end
    dates = dateshift(dates,'start','day');

    % Sort observations chronologically and reject duplicate calendar dates.
    [dates, order] = sort(dates);
    rain = rain(order);
    if numel(unique(dates)) ~= numel(dates)
        error('computeETCCDI:DuplicateDates', ...
              'Duplicate dates found. Resolve duplicates before ETCCDI calculation.');
    end

    yr = year(dates);
    uniqYr = unique(yr);
    nY = numel(uniqYr);

    RX1day  = NaN(nY,1); RX5day  = NaN(nY,1); SDII = NaN(nY,1);
    CDD     = NaN(nY,1); PRCPTOT = NaN(nY,1);
    R95p    = NaN(nY,1); R99p    = NaN(nY,1);
    nR95day = NaN(nY,1); nR99day = NaN(nY,1);
    nValid  = zeros(nY,1); nExpected = zeros(nY,1);
    isComplete = false(nY,1);

    % Construct complete daily calendar. This makes gaps explicit, so a
    % 5-row moving window cannot silently span a missing calendar date.
    fullDates = (dates(1):caldays(1):dates(end))';
    [hasDate, loc] = ismember(dates, fullDates);
    if ~all(hasDate)
        error('computeETCCDI:CalendarMapping','Could not map observations to daily calendar.');
    end
    fullRain = NaN(size(fullDates));
    fullRain(loc) = rain;

    % Identify strictly complete calendar years.
    for i = 1:nY
        y = uniqYr(i);
        expectedDates = (datetime(y,1,1):caldays(1):datetime(y,12,31))';
        nExpected(i) = numel(expectedDates);
        [present, ix] = ismember(expectedDates, dates);
        validForYear = false(size(present));
        if any(present)
            validForYear(present) = isfinite(rain(ix(present)));
        end
        nValid(i) = sum(validForYear);
        isComplete(i) = all(present) && all(validForYear);
    end

    % Determine station-specific percentile thresholds from complete
    % baseline years only, unless fixed thresholds are supplied.
    if nargin >= 5 && ~isempty(fixedThresholds)
        p95 = fixedThresholds.p95;
        p99 = fixedThresholds.p99;
    else
        if ~isempty(baseYears)
            baseYearMask = uniqYr >= baseYears(1) & uniqYr <= baseYears(2);
        else
            baseYearMask = true(size(uniqYr));
        end
        baseYearsUse = uniqYr(baseYearMask & isComplete);
        baseMask = ismember(yr,baseYearsUse) & isfinite(rain) & rain >= wetThresh;
        wetBase = rain(baseMask);
        if isempty(wetBase)
            error('computeETCCDI:NoBaselineWetDays', ...
                  'No wet days in complete baseline years; percentile thresholds cannot be calculated.');
        end
        p95 = prctile(wetBase,95);
        p99 = prctile(wetBase,99);
    end

    % Five-day running totals on the actual daily calendar. RX5day is
    % assigned to the END date of each consecutive five-day window. A
    % window is valid only when all five calendar days have finite rainfall.
    rain5 = NaN(size(fullRain));
    for j = 5:numel(fullRain)
        w = fullRain(j-4:j);
        if all(isfinite(w))
            rain5(j) = sum(w);
        end
    end

    for i = 1:nY
        y = uniqYr(i);
        if ~isComplete(i)
            continue; % retain NaN for every annual index in an incomplete year
        end

        idx = (yr == y);
        r = rain(idx);
        wet = r(r >= wetThresh);

        RX1day(i)  = max(r);
        if isempty(wet)
            SDII(i) = NaN;
        else
            SDII(i) = mean(wet);
        end
        PRCPTOT(i) = sum(wet);
        R95p(i)    = sum(wet(wet > p95));
        R99p(i)    = sum(wet(wet > p99));
        nR95day(i) = sum(wet > p95);
        nR99day(i) = sum(wet > p99);

        % Missing values break dry spells; complete-year screening above
        % ensures there are no missing values within the annual sequence.
        CDD(i) = maxConsecutive(r < wetThresh);

        % Assign each valid trailing five-day window to its end year.
        endYear = year(fullDates);
        idx5 = (endYear == y) & isfinite(rain5);
        if any(idx5)
            RX5day(i) = max(rain5(idx5));
        end
    end

    out.Year = uniqYr;
    out.RX1day = RX1day; out.RX5day = RX5day;
    out.SDII = SDII; out.R95p = R95p; out.R99p = R99p;
    out.CDD = CDD; out.PRCPTOT = PRCPTOT;
    out.p95 = p95; out.p99 = p99;
    out.nR95day = nR95day; out.nR99day = nR99day;
    out.nValidDays = nValid; out.nExpectedDays = nExpected;
    out.isCompleteYear = isComplete;
    out.incompleteYears = uniqYr(~isComplete);

    if any(~isComplete)
        fprintf('  ETCCDI: excluded incomplete year(s): %s\n', ...
                mat2str(uniqYr(~isComplete)'));
    end
end


function n = maxConsecutive(boolVec)
    boolVec = boolVec(:);
    if isempty(boolVec); n = 0; return; end
    d = diff([0; boolVec; 0]);
    starts = find(d == 1); ends = find(d == -1) - 1;
    if isempty(starts); n = 0; else; n = max(ends - starts + 1); end
end


function printETCCDIrow(label, s)
    fprintf('%-8s %9.2f %9.2f %8.2f %9.2f %9.2f %8.2f %10.2f\n', label, ...
        mean(s.RX1day,'omitnan'), mean(s.RX5day,'omitnan'), mean(s.SDII,'omitnan'), ...
        mean(s.R95p,'omitnan'), mean(s.R99p,'omitnan'), mean(s.CDD,'omitnan'), ...
        mean(s.PRCPTOT,'omitnan'));
end


function y = local_anom(x, mo)
% Departure from the calendar-month mean over the whole record.
    y = NaN(size(x));
    for m = 1:12
        k = (mo == m);
        y(k) = x(k) - mean(x(k), 'omitnan');
    end
end


function w = local_window(i0, ANTE, dnum, n)
% Row indices for the antecedent window, returned ONLY for rows that really
% are the right number of days away. A record with a missing date would
% otherwise shift the window silently. Written as a plain loop: the vectorised
% form has to reconcile row and column orientations, and gets that wrong in
% ways that do not announce themselves.
    w = [];
    for k = 1:numel(ANTE)
        i = i0 + ANTE(k);
        if i < 1 || i > n,                  continue; end
        if dnum(i) - dnum(i0) ~= ANTE(k),   continue; end
        w(end+1) = i; %#ok<AGROW>
    end
end


function ev = local_pickevents(r, thr, minSep)
% Days at or above thr, thinned so no two are within minSep days. When
% several qualify inside one window the largest is kept, so an event is
% represented by its peak day rather than by whichever day came first.
    cand = find(r >= thr & isfinite(r));
    if isempty(cand), ev = []; return; end
    [~, ord] = sort(r(cand), 'descend');
    cand = cand(ord);
    ev = [];
    for k = 1:numel(cand)
        if isempty(ev) || all(abs(cand(k) - ev) >= minSep)
            ev(end+1) = cand(k); %#ok<AGROW>
        end
    end
    ev = sort(ev);
end


function s = local_ci(v)
% 95% CI of the mean, printed as a bracket.
    v = v(isfinite(v));
    if numel(v) < 2, s = '(n < 2)'; return; end
    m  = mean(v);
    se = std(v) / sqrt(numel(v));
    tc = 1.96;
    try
        tc = tinv(0.975, numel(v) - 1);
    catch
        % Statistics Toolbox unavailable; the normal approximation is used.
    end
    s = sprintf('[%+.3f, %+.3f]', m - tc*se, m + tc*se);
end


function pr = local_pctrank(anomSeries, mo, ANTE, i0, dnum, n)
% Percentile of the event's antecedent mean among all days of the same
% calendar month, using an identical antecedent definition for every day.
    m0   = mo(i0);
    vals = NaN(n,1);
    for i = (1 - min(ANTE)) : n
        if mo(i) ~= m0, continue; end
        w = local_window(i, ANTE, dnum, n);
        if numel(w) == numel(ANTE)
            vals(i) = mean(anomSeries(w), 'omitnan');
        end
    end
    vals = vals(isfinite(vals));
    w0 = local_window(i0, ANTE, dnum, n);
    x0 = mean(anomSeries(w0), 'omitnan');
    if isempty(vals) || ~isfinite(x0), pr = NaN; return; end
    pr = 100 * sum(vals <= x0) / numel(vals);
end


function s = local_pstr(p)
% Never prints "p = 0.000".
    if p < 0.001, s = 'p <= 0.001'; else, s = sprintf('p = %.3f', p); end
end


function w = local_word(v)
    if v > 0, w = 'above'; else, w = 'below'; end
end


function [b, ci, R2, p] = olsLocal(x, y)
    ok = isfinite(x) & isfinite(y);
    x = x(ok); y = y(ok);
    n = numel(y);
    X = [ones(n,1) x];
    b = X \ y;
    r = y - X*b;
    s2 = (r'*r) / (n-2);
    V = s2 * inv(X'*X);                                  %#ok<MINV>
    se = sqrt(diag(V));
    tcrit = tinv(0.975, n-2);
    ci = [b - tcrit*se, b + tcrit*se];
    R2 = 1 - (r'*r) / sum((y - mean(y)).^2);
    tstat = b(2) / se(2);
    p = 2 * (1 - tcdf(abs(tstat), n-2));
end


function d = deseason(x, moL)
% Remove the calendar-month climatology. NaN in, NaN out, so the caller's
% isfinite mask still selects the same months.
    d = x;
    for m = 1:12
        sel = (moL == m) & isfinite(x);
        if any(sel); d(sel) = x(sel) - mean(x(sel)); end
    end
end


function [Tmax, kBreak, p] = snhtLocalMC(x, nMC)
    x = x(:); x = x(isfinite(x));
    [Tmax, kBreak] = snhtStatLocal(x);
    n = numel(x); count = 0;
    for b = 1:nMC
        if snhtStatLocal(randn(n,1)) >= Tmax; count = count + 1; end
    end
    p = (count + 1) / (nMC + 1);
end


function [Tmax, kBreak] = snhtStatLocal(x)
    n  = numel(x);
    z  = (x - mean(x)) / std(x);
    cs = cumsum(z);
    total = cs(end);
    Tk = zeros(n-1,1);
    for k = 1:n-1
        Tk(k) = k*(cs(k)/k)^2 + (n-k)*((total-cs(k))/(n-k))^2;
    end
    [Tmax, kBreak] = max(Tk);
end


function [Umax, kBreak, p] = pettittLocal(x)
    x = x(:); x = x(isfinite(x)); n = numel(x);
    r = tiedrank(x);
    U = zeros(n,1);
    for k = 1:n
        U(k) = 2*sum(r(1:k)) - k*(n+1);
    end
    [Umax, kBreak] = max(abs(U));
    p = min(2 * exp(-6*Umax^2 / (n^3 + n^2)), 1);
end


function d = monthsBetween(ym1, ym2)
    y1 = floor(ym1/100); m1 = mod(ym1,100);
    y2 = floor(ym2/100); m2 = mod(ym2,100);
    d = (y1 - y2)*12 + (m1 - m2);
end


function c=pick(vn,cands)
for i=1:numel(cands), j=find(vn==cands(i),1); if ~isempty(j), c=char(vn(j)); return; end, end
error('Missing required column. Candidates: %s. Available: %s',strjoin(cellstr(cands),', '),strjoin(cellstr(vn),', '));
end


function d=toDailyDatetime(x)
if isdatetime(x), d=dateshift(x,'start','day');
elseif isnumeric(x), d=dateshift(datetime(x,'ConvertFrom','datenum'),'start','day');
else, d=dateshift(datetime(x),'start','day'); end
end


function [b,ci,p,R2,n]=fitRH(y,rh,mo,dn,k,L)
[tf,loc]=ismember(dn-k,dn); iy=find(tf); ip=loc(tf);
yy=y(iy); xx=rh(ip); mm=mo(iy); dates=dn(iy);
g=isfinite(yy)&isfinite(xx)&isfinite(mm);
yy=yy(g); xx=xx(g); mm=mm(g); dates=dates(g);
[dates,o]=sort(dates); yy=yy(o); xx=xx(o); mm=mm(o);
n=numel(yy); D=zeros(n,11);
for m=2:12, D(:,m-1)=double(mm==m); end
X=[ones(n,1),xx,D];
if rank(X)<size(X,2), error('Rank-deficient RH-only model at lag %d.',k); end
bvec=X\yy; e=yy-X*bvec; V=dateNW(X,e,dates,L);
se=sqrt(max(diag(V),0)); df=n-size(X,2); tc=tinv(.975,df);
b=bvec(2); ci=[b-tc*se(2),b+tc*se(2)]; p=2*tcdf(-abs(b/se(2)),df);
R2=1-sum(e.^2)/sum((yy-mean(yy)).^2);
end


function V=dateNW(X,e,dn,L)
scores=X.*e; S=scores'*scores;
for j=1:L
[tf,loc]=ismember(dn-j,dn); i=find(tf); ip=loc(tf);
if isempty(i), continue; end
G=scores(i,:)'*scores(ip,:); w=1-j/(L+1); S=S+w*(G+G');
end
A=pinv(X'*X); V=A*S*A; V=(V+V')/2;
end


function [b,ci,p,R2,n]=fitSeasonModel(y,rain,rh,mo,dn,sel,L)
    good=sel&isfinite(y)&isfinite(rain)&isfinite(rh)&isfinite(mo);
    yy=y(good); xr=rain(good); xh=rh(good); mm=mo(good); dates=dn(good);
    [dates,ix]=sort(dates); yy=yy(ix); xr=xr(ix); xh=xh(ix); mm=mm(ix);
    n=numel(yy);
    % Fixed effects for calendar months represented in this season.
    um=unique(mm);
    D=zeros(n,max(numel(um)-1,0));
    for j=2:numel(um), D(:,j-1)=double(mm==um(j)); end
    X=[ones(n,1),xr,xh,D];
    if rank(X)<size(X,2), error('Seasonal model matrix is rank deficient.'); end
    b=X\yy; e=yy-X*b;
    V=dateAwareNW(X,e,dates,L);
    se=sqrt(max(diag(V),0));
    df=n-size(X,2);
    tc=tinv(0.975,df);
    ci=[b-tc.*se,b+tc.*se];
    p=2*tcdf(-abs(b./se),df);
    R2=1-sum(e.^2)/sum((yy-mean(yy)).^2);
end


function V=dateAwareNW(X,e,dn,L)
    scores=X.*e;
    S=scores'*scores;
    for j=1:L
        [tf,loc]=ismember(dn-j,dn);
        i=find(tf); ip=loc(tf);
        if isempty(i), continue; end
        G=scores(i,:)'*scores(ip,:);
        w=1-j/(L+1);
        S=S+w*(G+G');
    end
    XtX_inv=pinv(X'*X);
    V=XtX_inv*S*XtX_inv;
    V=(V+V')/2;
end


function [b,ci,p,R2,n] = fitLaggedModel(y,rain,rh,mo,dn,k,L)
    % Pair response at date t with predictors at exact date t-k.
    [tf,loc] = ismember(dn-k,dn);
    iy = find(tf);
    ip = loc(tf);
    yy = y(iy); xr = rain(ip); xh = rh(ip); mm = mo(iy); dates = dn(iy);
    good = isfinite(yy)&isfinite(xr)&isfinite(xh)&isfinite(mm);
    yy=yy(good); xr=xr(good); xh=xh(good); mm=mm(good); dates=dates(good);
    [dates,ix]=sort(dates); yy=yy(ix); xr=xr(ix); xh=xh(ix); mm=mm(ix);
    n=numel(yy);
    D=zeros(n,11);
    for m=2:12, D(:,m-1)=double(mm==m); end
    X=[ones(n,1),xr,xh,D];
    if rank(X)<size(X,2), error('Model matrix is rank deficient at lag %d.',k); end
    b=X\yy; e=yy-X*b;
    V=dateAwareNW(X,e,dates,L);
    se=sqrt(max(diag(V),0));
    df=n-size(X,2);
    tcrit=tinv(0.975,df);
    ci=[b-tcrit.*se,b+tcrit.*se];
    p=2*tcdf(-abs(b./se),df);
    R2=1-sum(e.^2)/sum((yy-mean(yy)).^2);
end


function d=toDaily(x)
if isdatetime(x),d=dateshift(x,'start','day');elseif isnumeric(x),d=dateshift(datetime(x,'ConvertFrom','datenum'),'start','day');else,d=dateshift(datetime(x),'start','day');end
end


function [b,ci,p,R2,n]=fitSeason(y,rain,rh,mo,dn,sel,L)
g=sel&isfinite(y)&isfinite(rain)&isfinite(rh)&isfinite(mo);yy=y(g);xr=rain(g);xh=rh(g);mm=mo(g);dates=dn(g);
[dates,o]=sort(dates);yy=yy(o);xr=xr(o);xh=xh(o);mm=mm(o);n=numel(yy);
um=unique(mm);D=zeros(n,max(numel(um)-1,0));for j=2:numel(um),D(:,j-1)=double(mm==um(j));end
X=[ones(n,1),xr,xh,D];if rank(X)<size(X,2),error('Rank-deficient seasonal model.');end
b=X\yy;e=yy-X*b;V=dateNW(X,e,dates,L);se=sqrt(max(diag(V),0));df=n-size(X,2);tc=tinv(.975,df);
ci=[b-tc.*se,b+tc.*se];p=2*tcdf(-abs(b./se),df);R2=1-sum(e.^2)/sum((yy-mean(yy)).^2);
end

function T = standardizeDailyHeaders(T,stationName)
% Standardize station-sheet headers while preserving the original data.
    raw = string(T.Properties.VariableNames);
    norm = lower(regexprep(raw,'[^a-zA-Z0-9]',''));
    targets = { ...
        'Date',      ["date","tanggal","tgl","datetime","observationdate"]; ...
        'Tmean',     ["tmean","meantemperature","temperaturemean"]; ...
        'Tmin',      ["tmin","mintemperature","minimumtemperature"]; ...
        'Tmax',      ["tmax","maxtemperature","maximumtemperature"]; ...
        'Rain',      ["rain","rainfall","precipitation","rr","curahhujan"]; ...
        'Hummidity', ["hummidity","humidity","rh","relativehumidity", ...
                      "kelembapan","kelembaban"]};
    for k = 1:size(targets,1)
        target = targets{k,1};
        aliases = targets{k,2};
        ix = find(ismember(norm,aliases),1,'first');
        if isempty(ix)
            error('%s sheet: cannot identify column %s. Available headers: %s', ...
                stationName,target,strjoin(cellstr(raw),', '));
        end
        T.Properties.VariableNames{ix} = target;
        raw(ix) = target;
        norm(ix) = lower(regexprep(string(target),'[^a-zA-Z0-9]',''));
    end
    T.Date = parseDailyDate(T.Date,stationName);
end


function d = parseDailyDate(x,stationName)
% Parse Excel serial dates, MATLAB datetimes, and common Indonesian date text.
    if isdatetime(x)
        d = dateshift(x(:),'start','day');
        return
    end
    if isnumeric(x)
        d = dateshift(datetime(x(:),'ConvertFrom','excel'),'start','day');
        return
    end
    sx = strtrim(string(x(:)));
    formats = {'dd/MM/yyyy','d/M/yyyy','yyyy-MM-dd','dd-MM-yyyy', ...
        'd-M-yyyy','MM/dd/yyyy','dd/MM/yyyy HH:mm:ss', ...
        'yyyy-MM-dd HH:mm:ss','yyyy-MM-dd''T''HH:mm:ss'};
    d = NaT(size(sx));
    for k = 1:numel(formats)
        try
            candidate = datetime(sx,'InputFormat',formats{k});
            if all(~isnat(candidate))
                d = dateshift(candidate,'start','day');
                return
            end
        catch
        end
    end
    try
        candidate = datetime(sx);
        if all(~isnat(candidate))
            d = dateshift(candidate,'start','day');
            return
        end
    catch
    end
    error('%s sheet: date values could not be parsed. Check the date format.', ...
        stationName);
end


function [n,b,se,p,ci,R2,pOLS,usedDates] = fitRHOnlyIntegrated( ...
    y,rh,monthVec,dn,lagDay,HAC_LAG,selection)
% RH-only OLS with calendar-month fixed effects and date-aware Newey-West.
% Pair response at day t with RH at the exact calendar date t-lagDay.
    if nargin < 7 || isempty(selection)
        selection = true(size(y));
    end
    y = y(:); rh = rh(:); monthVec = monthVec(:); dn = dn(:);
    selection = selection(:);
    [tf,loc] = ismember(dn-lagDay,dn);
    iy = find(tf);
    ip = loc(tf);
    yy = y(iy);
    xx = rh(ip);
    mm = monthVec(iy);
    dates = dn(iy);
    ss = selection(iy);
    good = isfinite(yy)&isfinite(xx)&isfinite(mm)&ss;
    yy=yy(good); xx=xx(good); mm=mm(good); dates=dates(good);
    [dates,ord] = sort(dates);
    yy=yy(ord); xx=xx(ord); mm=mm(ord);
    n = numel(yy);
    if n < 20
        error('Too few valid observations (%d) for RH-only model.',n);
    end

    % Use only month indicators represented in the selected sample.
    um = unique(mm);
    D = zeros(n,max(numel(um)-1,0));
    for j = 2:numel(um)
        D(:,j-1) = double(mm==um(j));
    end
    X = [ones(n,1),xx,D];
    k = size(X,2);
    if rank(X) < k
        error('RH-only design matrix is rank deficient at lag %d.',lagDay);
    end

    bhat = X\yy;
    resid = yy-X*bhat;
    XtXi = pinv(X'*X);
    df = n-k;
    if df <= 0
        error('Non-positive residual degrees of freedom in RH-only model.');
    end

    % OLS standard error, retained as a diagnostic.
    s2 = sum(resid.^2)/df;
    seOLS = sqrt(max(s2*XtXi(2,2),0));
    pOLS = 2*tcdf(-abs(bhat(2)/seOLS),df);

    % Newey-West Bartlett HAC, matched by exact calendar-day separation.
    meat = dateAwareMeat(X,resid,dates,HAC_LAG);
    V = XtXi*meat*XtXi*(n/df);
    V = (V+V')/2;
    seAll = sqrt(max(diag(V),0));
    se = seAll(2);
    if ~isfinite(se) || se<=0
        error('Invalid HAC standard error in RH-only model at lag %d.',lagDay);
    end
    p = 2*tcdf(-abs(bhat(2)/se),df);
    tcrit = tinv(0.975,df);
    b = bhat(2);
    ci = [b-tcrit*se,b+tcrit*se];
    R2 = 1-sum(resid.^2)/sum((yy-mean(yy)).^2);
    usedDates = dates;
end


function meat = dateAwareMeat(X,resid,dn,L)
% Bartlett HAC "meat"; only pairs exactly L calendar days apart contribute.
    scores = X.*resid;
    meat = scores'*scores;
    for lag = 1:L
        [tf,loc] = ismember(dn-lag,dn);
        i = find(tf);
        ip = loc(tf);
        if isempty(i), continue; end
        G = scores(i,:)'*scores(ip,:);
        w = 1-lag/(L+1);
        meat = meat+w*(G+G');
    end
end


function [xAnom,yAnom,usedDates] = monthAnomalyPairs(y,rh,monthVec,dn,lagDay)
% Matched response/predictor pairs, demeaned within response calendar month.
    [tf,loc] = ismember(dn-lagDay,dn);
    iy = find(tf);
    ip = loc(tf);
    yy = y(iy);
    xx = rh(ip);
    mm = monthVec(iy);
    usedDates = dn(iy);
    good = isfinite(yy)&isfinite(xx)&isfinite(mm);
    yy=yy(good); xx=xx(good); mm=mm(good); usedDates=usedDates(good);
    yAnom=yy; xAnom=xx;
    for m = 1:12
        q = mm==m;
        if any(q)
            yAnom(q)=yy(q)-mean(yy(q));
            xAnom(q)=xx(q)-mean(xx(q));
        end
    end
end


function [mu,ciHalf] = meanHACCI(x,dn,L)
% Newey-West confidence interval for a binned mean, using exact calendar-day lags.
    x=x(:); dn=dn(:); good=isfinite(x)&isfinite(dn); x=x(good); dn=dn(good);
    [dn,ord]=sort(dn); x=x(ord); n=numel(x); mu=NaN; ciHalf=NaN;
    if n<3, return; end
    e=x-mean(x); X=ones(n,1); V=dateAwareNW(X,e,dn,L);
    se=sqrt(max(V(1,1),0));
    ciHalf=tinv(0.975,n-1)*se; mu=mean(x);
end


function [xb,yb,eb,nb] = quantileBinSummary(x,y,nBins,minN)
% Quantile-binned means and approximate 95% confidence intervals.
    x=x(:); y=y(:);
    good=isfinite(x)&isfinite(y);
    x=x(good); y=y(good);
    xb=NaN(0,1); yb=NaN(0,1); eb=NaN(0,1); nb=zeros(0,1);
    if numel(x)<minN || numel(unique(x))<2, return; end
    edges=unique(quantile(x,linspace(0,1,nBins+1)));
    if numel(edges)<2, return; end
    bin=discretize(x,edges);
    % Include observations equal to the final edge in the last bin.
    bin(x==edges(end))=numel(edges)-1;
    nActual=numel(edges)-1;
    xb=NaN(nActual,1); yb=xb; eb=xb; nb=zeros(nActual,1);
    for j=1:nActual
        q=bin==j;
        nb(j)=sum(q);
        if nb(j)>=minN
            xb(j)=mean(x(q));
            yb(j)=mean(y(q));
            eb(j)=1.96*std(y(q))/sqrt(nb(j));
        end
    end
end

function out = subsetETCCDIYears(in, keepYears)
% Keep selected annual ETCCDI rows while preserving station thresholds/metadata.
% This helper is used to enforce a common set of complete pre-break years.
    annualFields = {'Year','RX1day','RX5day','SDII','CDD','PRCPTOT', ...
                    'R95p','R99p','nR95day','nR99day', ...
                    'nValidDays','nExpectedDays','isCompleteYear'};
    out = in;
    keepYears = keepYears(:);
    rowMask = ismember(in.Year, keepYears);
    for k = 1:numel(annualFields)
        fieldName = annualFields{k};
        if isfield(in, fieldName)
            values = in.(fieldName);
            out.(fieldName) = values(rowMask,:);
        end
    end
    out.incompleteYears = [];
end

