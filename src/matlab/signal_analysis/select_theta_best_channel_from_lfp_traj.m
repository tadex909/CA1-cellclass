%% Select the best theta channel from one LFP file and its trajdata file
%
% This script reproduces the first useful part of theta.m: it finds the LFP
% channel with the strongest theta/delta envelope ratio. The important
% difference is that it keeps the LFP at its native sampling rate and uses
% the behavior sample index clock only to choose the LFP interval.
%
% Alignment convention used here:
%   traj__start and traj__stop are 0-based behavior sample indices.
%   behavior index 0 corresponds to LFP sample 1.
%
% Therefore:
%   lfp_sample_index = round((behavior_index / Fs_beh) * Fs_lfp) + 1

clearvars;

%% Paths and parameters

repoRoot = 'C:\Users\tadse\OneDrive\Documenti\GitHub\CA1-cellclass';

lfpPath = fullfile(repoRoot, 'data', 'raw', 'lfp', ...
    'VS103_2024-11-10_14-58-13_lfp1250Hz.mat');

trajPath = fullfile(repoRoot, 'data', 'interim', 'VS103', '2024-11-10', ...
    'VS103_2024-11-10_14-58-13_trajdata.npz');

outputPath = fullfile(repoRoot, 'data', 'interim', 'VS103', '2024-11-10', ...
    'VS103_2024-11-10_14-58-13_theta_best_channel.mat');

Fs_lfp = 1250;
Fs_beh = 1000;
Fs_analysis = 1000;

thetaBand = [4 9];
deltaBand = [0.5 4];

% theta.m chooses the best channel from prm.win1, i.e. condition 1.
% Set params.useConditionOnly = false if you prefer to use the full session.
params.useConditionOnly = true;
params.conditionToUse = 1;
params.plotResult = true;
params.plotRandomTrial = true;
params.saveResult = true;
params.randomSeed = 3;
params.randomTrial = [];
params.computeTimeFrequency = true;
params.timeFrequencyRange = 0.5:0.5:20;
params.timeFrequencyNorm = 'freq';
params.timeFrequencySmoothSec = 0.10;

%% Add the local NPY reader used by the legacy MATLAB pipeline

npyReaderDir = fullfile(repoRoot, 'src', 'placefields', 'matlab', ...
    'FINAL_TEAMPIPELINE');
addpath(npyReaderDir);

assert(exist(lfpPath, 'file') == 2, 'LFP file not found: %s', lfpPath);
assert(exist(trajPath, 'file') == 2, 'trajdata file not found: %s', trajPath);
assert(exist('readNPY', 'file') == 2, ...
    'Could not find readNPY.m. Expected it under: %s', npyReaderDir);

%% Load the minimal trajdata fields

traj = read_minimal_traj_npz(trajPath);

startsBeh = double(traj.start(:));
stopsBeh = double(traj.stop(:));
cond = double(traj.cond(:));

if params.useConditionOnly
    selectedTrials = cond == params.conditionToUse;
else
    selectedTrials = true(size(cond));
end

assert(any(selectedTrials), 'No trials selected for best-channel analysis.');

selectedStarts = startsBeh(selectedTrials);
selectedStops = stopsBeh(selectedTrials);

%% Convert behavior sample indices into LFP sample indices

trialStartLfpIdx = floor((selectedStarts ./ Fs_beh) .* Fs_lfp) + 1;
trialStopLfpIdx = ceil((selectedStops ./ Fs_beh) .* Fs_lfp) + 1;

lfpInfo = whos('-file', lfpPath, 'allfp_ds');
assert(~isempty(lfpInfo), 'Variable allfp_ds was not found in: %s', lfpPath);

nChannels = lfpInfo.size(1);
nLfpSamples = lfpInfo.size(2);

trialStartLfpIdx = max(1, min(nLfpSamples, trialStartLfpIdx));
trialStopLfpIdx = max(1, min(nLfpSamples, trialStopLfpIdx));

idxStart = min(trialStartLfpIdx);
idxStop = max(trialStopLfpIdx);

assert(idxStop > idxStart, 'Selected LFP interval is empty.');

fprintf('Using LFP samples %d:%d (%.3f to %.3f s in LFP clock).\n', ...
    idxStart, idxStop, (idxStart - 1) / Fs_lfp, (idxStop - 1) / Fs_lfp);
fprintf('Selected %d/%d trials', nnz(selectedTrials), numel(selectedTrials));
if params.useConditionOnly
    fprintf(' from condition %d', params.conditionToUse);
end
fprintf('.\n');

%% Load only the selected LFP interval if possible

channelLabels = [];

try
    lfpMat = matfile(lfpPath);
    lfpSegment = double(lfpMat.allfp_ds(:, idxStart:idxStop));
    if ismember('channels', who(lfpMat))
        channelLabels = lfpMat.channels;
    end
catch
    loaded = load(lfpPath, 'allfp_ds', 'channels');
    lfpSegment = double(loaded.allfp_ds(:, idxStart:idxStop));
    if isfield(loaded, 'channels')
        channelLabels = loaded.channels;
    end
end

if isempty(channelLabels)
    channelLabels = (1:nChannels).';
else
    channelLabels = channelLabels(:);
end

%% Compute theta/delta envelope ratio per channel

[Btheta, Atheta] = butter(2, thetaBand ./ (Fs_lfp / 2), 'bandpass');
[Bdelta, Adelta] = butter(2, deltaBand ./ (Fs_lfp / 2), 'bandpass');

thetaEnvelopeMedian = nan(nChannels, 1);
deltaEnvelopeMedian = nan(nChannels, 1);

for ch = 1:nChannels
    x = lfpSegment(ch, :).';
    x = x - median(x, 'omitnan');

    thetaEnvelope = abs(hilbert(filtfilt(Btheta, Atheta, x)));
    deltaEnvelope = abs(hilbert(filtfilt(Bdelta, Adelta, x)));

    thetaEnvelopeMedian(ch) = median(thetaEnvelope, 'omitnan');
    deltaEnvelopeMedian(ch) = median(deltaEnvelope, 'omitnan');
end

thetaDeltaRatio = thetaEnvelopeMedian ./ deltaEnvelopeMedian;
[bestRatio, bestChannelRow] = max(thetaDeltaRatio);
bestChannelLabel = channelLabels(bestChannelRow);

fprintf('Best channel row in allfp_ds: %d\n', bestChannelRow);
fprintf('Best channel label from channels: %g\n', bestChannelLabel);
fprintf('Theta/delta median-envelope ratio: %.4f\n', bestRatio);

%% Compute theta-filtered LFP and theta amplitude on the best channel

prm = struct();
prm.fq_lfp = Fs_lfp;
prm.fd = Fs_analysis;

v = lfpSegment(bestChannelRow, :).';
v = v - median(v, 'omitnan');

[B_theta, A_theta] = butter(2, thetaBand ./ (prm.fq_lfp / 2), 'bandpass');
lfp_theta = filtfilt(B_theta, A_theta, v);

[~, ~, lfp_theta_mod] = theta_hilbert_envelope(lfp_theta, prm.fq_lfp);

theta = struct();
theta.lfpd = downsample_theta_signal(v, prm.fq_lfp, prm.fd);
theta.lfp_thetad = downsample_theta_signal(lfp_theta, prm.fq_lfp, prm.fd);
theta.lfp_thetad_mod = downsample_theta_signal(lfp_theta_mod, prm.fq_lfp, prm.fd);
theta.time = ((0:numel(theta.lfp_thetad_mod)-1).' ./ prm.fd) + ...
    ((idxStart - 1) ./ prm.fq_lfp);

fprintf('Computed theta amplitude on best channel at %.0f Hz.\n', prm.fd);

%% Pick one trial and plot theta amplitude, position, and speed

plotTrial = struct();

if params.plotRandomTrial
    candidateTrials = find(selectedTrials);

    if isempty(params.randomTrial)
        if isempty(params.randomSeed)
            rng('shuffle');
        else
            rng(params.randomSeed);
        end
        trialNumber = candidateTrials(randi(numel(candidateTrials)));
    else
        trialNumber = params.randomTrial;
        assert(ismember(trialNumber, candidateTrials), ...
            'params.randomTrial must be one of the selected trials.');
    end

    trialTime = get_trial_vector(traj.time, trialNumber);
    trialPosition = get_trial_vector(traj.position, trialNumber);
    trialSpeed = get_trial_vector(traj.speed, trialNumber);

    nTrial = min([numel(trialTime), numel(trialPosition), numel(trialSpeed)]);
    trialTime = trialTime(1:nTrial);
    trialPosition = trialPosition(1:nTrial);
    trialSpeed = trialSpeed(1:nTrial);

    trialLfpClockTime = startsBeh(trialNumber) ./ Fs_beh + trialTime;
    trialThetaAmplitude = interp1(theta.time, theta.lfp_thetad_mod, ...
        trialLfpClockTime, 'linear', NaN);

    plotTrial.trialNumber = trialNumber;
    plotTrial.condition = cond(trialNumber);
    plotTrial.startBehaviorIndex = startsBeh(trialNumber);
    plotTrial.stopBehaviorIndex = stopsBeh(trialNumber);
    plotTrial.timeFromTrialStart = trialTime;
    plotTrial.lfpClockTime = trialLfpClockTime;
    plotTrial.thetaAmplitude = trialThetaAmplitude;
    plotTrial.position = trialPosition;
    plotTrial.speed = trialSpeed;

    fprintf('Plotting trial %d, condition %g.\n', trialNumber, cond(trialNumber));

    if params.computeTimeFrequency
        trialStartLfpIdx = floor((startsBeh(trialNumber) ./ Fs_beh) .* Fs_lfp) + 1;
        trialStopLfpIdx = ceil((stopsBeh(trialNumber) ./ Fs_beh) .* Fs_lfp) + 1;
        trialStartLfpIdx = max(1, min(nLfpSamples, trialStartLfpIdx));
        trialStopLfpIdx = max(1, min(nLfpSamples, trialStopLfpIdx));

        relStart = max(1, trialStartLfpIdx - idxStart + 1);
        relStop = min(numel(v), trialStopLfpIdx - idxStart + 1);
        assert(relStop > relStart, 'Selected trial LFP segment is empty.');

        trialLfp = v(relStart:relStop);
        trialLfpd = downsample_theta_signal(trialLfp, prm.fq_lfp, prm.fd);
        trialLfpStartTime = ((trialStartLfpIdx - 1) ./ prm.fq_lfp) - ...
            (startsBeh(trialNumber) ./ Fs_beh);
        timeFrequencyTime = trialLfpStartTime + ...
            ((0:numel(trialLfpd)-1).' ./ prm.fd);

        [timeFrequencyPower, timeFrequencyFreq] = compute_time_frequency( ...
            trialLfpd, prm.fd, params.timeFrequencyRange, ...
            params.timeFrequencyNorm);
        timeFrequencyPower = smooth_matrix_time(timeFrequencyPower, ...
            round(params.timeFrequencySmoothSec .* prm.fd));

        plotTrial.timeFrequencyTime = timeFrequencyTime;
        plotTrial.timeFrequencyFreq = timeFrequencyFreq;
        plotTrial.timeFrequencyPower = timeFrequencyPower;
        plotTrial.timeFrequencyNorm = params.timeFrequencyNorm;

        fprintf('Computed trial time-frequency map from %.1f to %.1f Hz.\n', ...
            min(timeFrequencyFreq), max(timeFrequencyFreq));
    end

    figure('Color', 'w');
    if isfield(plotTrial, 'timeFrequencyPower')
        tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

        nexttile;
        imagesc(plotTrial.timeFrequencyTime, plotTrial.timeFrequencyFreq, ...
            plotTrial.timeFrequencyPower);
        axis xy;
        hold on;
        xl = xlim;
        plot(xl, [thetaBand(1) thetaBand(1)], 'w--', 'LineWidth', 1);
        plot(xl, [thetaBand(2) thetaBand(2)], 'w--', 'LineWidth', 1);
        xlim(xl);
        colorLimits = robust_color_limits(plotTrial.timeFrequencyPower);
        if all(isfinite(colorLimits))
            clim(colorLimits);
        end
        ylabel('freq. (Hz)');
        title(sprintf('Trial %d, condition %g, best channel row %d', ...
            trialNumber, cond(trialNumber), bestChannelRow));
        colorbar;
    else
        tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    end

    nexttile;
    plot(trialTime, trialThetaAmplitude, 'k', 'LineWidth', 1);
    ylabel('theta amp.');
    if ~isfield(plotTrial, 'timeFrequencyPower')
        title(sprintf('Trial %d, condition %g, best channel row %d', ...
            trialNumber, cond(trialNumber), bestChannelRow));
    end
    box off;

    nexttile;
    plot(trialTime, trialPosition, 'Color', [0.1 0.35 0.8], 'LineWidth', 1);
    ylabel('position');
    box off;

    nexttile;
    plot(trialTime, trialSpeed, 'Color', [0.75 0.2 0.1], 'LineWidth', 1);
    ylabel('speed');
    xlabel('time from trial start (s)');
    box off;
end

%% Save and plot

result = struct();
result.lfpPath = lfpPath;
result.trajPath = trajPath;
result.Fs_lfp = Fs_lfp;
result.Fs_beh = Fs_beh;
result.Fs_analysis = Fs_analysis;
result.thetaBand = thetaBand;
result.deltaBand = deltaBand;
result.useConditionOnly = params.useConditionOnly;
result.conditionToUse = params.conditionToUse;
result.selectedTrials = selectedTrials;
result.trialStartLfpIdx = trialStartLfpIdx;
result.trialStopLfpIdx = trialStopLfpIdx;
result.idxStart = idxStart;
result.idxStop = idxStop;
result.channelLabels = channelLabels;
result.thetaEnvelopeMedian = thetaEnvelopeMedian;
result.deltaEnvelopeMedian = deltaEnvelopeMedian;
result.thetaDeltaRatio = thetaDeltaRatio;
result.bestChannelRow = bestChannelRow;
result.bestChannelLabel = bestChannelLabel;
result.bestRatio = bestRatio;
result.theta = theta;
result.plotTrial = plotTrial;

if params.saveResult
    outputDir = fileparts(outputPath);
    if ~exist(outputDir, 'dir')
        mkdir(outputDir);
    end
    save(outputPath, 'result');
    fprintf('Saved result to: %s\n', outputPath);
end

if params.plotResult
    figure('Color', 'w');
    bar(thetaDeltaRatio);
    hold on;
    yl = ylim;
    plot([bestChannelRow bestChannelRow], yl, 'r-', 'LineWidth', 1.5);
    ylim(yl);
    xlabel('row in allfp\_ds');
    ylabel('theta/delta median-envelope ratio');
    title(sprintf('Best theta channel: row %d, label %g', ...
        bestChannelRow, bestChannelLabel));
    box off;
end

%% Local helpers

function traj = read_minimal_traj_npz(npzPath)
    tmpDir = tempname;
    mkdir(tmpDir);
    cleanupTmp = onCleanup(@() cleanup_dir(tmpDir));

    unzip(npzPath, tmpDir);

    traj = struct();
    traj.cond = readNPY(fullfile(tmpDir, 'traj__Cond.npy'));
    traj.start = readNPY(fullfile(tmpDir, 'traj__start.npy'));
    traj.stop = readNPY(fullfile(tmpDir, 'traj__stop.npy'));
    traj.time = read_npz_json_field(tmpDir, 'traj__time__json.npy');
    traj.position = read_npz_json_field(tmpDir, 'traj__VRtraj__json.npy');

    xSpeedFile = fullfile(tmpDir, 'traj__XSpeed__json.npy');
    speedFile = fullfile(tmpDir, 'traj__Speed__json.npy');
    if exist(xSpeedFile, 'file')
        traj.speed = read_npz_json_field(tmpDir, 'traj__XSpeed__json.npy');
    elseif exist(speedFile, 'file')
        traj.speed = read_npz_json_field(tmpDir, 'traj__Speed__json.npy');
    else
        error('Could not find traj__XSpeed__json.npy or traj__Speed__json.npy.');
    end

    clear cleanupTmp;
end

function values = read_npz_json_field(tmpDir, filename)
    jsonText = read_npy_string_scalar(fullfile(tmpDir, filename));
    values = jsondecode(jsonText);
end

function text = read_npy_string_scalar(filename)
    fid = fopen(filename, 'r', 'ieee-le');
    assert(fid ~= -1, 'Could not open file: %s', filename);
    cleanupFile = onCleanup(@() fclose(fid));

    magicString = fread(fid, [1 6], 'uint8=>uint8');
    assert(all(magicString == [147 78 85 77 80 89]), ...
        'File is not a valid NPY file: %s', filename);

    majorVersion = fread(fid, 1, 'uint8=>double');
    fread(fid, 1, 'uint8=>double');

    if majorVersion == 1
        headerLength = fread(fid, 1, 'uint16=>double');
    else
        headerLength = fread(fid, 1, 'uint32=>double');
    end

    header = fread(fid, [1 headerLength], '*char');
    token = regexp(header, '''descr''\s*:\s*''\|S(\d+)''', ...
        'tokens', 'once');
    assert(~isempty(token), 'Expected a scalar NPY byte string in: %s', filename);

    nBytes = str2double(token{1});
    raw = fread(fid, [1 nBytes], 'uint8=>uint8');
    raw = raw(raw ~= 0);

    try
        text = native2unicode(raw, 'UTF-8');
    catch
        text = char(raw);
    end

    clear cleanupFile;
end

function values = get_trial_vector(fieldValues, trialNumber)
    if iscell(fieldValues)
        values = fieldValues{trialNumber};
    elseif isnumeric(fieldValues)
        if isvector(fieldValues)
            values = fieldValues(:);
        else
            values = fieldValues(trialNumber, :).';
        end
    else
        error('Unsupported trial field type: %s', class(fieldValues));
    end

    values = double(values(:));
    values = values(~isnan(values));
end

function [phase, frequency, amplitude] = theta_hilbert_envelope(x, fs)
    if exist('fct_thetahilb', 'file') == 2
        [phase, frequency, amplitude] = fct_thetahilb(x, fs);
    else
        analyticSignal = hilbert(x);
        phase = angle(analyticSignal);
        amplitude = abs(analyticSignal);
        unwrappedPhase = unwrap(phase);
        frequency = [NaN; diff(unwrappedPhase) .* fs ./ (2*pi)];
    end
end

function y = downsample_theta_signal(x, fsIn, fsOut)
    if exist('fct_dsampling', 'file') == 2
        y = fct_dsampling(x, fsIn, fsOut);
    elseif fsIn == fsOut
        y = x;
    else
        [p, q] = rat(fsOut ./ fsIn);
        y = resample(x, p, q);
    end

    y = y(:);
end

function [tfPower, freq] = compute_time_frequency(x, fs, freqRange, normMode)
    x = double(x(:));
    x = x - median(x, 'omitnan');

    freqRange = freqRange(:);
    freqRange = freqRange(freqRange > 0 & freqRange < fs/2);
    assert(~isempty(freqRange), 'No valid frequencies requested.');

    if exist('fct_cwt', 'file') == 2
        [tfRaw, freq] = fct_cwt(x, 'fs', fs, 'frange', freqRange(:).');
        tfPower = abs(tfRaw);
        freq = freq(:);
    elseif exist('cwt', 'file') ~= 0
        [wt, freq] = cwt(x, fs, 'FrequencyLimits', ...
            [min(freqRange) max(freqRange)]);
        tfPower = abs(wt);
        [freq, order] = sort(freq(:));
        tfPower = tfPower(order, :);
        tfPower = interp1(freq, tfPower, freqRange, 'linear', NaN);
        freq = freqRange;
    else
        [tfPower, freq] = bandpass_hilbert_time_frequency(x, fs, freqRange);
    end

    tfPower = normalize_time_frequency(tfPower, freq, normMode);
end

function [tfPower, freq] = bandpass_hilbert_time_frequency(x, fs, freqRange)
    freq = freqRange(:);
    tfPower = nan(numel(freq), numel(x));

    if numel(freq) > 1
        halfWidth = median(diff(freq)) ./ 2;
    else
        halfWidth = 0.25;
    end

    for k = 1:numel(freq)
        lo = max(0.05, freq(k) - halfWidth);
        hi = min((fs/2) * 0.99, freq(k) + halfWidth);

        if lo <= 0.05
            [B, A] = butter(2, hi ./ (fs/2), 'low');
        else
            [B, A] = butter(2, [lo hi] ./ (fs/2), 'bandpass');
        end

        tfPower(k, :) = abs(hilbert(filtfilt(B, A, x)));
    end
end

function tfPower = normalize_time_frequency(tfPower, freq, normMode)
    switch lower(normMode)
        case 'freq'
            tfPower = tfPower .* sqrt(freq(:));
        case 'zscore'
            mu = mean(tfPower, 2, 'omitnan');
            sigma = std(tfPower, 0, 2, 'omitnan');
            sigma(sigma == 0 | isnan(sigma)) = 1;
            tfPower = (tfPower - mu) ./ sigma;
        case 'none'
        otherwise
            error('Unknown time-frequency normalization: %s', normMode);
    end
end

function dataSmooth = smooth_matrix_time(data, halfWindow)
    if halfWindow <= 1
        dataSmooth = data;
        return;
    end

    win = -halfWindow:halfWindow;
    kernel = exp(-(win .^ 2) ./ (halfWindow / 2) ^ 2);
    kernel = kernel ./ sum(kernel);
    dataSmooth = conv2(data, kernel, 'same');
end

function limits = robust_color_limits(data)
    values = data(:);
    values = values(isfinite(values));

    if isempty(values)
        limits = [NaN NaN];
        return;
    end

    values = sort(values);
    loIdx = max(1, round(0.02 * numel(values)));
    hiIdx = min(numel(values), round(0.98 * numel(values)));
    limits = [values(loIdx) values(hiIdx)];

    if limits(1) == limits(2)
        limits = [min(values) max(values)];
    end
end

function cleanup_dir(tmpDir)
    if exist(tmpDir, 'dir')
        rmdir(tmpDir, 's');
    end
end
