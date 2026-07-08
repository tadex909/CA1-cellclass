function result = detect_nonactive_theta_epochs(varargin)
%DETECT_NONACTIVE_THETA_EPOCHS Detect awake low-speed, low-theta intervals.
%
% This function identifies candidate non-active on-track periods for awake
% reactivation analyses. The detector is deliberately conservative:
%
%   1. pick the best theta channel from the raw LFP,
%   2. compute a wavelet theta/delta ratio on that channel,
%   3. estimate the low-theta threshold only from low-speed samples,
%   4. keep continuous periods that are both low-speed and low-theta,
%   5. split those periods into 20 ms bins for ReactivationStrength.
%
% Default usage, for the VS103 example session:
%
%   result = detect_nonactive_theta_epochs;
%
% Common overrides:
%
%   result = detect_nonactive_theta_epochs( ...
%       'thetaBand', [4 9], ...
%       'highpassCutoffHz', 0.1, ...
%       'speedThreshold', 2, ...
%       'thresholdPercentile', 20, ...
%       'analysisFs', 250);
%
% Important alignment convention:
%   traj__start and traj__stop are 0-based behavior sample indices, and
%   behavior index 0 corresponds to LFP sample 1.
%   By default, each trial is analyzed only from traj__start to endVR.

%% Parse inputs and add local project paths

scriptDir = fileparts(mfilename('fullpath'));
repoRoot = find_repo_root(scriptDir);
add_project_paths(repoRoot);

defaultLfpPath = fullfile(repoRoot, 'data', 'raw', 'lfp', ...
    'VS103_2024-11-10_14-58-13_lfp1250Hz.mat'); % Directory of downsampled LFP data
defaultTrajPath = fullfile(repoRoot, 'data', 'interim', 'VS103', ...
    '2024-11-10', 'VS103_2024-11-10_14-58-13_trajdata.npz'); % Directory of trajectory data

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'lfpPath', defaultLfpPath, @is_text_scalar);
addParameter(p, 'trajPath', defaultTrajPath, @is_text_scalar);
addParameter(p, 'outputPath', '', @is_text_scalar);
addParameter(p, 'Fs_lfp', 1250, @is_positive_scalar); %Argument is the sampling frequency of the LFP data in Hz
addParameter(p, 'Fs_beh', 1000, @is_positive_scalar); %Argument is the sampling frequency of the behavior data in Hz
addParameter(p, 'analysisFs', 250, @is_positive_scalar); %Argument is the sampling frequency of the analysis grid in Hz
addParameter(p, 'thetaBand', [6 9], @is_two_element_numeric); %Theta band frequencies in Hz
addParameter(p, 'deltaBand', [1 4], @is_two_element_numeric); %Delta band frequencies in Hz
addParameter(p, 'highpassCutoffHz', 0.1, ...
    @(x) isempty(x) || is_nonnegative_scalar(x)); %Frequency cutoff for hig-pass filter to remove slow drifts from the LFP
addParameter(p, 'highpassOrder', 2, @is_positive_scalar);
addParameter(p, 'waveletFreqs', 0.5:0.5:20, @(x) isnumeric(x) && isvector(x)); %Frequencies used for the spectrogram
addParameter(p, 'speedThreshold', 2, @is_positive_scalar); %Speed threshold in cm/s for low-speed detection
addParameter(p, 'thresholdPercentile', 20, @is_percent_scalar); %Percentile of low-speed theta/delta values used to define the low-theta threshold
addParameter(p, 'minEpochDuration', 0.5, @is_nonnegative_scalar);
addParameter(p, 'mergeGap', 0.25, @is_nonnegative_scalar); %If two intervals are separated by less than this many seconds, they will be merged into one interval
addParameter(p, 'binSize', 0.020, @is_positive_scalar); %Size of the bins used for ReactivationStrength, in seconds
addParameter(p, 'analyzeConditions', [], @(x) isempty(x) || isnumeric(x));
addParameter(p, 'limitTrialsToEndVR', true, @is_logical_scalar);
addParameter(p, 'bestChannelCondition', [], @(x) isempty(x) || isnumeric(x)); % Empty means use all analyzed conditions for best-channel selection
addParameter(p, 'bestChannelRow', [], @(x) isempty(x) || is_positive_scalar(x));
addParameter(p, 'bestChannelMaxDuration', 300, @is_positive_or_inf_scalar);
addParameter(p, 'useAbsoluteSpeed', true, @is_logical_scalar);
addParameter(p, 'computeMultitaper', false, @is_logical_scalar);
addParameter(p, 'multitaperMovingWin', [1.0 0.05], @is_two_element_numeric);
addParameter(p, 'multitaperTapers', [3 5], @(x) isnumeric(x) && isvector(x));
addParameter(p, 'multitaperFpass', [0.5 20], @is_two_element_numeric);
addParameter(p, 'makeFigures', true, @is_logical_scalar);
addParameter(p, 'saveFigures', true, @is_logical_scalar);
addParameter(p, 'saveResult', true, @is_logical_scalar);
addParameter(p, 'exampleTrial', [], @(x) isempty(x) || is_positive_scalar(x));
parse(p, varargin{:});

cfg = p.Results;
cfg.lfpPath = char(cfg.lfpPath);
cfg.trajPath = char(cfg.trajPath);
cfg.outputPath = char(cfg.outputPath);
cfg.thetaBand = sort(double(cfg.thetaBand(:).'));
cfg.deltaBand = sort(double(cfg.deltaBand(:).'));
if ~isempty(cfg.highpassCutoffHz)
    cfg.highpassCutoffHz = double(cfg.highpassCutoffHz);
end
cfg.highpassOrder = round(double(cfg.highpassOrder));
cfg.waveletFreqs = sort(double(cfg.waveletFreqs(:).'));
cfg.multitaperMovingWin = double(cfg.multitaperMovingWin(:).');
cfg.multitaperTapers = double(cfg.multitaperTapers(:).');
cfg.multitaperFpass = sort(double(cfg.multitaperFpass(:).'));

assert(exist(cfg.lfpPath, 'file') == 2, 'LFP file not found: %s', cfg.lfpPath);
assert(exist(cfg.trajPath, 'file') == 2, ...
    'Trajectory file not found: %s', cfg.trajPath);
assert(cfg.thetaBand(1) > cfg.deltaBand(1), ...
    'Expected thetaBand to sit above deltaBand.');

if isempty(cfg.outputPath)
    cfg.outputPath = default_output_path(cfg.lfpPath, cfg.trajPath);
end

%% Load trial metadata and define the LFP analysis span

traj = read_traj_file(cfg.trajPath);

startsBeh = double(traj.start(:));
stopsBeh = double(traj.stop(:));
cond = double(traj.cond(:));
nTrials = numel(cond);

assert(numel(startsBeh) == nTrials && numel(stopsBeh) == nTrials, ...
    'Trajectory start/stop/condition vectors have inconsistent lengths.');

endVrBeh = nan(size(startsBeh));
endVrFraction = nan(size(startsBeh));
if cfg.limitTrialsToEndVR
    [endVrBeh, endVrFraction] = endvr_behavior_index(traj, startsBeh, stopsBeh);
    analysisStopsBeh = endVrBeh;
else
    analysisStopsBeh = stopsBeh;
end

selectedTrials = select_trials(cond, cfg.analyzeConditions);
selectedTrials = selectedTrials & isfinite(startsBeh) & isfinite(stopsBeh) & ...
    isfinite(analysisStopsBeh) & analysisStopsBeh > startsBeh;
assert(any(selectedTrials), 'No valid trials selected for analysis.');

lfpInfo = whos('-file', cfg.lfpPath, 'allfp_ds');
assert(~isempty(lfpInfo), 'Variable allfp_ds was not found in: %s', cfg.lfpPath);
nChannels = lfpInfo.size(1);
nLfpSamples = lfpInfo.size(2);

channelLabels = load_channel_labels(cfg.lfpPath, nChannels);

trialStartLfpIdx = floor((startsBeh ./ cfg.Fs_beh) .* cfg.Fs_lfp) + 1;
trialStopLfpIdx = ceil((analysisStopsBeh ./ cfg.Fs_beh) .* cfg.Fs_lfp) + 1;
trialStartLfpIdx = clamp_index(trialStartLfpIdx, nLfpSamples);
trialStopLfpIdx = clamp_index(trialStopLfpIdx, nLfpSamples);

analysisStartIdx = min(trialStartLfpIdx(selectedTrials));
analysisStopIdx = max(trialStopLfpIdx(selectedTrials));
assert(analysisStopIdx > analysisStartIdx, 'The selected LFP interval is empty.');

fprintf('Analysis span: LFP samples %d:%d (%.3f to %.3f s).\n', ...
    analysisStartIdx, analysisStopIdx, ...
    (analysisStartIdx - 1) / cfg.Fs_lfp, ...
    (analysisStopIdx - 1) / cfg.Fs_lfp);
if cfg.limitTrialsToEndVR
    fprintf('Trial pooling is limited to start -> endVR for %d selected trials.\n', ...
        nnz(selectedTrials));
end

%% Select the best theta channel

lfpSource = open_lfp_source(cfg.lfpPath);

bestChannelInfo = struct();
if isempty(cfg.bestChannelRow)
    bestChannelTrials = select_trials(cond, cfg.bestChannelCondition);
    bestChannelTrials = bestChannelTrials & selectedTrials;
    if ~any(bestChannelTrials)
        warning(['No trials matched bestChannelCondition. Using all selected ' ...
            'trials for best-channel scoring.']);
        bestChannelTrials = selectedTrials;
    end

    bestChannelInfo = select_best_theta_channel( ...
        lfpSource, nChannels, channelLabels, ...
        trialStartLfpIdx(bestChannelTrials), ...
        trialStopLfpIdx(bestChannelTrials), cfg);
else
    assert(cfg.bestChannelRow <= nChannels, ...
        'bestChannelRow exceeds the number of LFP channels.');
    bestChannelInfo.row = round(cfg.bestChannelRow);
    bestChannelInfo.label = channelLabels(bestChannelInfo.row);
    bestChannelInfo.ratio = NaN;
    bestChannelInfo.thetaEnvelopeMedian = NaN(nChannels, 1);
    bestChannelInfo.deltaEnvelopeMedian = NaN(nChannels, 1);
    bestChannelInfo.thetaDeltaRatio = NaN(nChannels, 1);
end

fprintf('Best theta channel row: %d', bestChannelInfo.row);
if isfinite(bestChannelInfo.label)
    fprintf(' (label %g)', bestChannelInfo.label);
end
if isfinite(bestChannelInfo.ratio)
    fprintf(', theta/delta envelope ratio %.4f', bestChannelInfo.ratio);
end
fprintf('.\n');

%% Load best-channel LFP and put it on the analysis grid

bestLfpNative = load_lfp_contiguous( ...
    lfpSource, bestChannelInfo.row, analysisStartIdx, analysisStopIdx);
bestLfpNative = double(bestLfpNative(:));
bestLfpNative = bestLfpNative - median_omitnan(bestLfpNative);
bestLfpNative = highpass_lfp_if_requested( ...
    bestLfpNative, cfg.Fs_lfp, cfg.highpassCutoffHz, cfg.highpassOrder);

bestLfpAnalysis = resample_signal(bestLfpNative, cfg.Fs_lfp, cfg.analysisFs);
bestLfpAnalysis = bestLfpAnalysis(:);

time = ((0:numel(bestLfpAnalysis)-1).' ./ cfg.analysisFs) + ...
    ((analysisStartIdx - 1) ./ cfg.Fs_lfp);

fprintf('Spectral analysis grid: %.1f Hz, %d samples.\n', ...
    cfg.analysisFs, numel(time));

%% Align behavior to the same time base

[speed, position, trialNumberAtTime] = build_behavior_traces( ...
    traj, selectedTrials, time, cfg.Fs_beh, cfg.useAbsoluteSpeed, ...
    analysisStopsBeh);

validBehavior = isfinite(speed);
fprintf('Behavior covers %.1f s of the %.1f s analysis span.\n', ...
    nnz(validBehavior) / cfg.analysisFs, numel(time) / cfg.analysisFs);

%% Compute wavelet theta/delta ratio

assert(all(cfg.waveletFreqs > 0) && all(cfg.waveletFreqs < cfg.analysisFs / 2), ...
    'waveletFreqs must be between 0 and analysisFs/2.');

[cwtd, ffreq] = fct_cwt(bestLfpAnalysis, ...
    'fs', cfg.analysisFs, 'frange', cfg.waveletFreqs);
ffreq = ffreq(:);

waveletPower = single(abs(cwtd)).^2;
clear cwtd;

thetaIdx = ffreq >= cfg.thetaBand(1) & ffreq <= cfg.thetaBand(2);
deltaIdx = ffreq >= cfg.deltaBand(1) & ffreq <= cfg.deltaBand(2);
assert(any(thetaIdx), 'No wavelet frequencies fall inside thetaBand.');
assert(any(deltaIdx), 'No wavelet frequencies fall inside deltaBand.');

thetaPower = mean(double(waveletPower(thetaIdx, :)), 1).';
deltaPower = mean(double(waveletPower(deltaIdx, :)), 1).';
thetaDeltaDb = 10 .* log10(thetaPower ./ (deltaPower + eps));

%% Threshold low-speed theta/delta values and find continuous periods

validMask = validBehavior & isfinite(thetaDeltaDb) & isfinite(thetaPower) & ...
    isfinite(deltaPower);
lowSpeedMask = validMask & speed < cfg.speedThreshold;
assert(any(lowSpeedMask), ...
    'No low-speed samples found. Increase speedThreshold or check speed units.');

thresholdDb = percentile_omitnan(thetaDeltaDb(lowSpeedMask), ...
    cfg.thresholdPercentile);
candidateMask = lowSpeedMask & thetaDeltaDb < thresholdDb;

rawSeq = mask_to_sequences(candidateMask);
mergeGapSamples = round(cfg.mergeGap * cfg.analysisFs);
minEpochSamples = max(1, ceil(cfg.minEpochDuration * cfg.analysisFs));

mergedSeq = merge_sequences(rawSeq, mergeGapSamples);
finalSeq = remove_short_sequences(mergedSeq, minEpochSamples);

nonActiveIntervals = sequences_to_intervals(finalSeq, time, cfg.analysisFs);
binsNonActive = split_intervals_fixed(nonActiveIntervals, cfg.binSize);

if isempty(nonActiveIntervals)
    warning(['No non-active intervals remained after merging gaps and ' ...
        'removing short events. Diagnostics will still be saved.']);
end

fprintf('Low-speed samples: %.1f s. Threshold: %.3f dB (%gth percentile).\n', ...
    nnz(lowSpeedMask) / cfg.analysisFs, thresholdDb, cfg.thresholdPercentile);
fprintf('Detected %d non-active intervals and %d %.0f-ms bins.\n', ...
    size(nonActiveIntervals, 1), size(binsNonActive, 1), cfg.binSize * 1000);

%% Compute optional multitaper diagnostic

multitaper = struct();
if cfg.computeMultitaper
    multitaper = compute_multitaper_ratio(bestLfpAnalysis, time, cfg);
end

%% Pack and save results

trialInfo = struct();
trialInfo.condition = cond;
trialInfo.selected = selectedTrials;
trialInfo.startBehaviorIndex = startsBeh;
trialInfo.stopBehaviorIndex = analysisStopsBeh;
trialInfo.fullStopBehaviorIndex = stopsBeh;
trialInfo.endVRBehaviorIndex = endVrBeh;
trialInfo.endVRFraction = endVrFraction;
trialInfo.startTime = startsBeh ./ cfg.Fs_beh;
trialInfo.stopTime = analysisStopsBeh ./ cfg.Fs_beh;
trialInfo.fullStopTime = stopsBeh ./ cfg.Fs_beh;
if isfield(traj, 'tstart')
    trialInfo.tstart = double(traj.tstart(:));
end
if isfield(traj, 'tstop')
    trialInfo.tstop = double(traj.tstop(:));
end
if isfield(traj, 'endVR')
    trialInfo.endVR = double(traj.endVR(:));
end
if isfield(traj, 'WB')
    trialInfo.WB = traj.WB;
end
trialInfo.startLfpIndex = trialStartLfpIdx;
trialInfo.stopLfpIndex = trialStopLfpIdx;

result = struct();
result.lfpPath = cfg.lfpPath; %Path to the LFP file used for analysis
result.trajPath = cfg.trajPath; %Path to the trajectory file used for analysis
result.nonActiveIntervals = nonActiveIntervals; %Intervals of low-speed, low-theta epochs measured in seconds
result.binsNonActive = binsNonActive; 
result.candidateMask = candidateMask; %Sampled with frequency cfg.analysisFs, this is a logical mask of candidate low-speed, low-theta samples
result.lowSpeedMask = lowSpeedMask;
result.validMask = validMask;
result.time = time;
result.speed = speed;
result.position = position;
result.trialNumberAtTime = trialNumberAtTime;
result.thetaPower = thetaPower;
result.deltaPower = deltaPower;
result.thetaDeltaDb = thetaDeltaDb;
result.thresholdDb = thresholdDb;
result.waveletFreq = ffreq;
result.bestChannel = bestChannelInfo;
result.bestChannelRow = bestChannelInfo.row;
result.bestChannelLabel = bestChannelInfo.label;
result.trials = trialInfo;
result.multitaper = multitaper;
result.params = cfg;

if cfg.makeFigures
    result.figurePaths = make_diagnostic_figures( ...
        result, waveletPower, cfg);
else
    result.figurePaths = struct();
end

if cfg.saveResult
    outputDir = fileparts(cfg.outputPath);
    if ~isempty(outputDir) && exist(outputDir, 'dir') ~= 7
        mkdir(outputDir);
    end
    save(cfg.outputPath, 'result', '-v7.3');
    fprintf('Saved non-active epoch result to: %s\n', cfg.outputPath);
end

end

%% Path helpers

function repoRoot = find_repo_root(startDir)
repoRoot = startDir;
while true
    if exist(fullfile(repoRoot, 'CA1-cellclass.prj'), 'file') == 2 || ...
            exist(fullfile(repoRoot, '.git'), 'dir') == 7
        return;
    end

    parentDir = fileparts(repoRoot);
    if isempty(parentDir) || strcmp(parentDir, repoRoot)
        repoRoot = fileparts(fileparts(fileparts(startDir)));
        return;
    end
    repoRoot = parentDir;
end
end

function add_project_paths(repoRoot)
safe_addpath(fullfile(repoRoot, 'src', 'matlab', 'signal_analysis'));
safe_addpath(fullfile(repoRoot, 'src', 'matlab', 'FINAL_TEAMPIPELINE'));
safe_addpath(fullfile(repoRoot, 'src', 'matlab', 'LFP_analysis'));
safe_addpath(fullfile(repoRoot, 'src', 'neurocode', 'utilities', 'intervalsC+'));

chronuxRoot = fullfile(repoRoot, 'src', 'neurocode', ...
    'SpectralAnalyses', 'chronux', 'spectral_analysis');
safe_addpath(fullfile(chronuxRoot, 'helper'));
safe_addpath(fullfile(chronuxRoot, 'continuous'));
end

function safe_addpath(pathToAdd)
if exist(pathToAdd, 'dir') == 7
    addpath(pathToAdd);
end
end

function outputPath = default_output_path(lfpPath, trajPath)
outputDir = fileparts(trajPath);
if isempty(outputDir)
    outputDir = pwd;
end

[~, lfpName] = fileparts(lfpPath);
lfpName = regexprep(lfpName, '_lfp1250Hz$', '');
outputPath = fullfile(outputDir, [lfpName '_nonactive_theta_epochs.mat']);
end

%% Trajectory loading

function traj = read_traj_file(trajPath)
[~, ~, ext] = fileparts(trajPath);
switch lower(ext)
    case '.npz'
        traj = read_minimal_traj_npz(trajPath);
    case '.mat'
        traj = read_minimal_traj_mat(trajPath);
    otherwise
        error('Unsupported trajectory file type: %s', ext);
end
end

function traj = read_minimal_traj_npz(npzPath)
assert(exist('readNPY', 'file') == 2, ...
    ['Could not find readNPY.m. The function adds the expected repo path ' ...
    'automatically, so check that src/matlab/FINAL_TEAMPIPELINE exists.']);

tmpDir = tempname;
mkdir(tmpDir);
cleanupTmp = onCleanup(@() cleanup_dir(tmpDir));

unzip(npzPath, tmpDir);

traj = struct();
traj.cond = readNPY(fullfile(tmpDir, 'traj__Cond.npy'));
traj.start = readNPY(fullfile(tmpDir, 'traj__start.npy'));
traj.stop = readNPY(fullfile(tmpDir, 'traj__stop.npy'));
traj.tstart = readNPY(fullfile(tmpDir, 'traj__tstart.npy'));
traj.tstop = readNPY(fullfile(tmpDir, 'traj__tstop.npy'));
traj.endVR = readNPY(fullfile(tmpDir, 'traj__endVR.npy'));
wbFile = fullfile(tmpDir, 'traj__WB.npy');
wbJsonFile = fullfile(tmpDir, 'traj__WB__json.npy');
if exist(wbFile, 'file') == 2
    try
        traj.WB = readNPY(wbFile);
    catch
        traj.WB = read_npy_text_array(wbFile);
    end
elseif exist(wbJsonFile, 'file') == 2
    traj.WB = read_npz_json_field(tmpDir, 'traj__WB__json.npy');
end
traj.time = read_npz_json_field(tmpDir, 'traj__time__json.npy');
traj.position = read_npz_json_field(tmpDir, 'traj__VRtraj__json.npy');

xSpeedFile = fullfile(tmpDir, 'traj__XSpeed__json.npy');
speedFile = fullfile(tmpDir, 'traj__Speed__json.npy');
if exist(xSpeedFile, 'file') == 2
    traj.speed = read_npz_json_field(tmpDir, 'traj__XSpeed__json.npy');
elseif exist(speedFile, 'file') == 2
    traj.speed = read_npz_json_field(tmpDir, 'traj__Speed__json.npy');
else
    error('Could not find traj__XSpeed__json.npy or traj__Speed__json.npy.');
end

clear cleanupTmp;
end

function traj = read_minimal_traj_mat(matPath)
loaded = load(matPath);

if isfield(loaded, 'Traj')
    Traj = loaded.Traj;
    traj = struct();
    traj.cond = field_vector(Traj, 'Cond');
    traj.start = field_vector(Traj, 'start');
    traj.stop = field_vector(Traj, 'stop');
    traj.tstart = field_vector(Traj, 'tstart');
    traj.tstop = field_vector(Traj, 'tstop');
    traj.endVR = field_vector(Traj, 'endVR');
    if isfield(Traj, 'WB')
        traj.WB = field_optional_numeric_or_cell(Traj, 'WB');
    end
    traj.time = field_cell(Traj, 'time');

    if isfield(Traj, 'VRtraj')
        traj.position = field_cell(Traj, 'VRtraj');
    elseif isfield(loaded, 'X_ds_n')
        traj.position = split_continuous_by_trial(loaded.X_ds_n, traj);
    else
        error('MAT trajectory file does not contain VRtraj or X_ds_n.');
    end

    if isfield(Traj, 'XSpeed')
        traj.speed = field_cell(Traj, 'XSpeed');
    elseif isfield(Traj, 'Speed')
        traj.speed = field_cell(Traj, 'Speed');
    elseif isfield(loaded, 'XSpeed')
        traj.speed = split_continuous_by_trial(loaded.XSpeed, traj);
    else
        error('MAT trajectory file does not contain XSpeed or Speed.');
    end
else
    error('MAT trajectory file must contain a Traj variable.');
end
end

function values = field_vector(structArray, fieldName)
assert(isfield(structArray, fieldName), ...
    'Trajectory struct is missing required field: %s', fieldName);
values = double([structArray.(fieldName)]).';
end

function values = field_cell(structArray, fieldName)
assert(isfield(structArray, fieldName), ...
    'Trajectory struct is missing required field: %s', fieldName);
values = cell(numel(structArray), 1);
for k = 1:numel(structArray)
    values{k} = double(structArray(k).(fieldName)(:));
end
end

function values = field_optional_numeric_or_cell(structArray, fieldName)
raw = {structArray.(fieldName)}.';
if all(cellfun(@(x) isnumeric(x) && isscalar(x), raw))
    values = cellfun(@double, raw);
else
    values = raw;
end
end

function values = split_continuous_by_trial(vector, traj)
vector = double(vector(:));
values = cell(numel(traj.start), 1);
for k = 1:numel(values)
    idx = (round(traj.start(k)) + 1):round(traj.stop(k));
    idx = idx(idx >= 1 & idx <= numel(vector));
    values{k} = vector(idx);
end
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

function values = read_npy_text_array(filename)
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
descrToken = regexp(header, '''descr''\s*:\s*''([<|>]?[US]\d+)''', ...
    'tokens', 'once');
assert(~isempty(descrToken), 'Expected a text NPY array in: %s', filename);

dtype = descrToken{1};
charsPerItem = str2double(regexp(dtype, '\d+', 'match', 'once'));

shapeToken = regexp(header, '''shape''\s*:\s*\((.*?)\)', ...
    'tokens', 'once');
assert(~isempty(shapeToken), 'Could not read NPY shape from: %s', filename);
shapeText = regexprep(shapeToken{1}, '[L\s]', '');
shapeText = regexprep(shapeText, ',$', '');
shape = str2num(shapeText); %#ok<ST2NM>
if isempty(shape)
    shape = 1;
end
if isscalar(shape)
    shape = [shape 1];
end
nItems = prod(shape);

values = cell(nItems, 1);
if contains(dtype, 'U')
    raw = fread(fid, nItems * charsPerItem, 'uint32=>uint32');
    raw = reshape(raw, charsPerItem, nItems).';
    for k = 1:nItems
        codes = raw(k, :);
        codes = codes(codes ~= 0);
        values{k} = char(codes);
    end
else
    raw = fread(fid, [charsPerItem nItems], 'uint8=>uint8').';
    for k = 1:nItems
        bytes = raw(k, :);
        bytes = bytes(bytes ~= 0);
        values{k} = char(bytes);
    end
end

values = reshape(values, shape);
clear cleanupFile;
end

function cleanup_dir(tmpDir)
if exist(tmpDir, 'dir') == 7
    rmdir(tmpDir, 's');
end
end

%% LFP loading and best-channel selection

function labels = load_channel_labels(lfpPath, nChannels)
labels = (1:nChannels).';
vars = who('-file', lfpPath);
if any(strcmp(vars, 'channels'))
    loaded = load(lfpPath, 'channels');
    labels = double(loaded.channels(:));
end
if numel(labels) ~= nChannels
    labels = (1:nChannels).';
end
end

function source = open_lfp_source(lfpPath)
source = struct();
source.path = lfpPath;
try
    m = matfile(lfpPath);
    m.allfp_ds(1, 1);
    source.kind = 'matfile';
    source.obj = m;
catch
    loaded = load(lfpPath, 'allfp_ds');
    source.kind = 'memory';
    source.allfp_ds = loaded.allfp_ds;
end
end

function info = select_best_theta_channel( ...
    lfpSource, nChannels, channelLabels, startIdx, stopIdx, cfg)

[startIdx, stopIdx, scoreDuration] = limit_scoring_segments( ...
    startIdx, stopIdx, cfg.Fs_lfp, cfg.bestChannelMaxDuration);

fprintf('Best-channel scoring uses %.1f s of selected LFP data.\n', ...
    scoreDuration);

[Btheta, Atheta] = butter(2, cfg.thetaBand ./ (cfg.Fs_lfp / 2), 'bandpass');
[Bdelta, Adelta] = butter(2, cfg.deltaBand ./ (cfg.Fs_lfp / 2), 'bandpass');

thetaEnvelopeMedian = nan(nChannels, 1);
deltaEnvelopeMedian = nan(nChannels, 1);
scoreSamples = sum(stopIdx - startIdx + 1);
maxMatrixElements = 60e6;

if scoreSamples * nChannels <= maxMatrixElements
    scoreLfp = load_lfp_segments_all_channels( ...
        lfpSource, 1:nChannels, startIdx, stopIdx);
else
    scoreLfp = [];
    warning(['Best-channel scoring matrix would be large; reading one ' ...
        'channel at a time.']);
end

for ch = 1:nChannels
    if isempty(scoreLfp)
        x = load_lfp_segments(lfpSource, ch, startIdx, stopIdx);
    else
        x = scoreLfp(ch, :).';
    end
    x = double(x(:));
    x = x - median_omitnan(x);

    thetaEnvelope = abs(hilbert(filtfilt(Btheta, Atheta, x)));
    deltaEnvelope = abs(hilbert(filtfilt(Bdelta, Adelta, x)));

    thetaEnvelopeMedian(ch) = median_omitnan(thetaEnvelope);
    deltaEnvelopeMedian(ch) = median_omitnan(deltaEnvelope);
end

clear scoreLfp;

thetaDeltaRatio = thetaEnvelopeMedian ./ (deltaEnvelopeMedian + eps);
[bestRatio, bestRow] = max(thetaDeltaRatio);

info = struct();
info.row = bestRow;
info.label = channelLabels(bestRow);
info.ratio = bestRatio;
info.thetaEnvelopeMedian = thetaEnvelopeMedian;
info.deltaEnvelopeMedian = deltaEnvelopeMedian;
info.thetaDeltaRatio = thetaDeltaRatio;
info.scoreDuration = scoreDuration;
info.scoreStartLfpIndex = startIdx;
info.scoreStopLfpIndex = stopIdx;
end

function [startIdx, stopIdx, scoreDuration] = limit_scoring_segments( ...
    startIdx, stopIdx, fs, maxDuration)

startIdx = round(startIdx(:));
stopIdx = round(stopIdx(:));
lengths = stopIdx - startIdx + 1;
totalSamples = sum(lengths);

if isinf(maxDuration)
    scoreDuration = totalSamples ./ fs;
    return;
end

maxSamples = round(maxDuration .* fs);
if totalSamples <= maxSamples
    scoreDuration = totalSamples ./ fs;
    return;
end

fraction = maxSamples ./ totalSamples;
keepLengths = max(1, floor(lengths .* fraction));

while sum(keepLengths) > maxSamples
    [~, idx] = max(keepLengths);
    keepLengths(idx) = keepLengths(idx) - 1;
end

for k = 1:numel(startIdx)
    center = round((startIdx(k) + stopIdx(k)) ./ 2);
    halfWidth = floor((keepLengths(k) - 1) ./ 2);
    newStart = center - halfWidth;
    newStop = newStart + keepLengths(k) - 1;

    if newStart < startIdx(k)
        newStart = startIdx(k);
        newStop = newStart + keepLengths(k) - 1;
    end
    if newStop > stopIdx(k)
        newStop = stopIdx(k);
        newStart = newStop - keepLengths(k) + 1;
    end

    startIdx(k) = newStart;
    stopIdx(k) = newStop;
end

scoreDuration = sum(stopIdx - startIdx + 1) ./ fs;
end

function x = load_lfp_contiguous(lfpSource, channel, startIdx, stopIdx)
x = read_lfp(lfpSource, channel, startIdx, stopIdx);
x = x(:);
end

function x = load_lfp_segments(lfpSource, channel, startIdx, stopIdx)
startIdx = round(startIdx(:));
stopIdx = round(stopIdx(:));
lengths = stopIdx - startIdx + 1;
assert(all(lengths > 0), 'Invalid LFP segment boundaries.');

x = zeros(sum(lengths), 1);
cursor = 1;
for k = 1:numel(startIdx)
    segment = read_lfp(lfpSource, channel, startIdx(k), stopIdx(k));
    segment = segment(:);
    idx = cursor:(cursor + numel(segment) - 1);
    x(idx) = segment;
    cursor = cursor + numel(segment);
end
end

function x = load_lfp_segments_all_channels(lfpSource, channels, startIdx, stopIdx)
startIdx = round(startIdx(:));
stopIdx = round(stopIdx(:));
lengths = stopIdx - startIdx + 1;
assert(all(lengths > 0), 'Invalid LFP segment boundaries.');

x = zeros(numel(channels), sum(lengths));
cursor = 1;
for k = 1:numel(startIdx)
    segment = read_lfp(lfpSource, channels, startIdx(k), stopIdx(k));
    idx = cursor:(cursor + size(segment, 2) - 1);
    x(:, idx) = segment;
    cursor = cursor + size(segment, 2);
end
end

function x = read_lfp(lfpSource, channel, startIdx, stopIdx)
switch lfpSource.kind
    case 'matfile'
        x = double(lfpSource.obj.allfp_ds(channel, startIdx:stopIdx));
    case 'memory'
        x = double(lfpSource.allfp_ds(channel, startIdx:stopIdx));
    otherwise
        error('Unknown LFP source kind: %s', lfpSource.kind);
end
end

function idx = clamp_index(idx, nSamples)
idx = round(idx);
idx = max(1, min(nSamples, idx));
end

%% Behavior alignment

function selected = select_trials(cond, wantedConditions)
if isempty(wantedConditions)
    selected = true(size(cond));
else
    selected = ismember(cond, wantedConditions(:));
end
end

function [endVrBeh, frac] = endvr_behavior_index(traj, startsBeh, stopsBeh)
requiredFields = {'tstart', 'tstop', 'endVR'};
for k = 1:numel(requiredFields)
    assert(isfield(traj, requiredFields{k}), ...
        'Trajectory data is missing required field for endVR cutoff: %s', ...
        requiredFields{k});
end

tstart = double(traj.tstart(:));
tstop = double(traj.tstop(:));
endVR = double(traj.endVR(:));

assert(numel(tstart) == numel(startsBeh) && numel(tstop) == numel(startsBeh) && ...
    numel(endVR) == numel(startsBeh), ...
    'Trajectory tstart/tstop/endVR vectors have inconsistent lengths.');

duration = tstop - tstart;
frac = (endVR - tstart) ./ duration;
frac(~isfinite(frac) | duration <= 0) = NaN;
frac = max(0, min(1, frac));

endVrBeh = startsBeh + round(frac .* (stopsBeh - startsBeh));
endVrBeh = max(startsBeh, min(stopsBeh, endVrBeh));
end

function [speed, position, trialNumberAtTime] = build_behavior_traces( ...
    traj, selectedTrials, time, Fs_beh, useAbsoluteSpeed, analysisStopsBeh)

speed = nan(size(time));
position = nan(size(time));
trialNumberAtTime = nan(size(time));

trialNumbers = find(selectedTrials(:).');
for trialNumber = trialNumbers
    trialTime = get_trial_vector(traj.time, trialNumber);
    trialSpeed = get_trial_vector(traj.speed, trialNumber);
    trialPosition = get_trial_vector(traj.position, trialNumber);

    n = min([numel(trialTime), numel(trialSpeed), numel(trialPosition)]);
    if n < 2
        continue;
    end

    trialTime = trialTime(1:n);
    trialSpeed = trialSpeed(1:n);
    trialPosition = trialPosition(1:n);

    if useAbsoluteSpeed
        trialSpeed = abs(trialSpeed);
    end

    absoluteTime = double(traj.start(trialNumber)) ./ Fs_beh + trialTime;
    cutoffTime = double(analysisStopsBeh(trialNumber)) ./ Fs_beh;
    keep = absoluteTime <= cutoffTime;
    absoluteTime = absoluteTime(keep);
    trialSpeed = trialSpeed(keep);
    trialPosition = trialPosition(keep);

    if numel(absoluteTime) < 2
        continue;
    end

    [absoluteTime, uniqueIdx] = unique(absoluteTime, 'stable');
    trialSpeed = trialSpeed(uniqueIdx);
    trialPosition = trialPosition(uniqueIdx);

    idx = time >= absoluteTime(1) & time <= absoluteTime(end);
    if ~any(idx)
        continue;
    end

    speed(idx) = interp1(absoluteTime, trialSpeed, time(idx), 'linear', NaN);
    position(idx) = interp1(absoluteTime, trialPosition, time(idx), 'linear', NaN);
    trialNumberAtTime(idx) = trialNumber;
end
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
values = values(isfinite(values));
end

function y = resample_signal(x, fsIn, fsOut)
x = double(x(:));
if fsIn == fsOut
    y = x;
    return;
end

if exist('resample', 'file') == 2
    [p, q] = rat(fsOut ./ fsIn, 1e-12);
    y = resample(x, p, q);
else
    warning(['resample.m was not found. Falling back to linear ' ...
        'interpolation for LFP analysis resampling.']);
    sourceTime = (0:numel(x)-1).' ./ fsIn;
    targetTime = (0:(1/fsOut):sourceTime(end)).';
    y = interp1(sourceTime, x, targetTime, 'linear', 'extrap');
end
end

function y = highpass_lfp_if_requested(x, fs, cutoffHz, filterOrder)
y = double(x(:));

if isempty(cutoffHz) || cutoffHz <= 0
    return;
end

assert(cutoffHz < fs / 2, ...
    'highpassCutoffHz must be below Nyquist frequency.');

filterOrder = max(1, round(filterOrder));
[B, A] = butter(filterOrder, cutoffHz ./ (fs / 2), 'high');

minLength = 3 * (max(numel(A), numel(B)) - 1);
if numel(y) <= minLength
    warning(['LFP segment is too short for high-pass filtfilt. ' ...
        'Leaving it unfiltered.']);
    return;
end

y = filtfilt(B, A, y);
end

%% Interval operations

function seq = mask_to_sequences(mask)
mask = logical(mask(:));
edges = diff([false; mask; false]);
starts = find(edges == 1);
stops = find(edges == -1) - 1;
seq = [starts stops];
end

function merged = merge_sequences(seq, maxGapSamples)
if isempty(seq)
    merged = zeros(0, 2);
    return;
end

merged = seq(1, :);
for k = 2:size(seq, 1)
    gap = seq(k, 1) - merged(end, 2) - 1;
    if gap <= maxGapSamples
        merged(end, 2) = seq(k, 2);
    else
        merged(end + 1, :) = seq(k, :); %#ok<AGROW>
    end
end
end

function kept = remove_short_sequences(seq, minSamples)
if isempty(seq)
    kept = zeros(0, 2);
    return;
end

seqLength = seq(:, 2) - seq(:, 1) + 1;
kept = seq(seqLength >= minSamples, :);
end

function intervals = sequences_to_intervals(seq, time, fs)
if isempty(seq)
    intervals = zeros(0, 2);
    return;
end

startTimes = time(seq(:, 1));
stopTimes = time(seq(:, 2)) + 1 ./ fs;
intervals = [startTimes(:) stopTimes(:)];
end

function bins = split_intervals_fixed(intervals, binSize)
if isempty(intervals)
    bins = zeros(0, 2);
    return;
end

piecesPerInterval = floor((intervals(:, 2) - intervals(:, 1)) ./ binSize);
nBins = sum(piecesPerInterval);
bins = zeros(nBins, 2);

cursor = 1;
for k = 1:size(intervals, 1)
    nPieces = piecesPerInterval(k);
    if nPieces <= 0
        continue;
    end

    starts = intervals(k, 1) + (0:nPieces-1).' .* binSize;
    idx = cursor:(cursor + nPieces - 1);
    bins(idx, :) = [starts starts + binSize];
    cursor = cursor + nPieces;
end
end

%% Multitaper diagnostic

function multitaper = compute_multitaper_ratio(lfpAnalysis, time, cfg)
multitaper = struct();

try
    mtParams = struct();
    mtParams.Fs = cfg.analysisFs;
    mtParams.fpass = cfg.multitaperFpass;
    mtParams.tapers = cfg.multitaperTapers;
    mtParams.pad = 0;
    mtParams.err = 0;
    mtParams.trialave = 1;

    [S, mtTimeRelative, mtFreq] = mtspecgramc( ...
        lfpAnalysis, cfg.multitaperMovingWin, mtParams);

    mtFreq = mtFreq(:);
    mtThetaIdx = mtFreq >= cfg.thetaBand(1) & mtFreq <= cfg.thetaBand(2);
    mtDeltaIdx = mtFreq >= cfg.deltaBand(1) & mtFreq <= cfg.deltaBand(2);
    assert(any(mtThetaIdx), 'No multitaper frequencies fall inside thetaBand.');
    assert(any(mtDeltaIdx), 'No multitaper frequencies fall inside deltaBand.');

    thetaPowerMT = mean(double(S(:, mtThetaIdx)), 2);
    deltaPowerMT = mean(double(S(:, mtDeltaIdx)), 2);
    thetaDeltaDbMT = 10 .* log10(thetaPowerMT ./ (deltaPowerMT + eps));

    multitaper.time = time(1) + mtTimeRelative(:) - 1 ./ cfg.analysisFs;
    multitaper.freq = mtFreq;
    multitaper.thetaPower = thetaPowerMT;
    multitaper.deltaPower = deltaPowerMT;
    multitaper.thetaDeltaDb = thetaDeltaDbMT;
    multitaper.params = mtParams;
catch ME
    warning('detect_nonactive_theta_epochs:MultitaperFailed', ...
        'Multitaper diagnostic failed: %s', ME.message);
    multitaper.error = ME.message;
end
end

%% Diagnostic figures

function figurePaths = make_diagnostic_figures(result, waveletPower, cfg)
figurePaths = struct();

outputDir = fileparts(cfg.outputPath);
if isempty(outputDir)
    outputDir = pwd;
end
if cfg.saveFigures && exist(outputDir, 'dir') ~= 7
    mkdir(outputDir);
end

[~, outputBase] = fileparts(cfg.outputPath);

figurePaths.state = make_state_figure(result, cfg, outputDir, outputBase);
figurePaths.histogram = make_histogram_figure(result, cfg, outputDir, outputBase);
figurePaths.durationHistogram = make_duration_histogram_figure( ...
    result, cfg, outputDir, outputBase);
figurePaths.waveletExample = make_wavelet_example_figure( ...
    result, waveletPower, cfg, outputDir, outputBase);
end

function pathOut = make_state_figure(result, cfg, outputDir, outputBase)
fig = figure('Color', 'w', 'Name', 'Non-active theta epochs');

ax(1) = subplot(3, 1, 1);
plot(result.time, result.speed, 'Color', [0.75 0.2 0.1], 'LineWidth', 0.75);
hold on;
plot(xlim, [cfg.speedThreshold cfg.speedThreshold], 'k--', 'LineWidth', 1);
ylabel('speed');
title(sprintf('Low-speed, low-theta detector: channel row %d', ...
    result.bestChannelRow));
box off;

ax(2) = subplot(3, 1, 2);
plot(result.time, result.thetaDeltaDb, 'k', 'LineWidth', 0.75);
hold on;
plot(xlim, [result.thresholdDb result.thresholdDb], 'r--', 'LineWidth', 1);
ylabel('theta/delta (dB)');
box off;

ax(3) = subplot(3, 1, 3);
plot(result.time, double(result.candidateMask), 'Color', [0.15 0.35 0.8], ...
    'LineWidth', 0.75);
ylim([-0.05 1.05]);
ylabel('candidate');
xlabel('time in LFP clock (s)');
box off;

for k = 1:numel(ax)
    shade_intervals(ax(k), result.nonActiveIntervals);
end
linkaxes(ax, 'x');

pathOut = save_figure_if_requested(fig, outputDir, ...
    [outputBase '_state.png'], cfg.saveFigures);
end

function pathOut = make_histogram_figure(result, cfg, outputDir, outputBase)
fig = figure('Color', 'w', 'Name', 'Low-speed theta/delta histogram');
values = result.thetaDeltaDb(result.lowSpeedMask);
values = values(isfinite(values));

histogram(values, 60, 'FaceColor', [0.35 0.35 0.35], 'EdgeColor', 'none');
hold on;
yl = ylim;
plot([result.thresholdDb result.thresholdDb], yl, 'r--', 'LineWidth', 1.5);
ylim(yl);
xlabel('theta/delta among low-speed samples (dB)');
ylabel('count');
title(sprintf('Lower %gth percentile threshold = %.3f dB', ...
    cfg.thresholdPercentile, result.thresholdDb));
box off;

pathOut = save_figure_if_requested(fig, outputDir, ...
    [outputBase '_histogram.png'], cfg.saveFigures);
end

function pathOut = make_duration_histogram_figure(result, cfg, outputDir, outputBase)
fig = figure('Color', 'w', 'Name', 'Non-active interval durations');

if isempty(result.nonActiveIntervals)
    text(0.1, 0.5, 'No non-active intervals detected.');
    axis off;
else
    durations = result.nonActiveIntervals(:, 2) - result.nonActiveIntervals(:, 1);
    histogram(durations, 'FaceColor', [0.1 0.45 0.35], 'EdgeColor', 'none');
    hold on;
    yl = ylim;
    plot([cfg.minEpochDuration cfg.minEpochDuration], yl, 'k--', 'LineWidth', 1.5);
    ylim(yl);
    xlabel('non-active interval duration (s)');
    ylabel('count');
    title(sprintf('%d intervals, median %.3f s', ...
        numel(durations), median_omitnan(durations)));
    box off;
end

pathOut = save_figure_if_requested(fig, outputDir, ...
    [outputBase '_duration_histogram.png'], cfg.saveFigures);
end

function pathOut = make_wavelet_example_figure( ...
    result, waveletPower, cfg, outputDir, outputBase)

trialNumber = choose_example_trial(result, cfg);
if isempty(trialNumber)
    pathOut = '';
    return;
end

trialStart = result.trials.startTime(trialNumber);
trialStop = result.trials.stopTime(trialNumber);
idx = result.time >= trialStart & result.time <= trialStop;

fig = figure('Color', 'w', 'Name', 'Example wavelet theta/delta trial');

if ~any(idx)
    text(0.1, 0.5, sprintf('Trial %d has no samples on the analysis grid.', ...
        trialNumber));
    axis off;
else
    localTime = result.time(idx) - trialStart;
    trialPosition = result.position(idx);
    tfForDisplay = 10 .* log10(double(waveletPower(:, idx)) + eps);
    nonActiveTrialMask = interval_mask(result.time, result.nonActiveIntervals);
    nonActiveTrialMask = nonActiveTrialMask(idx);

    ax(1) = subplot(4, 1, 1);
    imagesc(localTime, result.waveletFreq, tfForDisplay);
    axis xy;
    hold on;
    xl = xlim;
    plot(xl, [cfg.thetaBand(1) cfg.thetaBand(1)], 'w--', 'LineWidth', 1);
    plot(xl, [cfg.thetaBand(2) cfg.thetaBand(2)], 'w--', 'LineWidth', 1);
    ylim([min(result.waveletFreq) max(result.waveletFreq)]);
    ylabel('freq. (Hz)');
    directionLabel = trial_direction_label(result, trialNumber);
    if isempty(directionLabel)
        title(sprintf('Trial %d, condition %g', trialNumber, ...
            result.trials.condition(trialNumber)));
    else
        title(sprintf('Trial %d, condition %g, direction %s', trialNumber, ...
            result.trials.condition(trialNumber), directionLabel));
    end
    colorbar;

    ax(2) = subplot(4, 1, 2);
    plot(localTime, result.thetaDeltaDb(idx), 'k', 'LineWidth', 0.75);
    hold on;
    plot(xlim, [result.thresholdDb result.thresholdDb], 'r--', 'LineWidth', 1);
    ylabel('theta/delta (dB)');
    box off;

    ax(3) = subplot(4, 1, 3);
    plot(localTime, result.speed(idx), 'Color', [0.75 0.2 0.1], 'LineWidth', 0.75);
    hold on;
    plot(xlim, [cfg.speedThreshold cfg.speedThreshold], 'k--', 'LineWidth', 1);
    ylabel('speed');
    box off;

    ax(4) = subplot(4, 1, 4);
    plot(localTime, trialPosition, 'Color', [0.15 0.35 0.8], ...
        'LineWidth', 0.75);
    hold on;
    plot(localTime(nonActiveTrialMask), trialPosition(nonActiveTrialMask), '.', ...
        'Color', [0.1 0.65 0.35], 'MarkerSize', 8);
    legend({'trajectory', 'non-active'}, 'Location', 'best');
    ylabel('position');
    xlabel('time from trial start (s)');
    box off;
    linkaxes(ax, 'x');
end

pathOut = save_figure_if_requested(fig, outputDir, ...
    [outputBase '_wavelet_example.png'], cfg.saveFigures);
end

function trialNumber = choose_example_trial(result, cfg)
if ~isempty(cfg.exampleTrial)
    trialNumber = round(cfg.exampleTrial);
    return;
end

trialNumber = [];
selectedTrials = find(result.trials.selected(:).');
for trial = selectedTrials
    trialRange = [result.trials.startTime(trial) result.trials.stopTime(trial)];
    if intervals_overlap_any(trialRange, result.nonActiveIntervals)
        trialNumber = trial;
        return;
    end
end

if ~isempty(selectedTrials)
    trialNumber = selectedTrials(1);
end
end

function label = trial_direction_label(result, trialNumber)
label = '';

if ~isfield(result.trials, 'WB') || numel(result.trials.WB) < trialNumber
    return;
end

value = result.trials.WB;
if iscell(value)
    value = value{trialNumber};
elseif isa(value, 'string')
    value = value(trialNumber);
elseif ischar(value)
    if size(value, 1) >= trialNumber
        value = value(trialNumber, :);
    end
elseif isnumeric(value) || islogical(value)
    value = value(trialNumber);
else
    return;
end

if isa(value, 'string')
    label = char(value);
elseif ischar(value)
    label = strtrim(value);
elseif isnumeric(value) || islogical(value)
    label = num2str(value);
end
end

function tf = intervals_overlap_any(interval, intervals)
if isempty(intervals)
    tf = false;
else
    tf = any(intervals(:, 1) < interval(2) & intervals(:, 2) > interval(1));
end
end

function mask = interval_mask(time, intervals)
mask = false(size(time));
if isempty(intervals)
    return;
end

for k = 1:size(intervals, 1)
    mask = mask | (time >= intervals(k, 1) & time <= intervals(k, 2));
end
end

function shade_intervals(ax, intervals)
if isempty(intervals)
    return;
end

axes(ax);
yl = ylim(ax);
hold(ax, 'on');
for k = 1:size(intervals, 1)
    patch(ax, [intervals(k, 1) intervals(k, 2) intervals(k, 2) intervals(k, 1)], ...
        [yl(1) yl(1) yl(2) yl(2)], [0.2 0.7 0.45], ...
        'EdgeColor', 'none', 'FaceAlpha', 0.18);
end
ylim(ax, yl);
end

function pathOut = save_figure_if_requested(fig, outputDir, fileName, saveFigures)
pathOut = '';
if saveFigures
    pathOut = fullfile(outputDir, fileName);
    saveas(fig, pathOut);
end
end

%% Numeric helpers and validators

function value = median_omitnan(values)
value = percentile_omitnan(values, 50);
end

function value = percentile_omitnan(values, pct)
values = double(values(:));
values = values(isfinite(values));
assert(~isempty(values), 'Cannot compute percentile of an empty vector.');
values = sort(values);

pct = max(0, min(100, pct));
rank = 1 + (numel(values) - 1) .* pct ./ 100;
lo = floor(rank);
hi = ceil(rank);

if lo == hi
    value = values(lo);
else
    weight = rank - lo;
    value = values(lo) .* (1 - weight) + values(hi) .* weight;
end
end

function tf = is_text_scalar(x)
tf = ischar(x) || isa(x, 'string');
if tf && ~ischar(x)
    tf = isscalar(x);
end
end

function tf = is_positive_scalar(x)
tf = isnumeric(x) && isscalar(x) && isfinite(x) && x > 0;
end

function tf = is_positive_or_inf_scalar(x)
tf = isnumeric(x) && isscalar(x) && (isinf(x) || (isfinite(x) && x > 0));
end

function tf = is_nonnegative_scalar(x)
tf = isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0;
end

function tf = is_percent_scalar(x)
tf = isnumeric(x) && isscalar(x) && isfinite(x) && x >= 0 && x <= 100;
end

function tf = is_two_element_numeric(x)
tf = isnumeric(x) && isvector(x) && numel(x) == 2 && all(isfinite(x));
end

function tf = is_logical_scalar(x)
tf = (islogical(x) || isnumeric(x)) && isscalar(x);
end
