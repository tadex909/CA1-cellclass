% Charge le signal CONTINOUS
% fpath='N:\Vinca\MATLAB\Data\Data_raw';
% [data, timestamps, info] = load_open_ephys_data('myfile.continuous');
clear all
savepath='\\Epsztein-nas02\TEAM\Tadeo\data\lfpdata';
%% DAT
sessname='VS103_2024-11-10_14-58-13';
dat_file = 'continuous.dat';
nChannels = 72;       % par ex. enregistrement 32 canaux
     % lire canal 5
fs = 25000;           % Hz
fs_new = 1250;        % Hz
allfp_ds=[];
channels=[];
channels_name=[];
for i= 1:64
    chan_to_read = i;
    [lfp_ds,channel,channel_name] = downsample_dat_fast(dat_file, nChannels, chan_to_read, fs, fs_new,strcat(sessname,'_lfp1250.mat'));
    allfp_ds=[allfp_ds;lfp_ds.'];
    channels=[channels;channel]
    channels_name=[channels_name;str2num(channel_name(3:end))]
end
save(strcat(savepath,'\',sessname,'_lfp1250Hz.mat'),'allfp_ds','channels','channels_name','-v7.3')
return
%% CONTINUOUS CLASSIQUE
fpath=('\\Epsztein-nas02\TEAM\Vinca\MATLAB\Data\Data_raw\VS11\VS11_2022-02-10_17-46-04');
sessname='VS11_2022-02-10_17-46-04'
allfp_ds=[];
channels=[];
channels_name=[];
for i = 1:64
fname = strcat(fpath,'\100_CH',string(i),'.continuous');
fs_orig = 25000;
fs_new = 1250;

block_skip = fs_orig / fs_new;  % ne lire qu'1 bloc sur 20

t_total = tic;
lfp_ds = read_continuous_downsample_block_fast(fname, block_skip);
   allfp_ds=[allfp_ds;lfp_ds.'];
    channels=[channels;i]
    channels_name=[channels_name;i];
end

save(strcat(savepath,'\',sessname,'_lfp1250Hz.mat'),'allfp_ds','channels', '-v7.3')

elapsed_total = toc(t_total);
fprintf('Total time: %.2f seconds (%.2f minutes)\n', ...
    elapsed_total, elapsed_total/60);

%% Filtered version    
fpath=('\\Epsztein-nas02\TEAM\Vinca\MATLAB\Data\Data_raw\VS47\VS47_2022-11-20_18-00-29');
sessname='VS47_2022-11-20_18-00-29'
allfp_ds=[];
channels=[];
channels_name=[];
for i = 1:64
fname = strcat(fpath,'\100_CH',string(i),'.continuous');
fs_orig = 25000;
fs_new = 1250;

block_skip = fs_orig / fs_new;  % ne lire qu'1 bloc sur 20

t_total = tic;
lfp_ds = read_continuous_resample_fast(fname, fs_orig, fs_new);
   allfp_ds=[allfp_ds;lfp_ds.'];
    channels=[channels;i]
    channels_name=[channels_name;i];
end

save(strcat(savepath,'\',sessname,'_lfp1250Hz.mat'),'allfp_ds','channels', '-v7.3')

elapsed_total = toc(t_total);
fprintf('Total time: %.2f seconds (%.2f minutes)\n', ...
    elapsed_total, elapsed_total/60);

%% CONTINUOUS 100_1
fpath=(cd);
sessname='VS61_2023-02-09_17-18-50'
allfp_ds=[];
channels=[];
channels_name=[];
for i = 1:64
fname = strcat(fpath,'\101_',string(i),'.continuous');
fs_orig = 25000;
fs_new = 1250;

block_skip = fs_orig / fs_new;  % ne lire qu'1 bloc sur 20

[lfp_ds , channel_name]= read_continuous_downsample_block_fast(fname, block_skip);
   allfp_ds=[allfp_ds;lfp_ds.'];
    channels=[channels;i]
    channels_name=[channels_name;channel_name];
end
save(strcat(savepath,'\',sessname,'_lfp1250Hz.mat'),'allfp_ds','channels','channels_name')

%%
function [lfp_ds,chan,chan_name] = downsample_dat_fast(dat_file, nChannels, chan_to_read, fs, fs_new, save_file)
    % Inputs
    % dat_file     : chemin du fichier .dat
    % nChannels    : nombre de canaux total
    % chan_to_read : index du canal � lire (1..nChannels)
    % fs           : fr�quence originale (ex. 25000)
    % fs_new       : fr�quence finale (ex. 1250)
    % save_file    : chemin pour sauvegarde (optionnel)

    % Facteur de d�cimation
    oebpath=pwd; 
            fname = 'structure.oebin';
            fid = fopen(fname);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            val = jsondecode(str);
            header=val.continuous(1);
            if startsWith(header.channels(chan_to_read).channel_name,'ADC')
                disp('skipping ADC')
                
                return
            end
    dec_factor = fs / fs_new;
    if mod(fs, fs_new) ~= 0
        error('Le ratio fs/fs_new doit �tre un entier (ici %f)', fs/fs_new);
    end

    % Ouvrir le fichier
    fid = fopen(dat_file, 'r');
    if fid == -1
        error('Impossible d''ouvrir le fichier %s', dat_file);
    end
    
    % On lit 1 �chantillon sur dec_factor
    % Chaque �chantillon = nChannels * int16 (2 bytes chacun)
    sample_size_bytes = nChannels * 2;
    skip_bytes = (dec_factor - 1) * sample_size_bytes;  % � sauter entre lectures

    % Positionner au canal voulu
    fseek(fid, (chan_to_read-1)*2, 'bof');

    % Lire directement un �chantillon sur dec_factor
    lfp_ds = fread(fid, Inf, 'int16', skip_bytes + (nChannels-1)*2);

    fclose(fid);

    % Convertir en double
    lfp_ds = double(lfp_ds);
    chan=chan_to_read;
    chan_name=header.channels(chan).channel_name
%     % Sauvegarde optionnelle
%     if nargin >= 6 && ~isempty(save_file)
%         save(save_file, 'lfp_ds', 'fs_new', 'chan_to_read');
%         disp(['Signal sauvegard� dans : ' save_file]);
%     end
end
function [lfp_ds , channel_name]= read_continuous_downsample_block_fast(fname, block_skip)
% Lire un fichier Open Ephys .continuous en ne prenant qu'1 bloc sur block_skip
% Inputs :
%   fname      : chemin complet du fichier .continuous
%   block_skip : nombre de blocs � sauter entre chaque bloc lu (ex: 20 = lire 1 bloc sur 20)
%
% Output :
%   lfp_ds : vecteur du signal downsampl�

fid = fopen(fname, 'r');
if fid == -1
    error('Impossible d''ouvrir le fichier %s', fname);
end

% Taille du fichier
fseek(fid, 0, 'eof');
file_size = ftell(fid);
fseek(fid, 0, 'bof');

% Lire header
hdr_bytes = 1024;
hdr = fread(fid, hdr_bytes, 'char*1');
eval(char(hdr'));
clear hdr
channel_name=str2num(header.channel(3:end))
lfp_ds = [];
block_count = 0;
record_marker_size = 10;

% Position initiale apr�s header
pos = hdr_bytes;
fseek(fid, pos, 'bof');

while pos < file_size
    block_count = block_count + 1;

    % Lire timestamp et nSamples
    timestamp = fread(fid, 1, 'int64', 0, 'l');
    nSamples = fread(fid, 1, 'uint16', 0, 'l');
    recordInd = fread(fid, 1, 'uint16');

    % Position apr�s en-t�te de record
    pos = ftell(fid);

    if mod(block_count-1, block_skip) == 0
        % Lire le bloc voulu
        data_block = fread(fid, nSamples, 'int16', 0, 'b');
        fread(fid, record_marker_size, 'char*1'); 
        lfp_ds = [lfp_ds; double(data_block)];
        pos = ftell(fid);
    else
        % Sauter le bloc non d�sir�
        skip_bytes = nSamples*2 + record_marker_size;
        if pos + skip_bytes > file_size
            skip_bytes = file_size - pos; % ajuster si dernier bloc
        end
        fseek(fid, skip_bytes, 'cof');
        pos = ftell(fid);
    end
end

fclose(fid);

% Conversion en microVolts
lfp_ds = lfp_ds * header.bitVolts;
end

function [lfp_ds, channel_name] = read_continuous_resample_fast(fname, fs_orig, fs_new)
% Read one Open Ephys .continuous file and downsample using MATLAB resample.
% resample applies an anti-aliasing filter before changing the sampling rate.
%
% Example:
%   [lfp_ds, channel_name] = read_continuous_resample_fast(fname, 25000, 1250);

if exist('resample', 'file') ~= 2
    error(['The function resample was not found. ', ...
        'Install Signal Processing Toolbox to use this reader.']);
end

[fid, msg] = fopen(fname, 'r');
if fid == -1
    error('Impossible d''ouvrir le fichier %s: %s', fname, msg);
end
cleanupObj = onCleanup(@() fclose(fid));

fseek(fid, 0, 'eof');
file_size = ftell(fid);
fseek(fid, 0, 'bof');

hdr_bytes = 1024;
hdr = fread(fid, hdr_bytes, 'char*1');
if numel(hdr) ~= hdr_bytes
    error('Header incomplet dans le fichier %s', fname);
end

eval(char(hdr'));
clear hdr

if isfield(header, 'channel')
    channel_token = regexp(header.channel, '\d+$', 'match', 'once');
    if isempty(channel_token)
        channel_name = NaN;
    else
        channel_name = str2double(channel_token);
    end
else
    channel_name = NaN;
end

if isfield(header, 'bitVolts')
    bit_volts = header.bitVolts;
else
    bit_volts = 1;
    warning('header.bitVolts not found in %s; returning raw int16 units.', fname);
end

records = {};
record_count = 0;
record_marker_size = 10;

while ftell(fid) < file_size
    timestamp = fread(fid, 1, 'int64', 0, 'l');
    nSamples = fread(fid, 1, 'uint16', 0, 'l');
    recordInd = fread(fid, 1, 'uint16', 0, 'l');

    if isempty(timestamp) || isempty(nSamples) || isempty(recordInd)
        break
    end

    data_block = fread(fid, nSamples, 'int16', 0, 'b');
    if numel(data_block) ~= nSamples
        warning('Dernier bloc incomplet dans %s; il sera ignore.', fname);
        break
    end

    marker = fread(fid, record_marker_size, 'char*1');
    if numel(marker) ~= record_marker_size
        warning('Marqueur de fin incomplet dans %s; arret de la lecture.', fname);
        break
    end

    record_count = record_count + 1;
    records{record_count, 1} = data_block; %#ok<AGROW>
end

if isempty(records)
    lfp_ds = [];
    warning('Aucun echantillon lu depuis %s.', fname);
    return
end

lfp = double(vertcat(records{:})) * bit_volts;
lfp_ds = resample(lfp, fs_new, fs_orig);
end

function [lfp_ds, channel_name] = read_continuous_firdecimator_block_fast(fname, fs_orig, fs_new, records_per_chunk)
% Read one Open Ephys .continuous file and downsample with dsp.FIRDecimator.
% This streams records in chunks, so it avoids loading the full channel into
% memory while keeping the FIR filter state between chunks.
%
% Example:
%   [lfp_ds, channel_name] = read_continuous_firdecimator_block_fast(fname, 25000, 1250);

if nargin < 4 || isempty(records_per_chunk)
    records_per_chunk = 128;
end

if records_per_chunk < 1 || records_per_chunk ~= floor(records_per_chunk)
    error('records_per_chunk must be a positive integer.');
end

if exist('dsp.FIRDecimator', 'class') ~= 8
    error(['The class dsp.FIRDecimator was not found. ', ...
        'Install DSP System Toolbox to use this reader.']);
end

if mod(fs_orig, fs_new) ~= 0
    error('Le ratio fs_orig/fs_new doit etre un entier (ici %f).', fs_orig/fs_new);
end

dec_factor = fs_orig / fs_new;
decim = dsp.FIRDecimator(dec_factor);

[fid, msg] = fopen(fname, 'r');
if fid == -1
    error('Impossible d''ouvrir le fichier %s: %s', fname, msg);
end
cleanupObj = onCleanup(@() fclose(fid));

fseek(fid, 0, 'eof');
file_size = ftell(fid);
fseek(fid, 0, 'bof');

hdr_bytes = 1024;
hdr = fread(fid, hdr_bytes, 'char*1');
if numel(hdr) ~= hdr_bytes
    error('Header incomplet dans le fichier %s', fname);
end

eval(char(hdr'));
clear hdr

if isfield(header, 'channel')
    channel_token = regexp(header.channel, '\d+$', 'match', 'once');
    if isempty(channel_token)
        channel_name = NaN;
    else
        channel_name = str2double(channel_token);
    end
else
    channel_name = NaN;
end

if isfield(header, 'bitVolts')
    bit_volts = header.bitVolts;
else
    bit_volts = 1;
    warning('header.bitVolts not found in %s; returning raw int16 units.', fname);
end

record_marker_size = 10;
record_buffer = cell(records_per_chunk, 1);
buffer_count = 0;
lfp_chunks = {};
chunk_count = 0;

while ftell(fid) < file_size
    timestamp = fread(fid, 1, 'int64', 0, 'l');
    nSamples = fread(fid, 1, 'uint16', 0, 'l');
    recordInd = fread(fid, 1, 'uint16', 0, 'l');

    if isempty(timestamp) || isempty(nSamples) || isempty(recordInd)
        break
    end

    data_block = fread(fid, nSamples, 'int16', 0, 'b');
    if numel(data_block) ~= nSamples
        warning('Dernier bloc incomplet dans %s; il sera ignore.', fname);
        break
    end

    marker = fread(fid, record_marker_size, 'char*1');
    if numel(marker) ~= record_marker_size
        warning('Marqueur de fin incomplet dans %s; arret de la lecture.', fname);
        break
    end

    buffer_count = buffer_count + 1;
    record_buffer{buffer_count} = data_block;

    if buffer_count == records_per_chunk
        x_chunk = double(vertcat(record_buffer{:})) * bit_volts;
        y_chunk = decim(x_chunk);
        if ~isempty(y_chunk)
            chunk_count = chunk_count + 1;
            lfp_chunks{chunk_count, 1} = y_chunk; %#ok<AGROW>
        end
        record_buffer(:) = {[]};
        buffer_count = 0;
    end
end

if buffer_count > 0
    x_chunk = double(vertcat(record_buffer{1:buffer_count})) * bit_volts;
    y_chunk = decim(x_chunk);
    if ~isempty(y_chunk)
        chunk_count = chunk_count + 1;
        lfp_chunks{chunk_count, 1} = y_chunk;
    end
end

if isempty(lfp_chunks)
    lfp_ds = [];
    warning('Aucun echantillon lu depuis %s.', fname);
else
    lfp_ds = vertcat(lfp_chunks{:});
end
end


