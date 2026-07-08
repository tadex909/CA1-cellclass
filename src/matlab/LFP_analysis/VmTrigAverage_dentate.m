% Addpath to where the program to use are
 addpath(genpath('C:\1_jerome\matlab\MATLAB_pp_221012'));
% pour demarer les analyses, on peut faire:
[abfname, abfpath]=uigetfile('*.abf','Select an ABF file','Multiselect','off');
cd (abfpath)
PatchChannel{1}='10Vmclamp';
[data,si]=abfload(abfname,'channels', PatchChannel ); %%%% edit chanel here
SF=1/si*1E6;
T=[0:si:(si)*length(data)];
T(1)=[];

% if want to analyse only part of the data (eg a flat part for example)
% determine the start_time in ms *20 and end_time in ms *20  (*1000 in
% microsec/si)
TShort = T (1, 4E6:6E6);% we start by dimension 2 because T is in line
dataShort = data (4E6:6E6, 1);% we start by di

% if want to concatenate some data 
TShort = cat(2, T(0.97E6:1.45E6), T(3.18E6:4.52E6), T(5.9E6:7.25E6), T(8.54E6:9.88E6), T (1.14E7:1.277E7));% we start by dimension 2 because T is in line
dataShort = cat(1,data(0.97E6:1.45E6), data(3.18E6:4.52E6), data(5.9E6:7.25E6), data(8.54E6:9.88E6), data (1.14E7:1.277E7));% we start by dimension 1 because data is in column

%FieldTs_Short = FieldTs (1,7200000:12125184);

% to compensate for slow drift of the membrane potential
[Compensated]= compensate_drift(dataShort,SF);
plot (Compensated)
% to get time of the spikes
SpikeThresh=-40;

% to remove spikes
[NoSpEEG]=removeSpikes2 (dataShort,T,SF, SpikeThresh);
% to verify good spike removal
figure
plot (dataShort,'r-')
hold on
plot (NoSpEEG)

% if want to analyze only part of the data:
% The time interval is the time in ms * 20
NoSpEEGShort = NoSpEEG (4E6:6E6, 1);
figure, plot (NoSpEEGShort)

% another option is to detrend the data (e.g. remove the trend from the
% data)
    %detrend_data=detrend(data);
    %trend = data-detrend_data;

%hold on
%plot (T,data)
%plot (T,trend,':r')
%plot (T,detrend_data,'m')
%plot (T,zeros(size(T)),':k')
% Then remove spikes
    %SpikeThresh=20;
    %[SpikeTimes,Samples,timInSamples]= Extunits (detrend_data,T, SpikeThresh);
% to remove spikes

    %[NoSpEEG]=removeSpikes2 (Compensated,T,SF, SpikeThresh);


%to Z score
Z = zscore(Compensated);% or Z = zscore (NoSpEEG);%
figure
plot (Z);

% then start UD transition trig  average Vm per cell

[filename, path]=uigetfile ('*.txt', 'Select transition times');
cd (path)
[TTimes]=dlmread(filename);
TSamp = (TTimes*20);% note times should be expressed in samples not in ms (20 samples per msec)
TSamp = round(TSamp);
nBefore=6E4;
nAfter=6E4;
[Avs, StdErr] = TriggeredAvdentate(NoSpEEG, nBefore, nAfter, TSamp);

TAvs = (-3:50e-6:3);
TAvs = TAvs';
figure;
plot (TAvs,Avs);
save Y:\DAVID\2-in_vivo_analysis\PSTH\DGVm_vs_DUT\pilo\DA120704\DU_trig_Vm_Av_DA120704_c1_0011_0-300s_Zscored...
Z...
Avs ...
TAvs 


