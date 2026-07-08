function  QuantifyUpAndDown05 
 %%% we use this to perform spectrogram of a given electrophysiological
 %%% signal note that 5 sec should be added at the time interval we want to
 %%% analyse to have a correspondance between the trace and the spectrogram

 %% Parameters for the spectrogram
params.Fs=20000;
params.fpass=[0 20];
params.err=[2 0.05];
params.trialave=1;


% first get the data to analyse:
% for Vm
if 1
[abfname, abfpath]=uigetfile('*.abf','Select an ABF file','Multiselect','off');
cd (abfpath)
PatchChannel{1}='IN 0';
[data,si]=abfload(abfname,'channels', PatchChannel ); %%%% edit chanel here
data = (data-100)/20;
SF=1/si*1E6;
T=[0:si:(si)*length(data)];
T(1)=[];
end

if 0
% et/ou for Field
[abfname, abfpath]=uigetfile('*.abf','Select an ABF file','Multiselect','off');
cd (abfpath)
FieldChan{1}='Field_pCT';
[Field,si]=abfload(abfname,'channels', FieldChan ); %%%% edit chanel here
FieldTs=[0:si:(si)*length(Field)];
FieldTs(1)=[];
SF=1/si*1E6;
T=[0:si:(si)*length(Field)];
T(1)=[];
end

%% have a look at data
figure, plot(data); % and/or figure, plot(Field);

%% if want to analyse only part of the data
if 1
TShortStart1 = 16 * 20000; % with the desired start time expressed in seconds
TShortEnd1 = 285 * 20000; % with the desired end time expressed in seconds
dataShort = (data(TShortStart1:TShortEnd1,1));
TShort=[0:si:(si)*(length(dataShort))-1];
T(1)=[];
T = TShort;
data = dataShort;
end

%% have a look at data
if 1
figure, plot(dataShort); % and/or figure, plot(Field);
end
%% if want to compensate for slowly drifting in the data
if 1
[Compensated]= compensate_drift(dataShort,SF);%or [Compensated]= compensate_drift(dataShort,SF); if only part of the data
figure, plot(Compensated);
end

%% remove spikes for time-frequency analysis if no compensation
 if 1
SpikeThresh=8;
[Compensated]= compensate_drift(dataShort,SF);
[SpikeTimes,Samples,timInSamples]= Extunits (dataShort,T, SpikeThresh);
[NoSpEEG]=removeSpikes2 (dataShort,T,SF, SpikeThresh);

% verify good spike removal
figure
plot (dataShort,'r-')
hold on
plot (NoSpEEG)
 end
 

%% remove spikes for time-frequency analysis after compensation
 if 0
SpikeThresh=8;
[Compensated]= compensate_drift(data,SF);
[SpikeTimes,Samples,timInSamples]= Extunits (Compensated,T, SpikeThresh);
[NoSpEEG]=removeSpikes2 (Compensated,T,SF, SpikeThresh);

% verify good spike removal
figure
plot (Compensated,'r-')
hold on
plot (NoSpEEG)
 end
 

%% If want to remove Vm portion corresponding to hyperpolarizing pulses
if 0
% get the Im signal to detect hyp steps
[abfname, abfpath]=uigetfile('*.abf','Select an ABF file','Multiselect','off');
cd (abfpath)
PatchChannel{1}='Imscclamp';
[Im,si]=abfload(abfname,'channels', PatchChannel ); %%%% edit chanel here
% if we analyze only part of the data
if 1
ImShort = (Im(TShortStart1:TShortEnd1,1));
Im = ImShort;
end
% set a threshold to detect the steps
threshold = -20; % here 20 pA
% set a marge after the step in samples tp be sure to remove all of it
ind_marge = 400;% here 200 sample so 20 ms at 20KHz

[Field_nh ind_remove] = fct_remove_hyp_step(NoSpEEG, Im, threshold, ind_marge);
% to verify good removal
figure, plot (NoSpEEG)
hold on
plot (Field_nh, '-g')
plot (Im-20, '-r')
NoSpEEG = Field_nh;
end
%% Perform spectrogram of the FIELD and determine threshold for SO detection
 params.pad=1;
 params.tapers=[1 1];
 
WhitenNormalizeddata=locdetrend(NoSpEEG,SF,[0.3 0.1]);
[SpecIntra,t,f]=mtspecgramc(Field,[5 0.2],params);
figure;imagesc(t/60,f,SpecIntra',[0 3E5])% or [0 3E5] for extra

%% Perform spectrogram of the Vm and determine threshold for SO detection
 params.pad=1;
 params.tapers=[1 1];
 WhitenNormalizeddata=locdetrend(NoSpEEG,SF,[5 0.2]);
[SpecIntra,t,f]=mtspecgramc(WhitenNormalizeddata,[5 0.2],params);
figure;imagesc(t/60,f,SpecIntra',[0 20])% or [0 3E5] for extra

%% to plot locally detrended signal
figure, plot (NoSpEEG)
hold on
plot (WhitenNormalizeddata,'-r')
%% calculate power SWO band
PowerSOBand=sum(SpecIntra(:,f>0.1&f<=2),2);
figure;plot(PowerSOBand)
PowerSOBandSc=sum(SpecIntra(:,f>0.1&f<=2),2)/60;% not real value just to have it at the good scale for fig

MeanPowerSOBand = mean (PowerSOBand);


%% LowPowerSOBand = PowerSOBand (900:1200,1);
%figure, hist(PowerSOBand, 0:0.5:35);
ThresholdDetectionSO_90= prctile(PowerSOBand,90);
ThresholdDetectionSO_95= prctile(PowerSOBand,95);
figure, hist(PowerSOBand, 0:0.5:200);
ThresholdDetectionSO = 48.59436;
SOIndex=PowerSOBand>ThresholdDetectionSO;

%% to plot result
t = [1 1501];
figure, imagesc(t,f,SpecIntra',[0 20])% or [0 3E5] for extra
hold on
plot(PowerSOBandSc,'w')
hold on
plot(SOIndex, 'r')


MeanPowerSOBand = mean (PowerSOBand);

%% find long periods of slow oscillations
Beginnings=(find(diff(SOIndex)>=1))*0.2;% start time of SO in sec
Ends=(find(diff(SOIndex)<=-1))*0.2;% end time of SO in sec

%% check that beginnings and ends match !
PeriodsSO =Ends - Beginnings;
LongDetectedSO = PeriodsSO>4; % in sec
BegLongSO = Beginnings(LongDetectedSO);
EndLongSO = Ends(LongDetectedSO);
PeriodLongSO = EndLongSO - BegLongSO;
TotLongSO = sum (PeriodLongSO);
MeanLongSO = mean (PeriodLongSO);

%% save results
save '/Users/jeromeepsztein/Documents/jerome/Articles/1-Ouedraogo_et_al/eNeuro/results/Fig1/duration_power_SWO_period/threshold_90_percent_4s_with_spikes/Pilo/DA120726_c1_0008'...save '/Volumes/EqpCrepel/DAVID_2/2-in_vivo_analysis/Oscillation_Analysis/Spectrograms_and_Threshold_Analysis/SO_periods_intra/Cont/threshold_90_percent_4s_depol/DA110928_c2_0004'...
Beginnings...
Ends...
BegLongSO...
EndLongSO...
PeriodLongSO...
TotLongSO...
MeanLongSO...
MeanPowerSOBand



