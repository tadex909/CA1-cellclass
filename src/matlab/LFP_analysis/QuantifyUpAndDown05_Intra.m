function  QuantifyUpAndDown05 
 %%% we use this to perform spectrogram of a given electrophysiological
 %%% signal note that 5 sec should be added at the time interval we want to
 %%% analyse to have a correspondance between the trace and the spectrogram

 % Parameters for the spectrogram
params.Fs=20000;
params.fpass=[0 20];
params.err=[2 0.05];
params.trialave=1;

% first get the data to analyse:
[abfname, abfpath]=uigetfile('*.abf','Select an ABF file','Multiselect','off');
cd (abfpath)
PatchChannel{1}='10Vmclamp';
[data,si]=abfload(abfname,'channels', PatchChannel ); %%%% edit chanel here
SF=1/si*1E6;
T=[0:si:(si)*length(data)];
T(1)=[];
% et/ou
% FieldChan{1}='Field';
% [Field,si]=abfload(abfname,'channels', FieldChan ); %%%% edit chanel here
% FieldTs=[0:si:(si)*length(Field)];
% FieldTs(1)=[];
% SF=1/si*1E6;
% T=[0:si:(si)*length(Field)];
% T(1)=[];
[Compensated]= compensate_drift(data,SF);
figure, plot(Compensated);
% if want to analyse only part of the data
if 1
TShortStart1 = 1; % with the desired start time expressed in seconds
TShortEnd1 = 305 * 20000; % with the desired end time expressed in seconds
dataShort = (Compensated(TShortStart1:TShortEnd1,1));
TShort=[0:si:(si)*(length(dataShort))-1];
T(1)=[];
T = TShort;
data = dataShort;
end

%dataShort = data (4918000:5500920 ,1);

%Field_pc_Short = Field_pc (7200000:12125184,1);
%FieldTs_Short = FieldTs (1,7200000:12125184);

if 0
 %remove spikes for time-frequency analysis
SpikeThresh=8;
%[Compensated]= compensate_drift(data,SF);
[SpikeTimes,Samples,timInSamples]= Extunits (dataShort,T, SpikeThresh);
[NoSpEEG]=removeSpikes2 (dataShort,T,SF, SpikeThresh);
end

% verify good spike removal
figure
plot (dataShort,'r-')
hold on
plot (NoSpEEG)


if 0
% if there is no spike
[Compensated]= compensate_drift(data,SF);
[NoSpEEG]= Compensated;
end


% Perform spectrogram of the data and determine threshold for SO detection
 params.pad=1;
 params.tapers=[1 1];
 
% WhitenNormalizeddata=locdetrend(NoSpEEG,SF,[0.3 0.1]);
[SpecIntra,t,f]=mtspecgramc(dataShort,[0.3 0.2],params);% do the spectrogram on 5sec data length with 0.2 sec taper
PowerSOBand=sum(SpecIntra(:,f>0.1&f<=2),2);
PowerSOBandSc=sum(SpecIntra(:,f>0.1&f<=2),2)/60;% not real value just to have it at the good scale for fig
ThresholdDetectionSO = 48.59436;
SOIndex=PowerSOBand>ThresholdDetectionSO;
% to plot result
t = [1 1501];
figure
imagesc(t,f,SpecIntra',[0 3E5])% or [0 3E5] for extra
hold on
plot(PowerSOBandSc,'w')
hold on
plot(SOIndex, 'r')
MeanPowerSOBand = mean (PowerSOBand);
%LowPowerSOBand = PowerSOBand (900:1200,1);
%figure, hist(PowerSOBand, 0:0.5:35);
%ThresholdDetectionSO= prctile(PowerSOBand,90);
%ThresholdDetectionSO=mean(LowPowerSOBand)+std(LowPowerSOBand)*2;
ThresholdDetectionSO = 48.59436;



%save 'Y:\DAVID\2-in_vivo_analysis\Oscillation Analysis\UP_DOWN_states\Intra\DA120726_cell1_0008'...
%MeanPowerSOBand ...
%PowerSOBand ...
%ThresholdDetectionSO ...
%SpecIntra
    
% Find times where index is more important than threshold
SOIndex=PowerSOBand>ThresholdDetectionSO;
figure;plot(SOIndex)

% find long periods of slow oscillations
Beginnings=(find(diff(SOIndex)>=1))*0.2;% start time of SO in sec
Ends=(find(diff(SOIndex)<=-1))*0.2;% end time of SO in sec
PeriodsSO =Ends - Beginnings;
LongDetectedSO = PeriodsSO>4; % in sec
BegLongSO = Beginnings(LongDetectedSO);
EndLongSO = Ends(LongDetectedSO);
PeriodLongSO = EndLongSO - BegLongSO;
TotLongSO = sum (PeriodLongSO);
MeanLongSO = mean (PeriodLongSO);
save 'Y:\DAVID\2-in_vivo_analysis\Oscillation_Analysis\Threshold_Analysis\SO_periods_intra\Cont\threshold_90_percent_4s\DA120412_c1_0008' ...
ThresholdDetectionSO...
MeanPowerSOBand...
BegLongSO...
EndLongSO...
PeriodLongSO...
MeanLongSO ...
TotLongSO


