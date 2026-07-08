% Addpath to where the program to use are
 %addpath(genpath('C:\1_jerome\matlab\MATLAB_pp_190911'));
 % My note: in the original program the dispersion delta gave strange
 % values. This is why I used "circ_std" from matlab circular statistics to calculate it

 % first open field recording for parietal cortex then entorhinal cortex
[abfname, abfpath]=uigetfile('*.abf','Select an ABF file','Multiselect','off');
cd (abfpath)
FieldChan{1}='Field_pCT';% or FieldChan{1}='Imscclamp'depending on the data
[Field_pc,si]=abfload(abfname,'channels', FieldChan ); %%%% edit chanel here
FieldTs=[0:si:(si)*(length(Field_pc))];
FieldTs(1)=[];
SF=1/si*1E6;
T=[0:si:(si)*length(Field_pc)];
T(1)=[];

% if want to analyse only part of the data (eg a flat part for example)
% determine the start_time in ms *200 and end_time in ms *200  (*1000 in
% microsec/si)

%TShort = T (1,7200000:12125184);
%Field_pc_Short = Field_pc (7200000:12125184,1);
%FieldTs_Short = FieldTs (1,7200000:12125184);
if 0
% get the times of the spikes if spike were detected in clampfit
[filename, path]=uigetfile ('*.txt', 'Select unit file');
cd (path)
[SpikeTimes_dg]=dlmread(filename);

SpikeTimes_dg = SpikeTimes_dg'*1000; % in microsec 

end

if 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to detect the Spikes in matlab
% design the filter to use
d=fdesign.bandpass(0.1, 500, 3000, 5000, 60, .5, 60, SF);%% Generation of a 
%butterworth bandpass filter 500-3000 Hz
% Specification: 'Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2'
%d.description: {'Filter Order';'First Stopband Frequency';'First Passband Frequency';'Second Passband Frequency';'Second Stopband Frequency';'First Stopband Attenuation (db)';...
%'Passband Ripple (dB)';'Second Stopband Attenuation (dB)'} SF is sampling frequency
% First Stopband Frequency'    'First Passband Frequency'    'Second Passband Frequency'    'Second Stopband Frequency'    'First Stopband Attenuation (dB)'    'Passband Ripple (dB)'    'Second Stopband Attenuation (dB)'';...

d2=design(d,'butter', 'MatchExactly', 'passband');%last argument is IMPORTANT!!!!!!!!!!!!!!
%fvtool(d2); %To visualize the shape of filter

% filter the data
FiltField=filtfilt(d2.sosMatrix,d2.ScaleValues,Field_pc);% 
% to plot the filtered data on to of the original data
figure
plot(FieldTs(1:round(end/10)),Field_pc(1:round(end/10)))% if want to plot in fragment so it is less heavy to handle
hold on
plot(FieldTs(1:round(end/10)),FiltField(1:round(end/10)),'r-')% to plot the filtered data on top of original data in segments.
% or to plot ALL the original data
%figure
%plot(FieldTs_Short,Field_pc)% 
%hold on% to plot ALL filtered data on top of ALL original data
%plot(FieldTs_Short,FiltField,'r-')


% get time of units using extunits
SpikeThresh = round (5 * std(FiltField));% threshold is set to 5 time std
 data=FiltField(:);

T=T(:);


CC=find(data>SpikeThresh);
CC((diff(CC)==1)+1)=[];

 SpikeTimes_pc=T(CC);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if 0
% then open field pc recording for parietal cortex then entorhinal cortex
[abfname, abfpath]=uigetfile('*.abf','Select an ABF file','Multiselect','off');
cd (abfpath)
FieldChan{1}='Field';% or FieldChan{1}='Imscclamp'depending on the data
[Field_pc,si]=abfload(abfname,'channels', FieldChan ); %%%% edit chanel here
FieldTs=[0:si:(si)*(length(Field_pc))];
FieldTs(1)=[];
SF=1/si*1E6;
T=[0:si:(si)*length(Field_pc)];
T(1)=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if 0
% low pass filter the EEG
FreqRange = [40];
[b,a]=cheby2 (2,20,FreqRange./SF/2,'low');
FiltField=filtfilt(b,a,Field_pc);
% plot result
figure, plot (Field_pc)
hold on
plot (FieldPhase,'g','LineWidth',2);

end

% get the phase of the oscillation
[FieldPhase]=getPhaseFromField (Field_pc, SF);

% verify phase
normField_pc = (Field_pc - min(Field_pc)) / ( max(Field_pc) - min(Field_pc) );
normFieldPhase = (FieldPhase - min(FieldPhase)) / ( max(FieldPhase) - min(FieldPhase) );
figure, plot (normField_pc)
hold on
plot (normFieldPhase,'g','LineWidth',2);

% get the phase of each spike
[PhaseBin, PhaseHisto,theta, rbar, delta, sygma]=PhaseHistSpikes (SpikeTimes_pc, FieldTs, FieldPhase, 1);

% plot the histogram of spike phase on the top of mean oscilation cycle
figure
bar(rad2ang(PhaseBin), PhaseHisto,'k')
hold on
bar(rad2ang(PhaseBin)+360, PhaseHisto,'k')
plot(rad2ang(PhaseBin), (cos(PhaseBin)+1)/2,'r')
plot(rad2ang(PhaseBin)+360, (cos(PhaseBin)+1)/2,'r')

title ('phase histogram')
xlabel('Field Phase')
ylabel('normalized # of spikes')

xlim([-180 540])

if 1
% to perform a Raileigh test on unit phase
 unitsPhase=interp1(FieldTs,FieldPhase,SpikeTimes_pc,'nearest');
 dataR = unitsPhase;% dataR = unitPhase' if spike are detected in clampfit and dataR = unitPhase if they are detetcted in matlab
 n = size(dataR, 1);
cols = size(dataR, 2);


Rbar = []; p = [];
for i=1:cols
  [theta, Rbar(i)] = circmean (dataR(:, i));
  Z = n*Rbar(i)^2;
  p(i) = exp(-Z) * (1 + (2*Z - Z^2) / (4*n) - (24*Z - 132*Z^2 + 76*Z^3 - 9*Z^4) / (288*n^2));
end
end
% verify  mean phase and phase dispersin and convert in deg
meanPhase = circ_mean (dataR);
meanDisp = circ_std (dataR);
meanPhaseDeg = rad2ang (meanPhase);
meanDispDeg = rad2ang (meanDisp);



% save results
save /Users/jeromeepsztein/Desktop/For_figS3/modDA120725_c1_0002_PCmua_vs_pcphase ...
PhaseBin ...
PhaseHisto ...
FieldTs ...
n ...
p ...
Rbar...% mean Rayleigh vector lenght
delta ...% dispersion
sygma ...% standard deviation
theta ... % mean phase (rad)
meanPhase ...% mean phase matlab (rad)
meanDisp ... % mean dispersion (rad)
meanPhaseDeg ... % mean phase (deg)
meanDispDeg  % mean dispersion (deg)Matlab
