% Addpath to where the program to use are
 %addpath(genpath('C:\1_jerome\matlab\MATLAB_pp_190911'));
 % My note: in the original program the dispersion delta gave strange
 % values. This is why I used "circ_std" from matlab circular statistics to calculate it

 % first open field recording for parietal cortex then entorhinal cortex
[abfname, abfpath]=uigetfile('*.abf','Select an ABF file','Multiselect','off');
cd (abfpath)
FieldChan{1}='Field';% or FieldChan{1}='Imscclamp'depending on the data
[Field_dg,si]=abfload(abfname,'channels', FieldChan ); %%%% edit chanel here
FieldTs=[0:si:(si)*(length(Field_dg))];
FieldTs(1)=[];
SF=1/si*1E6;
T=[0:si:(si)*length(Field_dg)];
T(1)=[];



% design the filter to use here we filter between 20 and 100 Hz
d=fdesign.bandpass(1, 20, 100, 5000, 60, .5, 60, SF);%% Generation of a 
%butterworth bandpass filter 500-3000 Hz
% Specification: 'Fst1,Fp1,Fp2,Fst2,Ast1,Ap,Ast2'
%d.description: {'Filter Order';'First Stopband Frequency';'First Passband Frequency';'Second Passband Frequency';'Second Stopband Frequency';'First Stopband Attenuation (db)';...
%'Passband Ripple (dB)';'Second Stopband Attenuation (dB)'} SF is sampling frequency
% First Stopband Frequency'    'First Passband Frequency'    'Second Passband Frequency'    'Second Stopband Frequency'    'First Stopband Attenuation (dB)'    'Passband Ripple (dB)'    'Second Stopband Attenuation (dB)'';...

d2=design(d,'butter', 'MatchExactly', 'passband');%last argument is IMPORTANT!!!!!!!!!!!!!!
fvtool(d2); %To visualize the shape of filter

% filter the data
FiltField=filtfilt(d2.sosMatrix,d2.ScaleValues,Field_dg);% 

% plot the filtered data on top of original data to check if dentate
% spikes are well detected
figure
plot(FieldTs(1:round(end/10)),Field_dg(1:round(end/10)))% if want to plot in fragment so it is less heavy to handle
hold on
plot(FieldTs(1:round(end/10)),FiltField(1:round(end/10)),'r-')% to plot the filtered data on top of original data in segments.
hold on
plot(FieldTs(1:round(end/10)),Compensated(1:round(end/10)),'g-')

% or to plot ALL the original data
%figure
%plot(FieldTs_Short,Field_pc)% 
%hold on% to plot ALL filtered data on top of ALL original data
%plot(FieldTs_Short,FiltField,'r-')


% get time of dentate spikes 
DSThresh = round (3 * std(FiltField));% threshold is set to 3 time std

CC=find(FiltField<(-1)*DSThresh);
%CC((diff(CC)==1)+1)=[];
ind_dentate_spike = fct_dentate_spike_detect(CC, Field_dg);

% plot dentatespiketimes to check for good detection
vec_ratio = 0.2; % ratio of the amm signal we want to plot
mode = 1; % mode = 1 to plot the dentate spike detected as circles/mode = 2 to plot the dentate spike detected as vertical bars
 fct_plot_dentate_spike(FieldTs, Field_dg, ind_dentate_spike, vec_ratio, mode)
DSTimes = FieldTs(ind_dentate_spike);
DSFreq = (numel (DSTimes)/ ((si)*(length(Field_dg)))*1E6);% in Hz


%figure
%plot(FieldTs,Field_dg)% if want to plot in fragment so it is less heavy to handle
%hold on
%plot(DSTimes, Field_dg(CC), 'ok')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

%%
if 0
% low pass filter the EEG
FreqRange = [40];
[b,a]=cheby2 (2,20,FreqRange./SF/2,'low');
SlowfiltField=filtfilt(b,a,Field_dg);
% plot result
figure, plot (Field_dg)
hold on
plot (SlowfiltField,'g','LineWidth',2);

end

% get the phase of the oscillation
[FieldPhase]=getPhaseFromField (SlowfiltField, SF);

% verify phase
normField_dg = (Field_dg - min(Field_dg)) / ( max(Field_dg) - min(Field_dg) );
normFieldPhase = (FieldPhase - min(FieldPhase)) / ( max(FieldPhase) - min(FieldPhase) );
figure, plot (normField_dg)
hold on
plot (normFieldPhase,'g','LineWidth',2);

% get the phase of each spike
[PhaseBin, PhaseHisto,theta, rbar, delta, sygma]=PhaseHistSpikes (DSTimes, FieldTs, FieldPhase, 1);

% plot the histogram of spike phase on the top of mean oscilation cycle
figure
bar(rad2ang(PhaseBin), PhaseHisto,'k')
hold on
bar(rad2ang(PhaseBin)+360, PhaseHisto,'k')
plot(rad2ang(PhaseBin), (cos(PhaseBin)+1)/2,'r')
plot(rad2ang(PhaseBin)+360, (cos(PhaseBin)+1)/2,'r')

title ('phase histogram')
xlabel('Field Phase')
ylabel('normalized # of DS')

xlim([-180 540])
%%
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
