% We use Hartigan's dip test to test the bimodality (or non unimodality) of dentate
% granule cells Vm
% set path to find analysis softwares
addpath(genpath(':/Users/jeromeepsztein/Documents/jerome/matlab'));% on JE inmed's Imac
% first open data
%% 
[abfname, abfpath]=uigetfile('*.abf','Select an ABF file','Multiselect','off');
cd (abfpath)
PatchChannel{1}='10Vmclamp';
[data,si]=abfload(abfname,'channels', PatchChannel ); %%%% edit chanel here
SF=1/si*1E6;
T=[0:si:(si)*length(data)];
T(1)=[];

%if want to compensate for slowly drifting baseline Vm

[Compensated]= compensate_drift(data,SF);

if 1
% Then remove spike if there are some 
SpikeThresh=8;
[NoSpEEG]=removeSpikes2 (Compensated,T,SF,SpikeThresh);
VCorr = mean(data(1:2000000,1));
NoSpEEG = NoSpEEG + VCorr;
end

if 0
%otherwise
VCorr = mean(data(1:2000000,1));
NoSpEEG = Compensated+VCorr;% if there is no spikes
end
% verify good spike removal

figure,plot (data)
hold on
plot (NoSpEEG,'r')

if 1
% low pass filter to avoid gamma biasing the result
FreqRange = [40];
[b,a]=cheby2 (2,20,FreqRange./SF/2,'low');
FiltNoSpEEG=filtfilt(b,a,NoSpEEG);
end

% If want to remove some part of the data such as obvious artifact or drift
TShortStart1 = 200 * 20000; % with the desired start time expressed in seconds
TShortEnd1 = 380 * 20000; % with the desired end time expressed in seconds
Comp_Short1 = (NoSpEEG(TShortStart1:TShortEnd1,1));

%TShortStart2 = 400 * 20000; % with the desired start time expressed in seconds
%TShortEnd2 = 480 * 20000; % with the desired end time expressed in seconds
%Comp_Short2 = (NoSpEEG(TShortStart2:TShortEnd2,1));

%Comp_Short = cat (1,Comp_Short1, Comp_Short2);
Comp_Short = Comp_Short1;
%Comp_Short = FiltNoSpEEG;% if we take the all trace
data_Short = (data(TShortStart1:TShortEnd1,1));

if 0
% to verify Comp_Short
%data_Short1 = (data(TShortStart1:TShortEnd1,1));
%data_Short2 = (data(TShortStart2:TShortEnd2,1));
data_Short = cat (1,data_Short1,data_Short2);
%DD = downsample(data_Short,20,1);
end
% verify good trace selection and/or compensation
figure, plot(data_Short);
hold on
plot (Comp_Short, 'g');

% to perform silverman test we need data
data = CompShort;
save '/Volumes/EqpEpsztein/Geoff/Cont/DA110914_c1_200_380'...
    data ...
    SF
% test skewness
[skew] = skewness(Comp_Short)


if 0
% calculate the probability density function of Comp_short using a Kernel
% function: pdComp_Short

% enter binning here
Bin = 0.2;% to get 300 points

pdComp_Short = fitdist(Comp_Short,'Kernel','BandWidth',1); %pdComp_Short = fitdist(Comp_Short,'Kernel','BandWidth',4);
x = -90:Bin:-30;
yComp_Short = pdf(pdComp_Short,x);

norm_yComp_Short = (yComp_Short - min(yComp_Short)) / ( max(yComp_Short) - min(yComp_Short) );

%perform the test on Comp_Shortpdf
nboot = 10000;
[dip, p_value, xlow,xup]=HartigansDipSignifTest(norm_yComp_Short,nboot);
end

if 1
% another way to calculate the probability density function of compShort
% without Kernel density estimate but based on real density
% enter binning here
Bin = 0.2;
x = -90:Bin:-30;
Binranges = -90:Bin:-30;  % does the bins
bincounts = histc(Comp_Short,Binranges,1) ; % Vm being the Vm you what to test. Dim is the direction of the array (column in 1, line is 2). Facultatif
pdf  =  bincounts/length(Comp_Short); %should give you the pdf
figure, plot(x,pdf,'g-','LineWidth',1)

%perform the test on Comp_Shortpdf
nboot = 10000;
[dip, p_value, xlow,xup]=HartigansDipSignifTest(pdf,nboot);

end

% plot the distribution
Bin = 0.2;
x = -90:Bin:-30;
anal_time = (length (Comp_Short))/20000
dist_Comp_Short = hist(Comp_Short,-90:0.2:-30);
norm_dist = (dist_Comp_Short - min(dist_Comp_Short)) / ( max(dist_Comp_Short) - min(dist_Comp_Short) );
figure, plot(x,norm_dist,'b-','LineWidth',1)
xlim([-90 -30])
title([abfname,'dip=',num2str(dip,3),'dur=', num2str(anal_time,3) 'sec','bin=',num2str(Bin,2), 'p=',num2str(p_value,3)])

if 0
% plot the histogram
figure,
hist(Comp_Short,-90:0.15:-20)
title([abfname,'Bin=',num2str(Bin,2),'dip=',num2str(dip,3), ', p=',num2str(p_value,3)])
xlim([-90 -20])
end