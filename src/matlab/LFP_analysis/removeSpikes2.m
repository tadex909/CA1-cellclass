function [NoSpEEG]=removeSpikes2 (data,T,SF, SpikeThresh);
%[NoSpEEG]=removeSpikes (EEG,T, SpikeThresh);
% 
% it takes 6ms (2ms before, 4 after) around spike and interpolates in between
% right now it does a line in between...

% I filter to detec spikes only
FreqRange = [300];
[b,a]=cheby2 (2,20,FreqRange./SF/2,'high');
FilteredForSpikes=filtfilt(b,a,data);



%SF=1/(T(2)-T(1))*1E6;
SamplesBefore=0.002*SF;
SamplesAfter=0.004*SF;
TotSamples=SamplesBefore+SamplesAfter;
   [SpikeTimes,Samples,timInSamples]= Extunits (FilteredForSpikes,T, SpikeThresh);
        
        SpikeAdresses=[timInSamples-SamplesBefore; timInSamples+SamplesAfter];
         if (SpikeAdresses(2,end))>length(data)
            SpikeAdresses(2,end)=length(data);
        end;
        ValuesAtAdresses=data(SpikeAdresses);
         LineBetween=zeros(TotSamples,length(SpikeTimes));
         ReplaceAdresses= LineBetween;
        for ind=1:length(SpikeTimes)
             LineBetween (:,ind)=linspace(ValuesAtAdresses(1,ind),ValuesAtAdresses(2,ind),TotSamples);
       
             ReplaceAdresses (:,ind)=(SpikeAdresses(1,ind):SpikeAdresses(2,ind)-1);
        end;
        NoSpEEG=data;
        NoSpEEG(ReplaceAdresses)=LineBetween;