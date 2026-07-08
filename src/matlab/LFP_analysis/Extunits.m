function [SpikeTimes,Samples,timInSamples]= Extunits (data,T, threshold);
%[SpikeTimes,Samples,timInSamples]= extractunits (data,threshold);
% data: a vector with the membrane potential data
%T: the time vector corresponding to data (must be same size)
% the threshold to be crossed to be considered an action potential
% outputs
% SpikeTimes: in ms
% Samples :a # APs* 250 samples array of voltage values
% timInSamples: where the spikes occured expressed in samples
%
data=data(:);
T=T(:);


CC=find(data>threshold);
CC(find(diff(CC)==1)+1)=[];


Samples=zeros(250,length(CC));


for ind=1:length(CC)
    if (CC(ind)-50)>1 && (CC(ind)+199)<length(data)
    Samples(:,ind)=data(CC(ind)-50:CC(ind)+199);
%     SamplesTimes(:,ind)=T(CC(ind)-50:CC(ind)+199);

    end;
end;


[AA,BB]=max(Samples);
DD=CC'+BB-50;
SpikeTimes=T(DD);
timInSamples=DD;
for ind=1:length(DD)
    if (DD(ind)-100)>1 && (DD(ind)+199)<length(data)
    Samples(:,ind)=data(DD(ind)-100:DD(ind)+149);
%   Samples(:,ind)=data(DD(ind)-50:DD(ind)+199);

    end;
    
end;