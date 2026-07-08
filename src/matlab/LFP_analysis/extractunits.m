function [SpikeTimes,Samples]= Extunits (abfpath, abfname,channel, threshold);
%[SpikeTimes,Samples]= extractunits (abfpath, abfname,channel);
%abfpath: where to find the file (it's a string)
%abfname: the name of the file (it's a string)
%channel: the name of the channel to use to find the strings (it ust be a
%cell of string
% outputs
% SpikeTimes: a # APs* 250 samples array of voltage values
%
%

cd (abfpath)
[data,si]=abfload(abfname,'channels', channel ); %%%% edit chanel here
SF=1/si*1E6;
% SF=SF;
T=[0:si:(si)*length(data)];
T(1)=[];


CC=find(data>threshold);
CC((diff(CC)==1)+1)=[];
times=T(CC);

Samples=zeros(250,length(CC));
SamplesTimes=Samples;
toadd=zeros(length(CC),1);

for ind=1:length(CC)
    Samples(:,ind)=data(CC(ind)-50:CC(ind)+199);
%     SamplesTimes(:,ind)=T(CC(ind)-50:CC(ind)+199);
end;


[AA,BB]=max(Samples);
DD=CC'+BB-50;
SpikeTimes=T(DD);

for ind=1:length(DD)
    Samples(:,ind)=data(DD(ind)-50:DD(ind)+199);
%     SamplesTimes(:,ind)=T(DD(ind)-50:DD(ind)+199);
end;