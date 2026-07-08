%function [peak1,peak2,trough,peak1t,peak2t,trought]=spikeinfo(waveform,smplrate)

%[min2,min2I]=min(waveform);
%[max1,max1I]=max(waveform(1:min2I));
%[max2,max2I]=max(waveform(min2I:end));
%max2I=max2I+min2I;
%[min1,min1I]=min(waveform(1:max1I));

%peak1=max1-min1;
%peak2=max2-min1;
%trough=min2-min1;

%peak1t=max1I/smplrate;
%peak2t=max2I/smplrate;
%trought=min2I/smplrate;

function [peak1,peak2,trough,peak1t,peak2t,trought]=spikeinfo_julie(waveform,smplrate)

[min2,min2I]=min(waveform); % trough....I for index (bin number)
[max1,max1I]=max(waveform(1:min2I)); %peak
% [max2,max2I]=max(waveform(min2I:end)); % second peak after the trough
max2Ix=find(waveform(min2I:end)>=max1)
if ~isempty(max2Ix)
max2I=max2Ix(1)+min2I - 1;  %%time to second peak
max2=waveform(max2I);
else
     [max2,max2I]=max(waveform(min2I:end)); 
     max2I=max2I+min2I - 1;
end
[min1,min1I]=min(waveform(1:max1I)); %%%first minimum beginning of waveform

peak1=max1-min1;
peak2=max2-min1;
trough=min2-min1;

%peak1=max1-min2;
%peak2=max2-min2;
%trough=min2-min1;

peak1t=max1I/smplrate;
peak2t=max2I/smplrate;
trought=min2I/smplrate;
% 
% figure(666)
% plot(waveform)
% hold on
% plot(max2I,max2,'o')
% plot(min2I,min2,'o')
% close(figure(666))
end
