function [WVMat]=DatWV (shank, AnlPath,SpikeTimes, config, header );
% to extract waveforms from dat file
%
%
%
%

minimunSpikes=50;
 SF=str2double(config.sampling_rate);
if length(SpikeTimes)>50
    NN=randperm(length(SpikeTimes),minimunSpikes);

    SpikeTimes=SpikeTimes(NN);
else
    SpikeTimes=SpikeTimes;    
end

FreqRange=[300, 5000];
[bb,aa]=fir1(3,FreqRange./(SF/2),'bandpass');

channels=[1:16; 17:32; 33:48; 49:64];
Ch2Keep=channels(shank,:);
cd (AnlPath)
AA=dir ('*.dat');
if ~isempty(AA)
    
    samples=AA(1).bytes/2/str2double(config.nb_channels);
    D=memmapfile(AA(1).name,'Format',{config.data_dtype [str2double(config.nb_channels) samples] 'mapped'});
else
    disp('cant find the dat file for waveform analysis')
    
end

try
    lfp=double(D.Data.mapped(Ch2Keep(1),1:end));
catch

    disp ('oops')
end

% recdur=length(lfp)/SF;
% filter the EEG
lfp=filtfilt(bb,aa,lfp);
%     lfp=double(D.Data.mapped(Channel,1:end));
% taking 2 ms before and 4 after
[WV, ~] = GetSegs(lfp, SpikeTimes-(0.002*SF), 0.006*SF);
MnWV=mean(WV,2);


WVMat=zeros(length(MnWV),length(Ch2Keep), minimunSpikes );
for ck=1: length(Ch2Keep)
    lfp=double(D.Data.mapped(Ch2Keep(ck),1:end));
    lfp=filtfilt(bb,aa,lfp');%probably useless in Rachel case since she 
    % taking 2 ms before and 4 after
    [WV,~] = GetSegs(lfp, SpikeTimes-(0.002*SF), 0.006*SF);
    
    
    for ww=1:size(WV,2) %normalize baseline
        bsln=mean(WV(1:0.001*SF, ww));
        WV(:,ww)=WV(:,ww)-bsln;
        
    end
    if length(SpikeTimes)<minimunSpikes
       disp('oups') 
       AddMe=zeros(size(WV,1),minimunSpikes-size(WV,2));
       WV=[WV, AddMe];
    end
    WVMat(:,ck,:)=(WV).*header.channels(1).bit_volts;
end


