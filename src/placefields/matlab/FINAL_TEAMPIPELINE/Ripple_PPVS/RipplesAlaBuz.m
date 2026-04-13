function [ripples,maps,data,stats]= RipplesAlaBuz (RecPath, bestChan)

if nargin <1 % pour moi a modifier dans la ver finale
    
    %     RecPath=('D:\Rachel\RDM09\RDM09_2022-07-25_09-02-18\RecordNode103');
    %     bestChan=26;
    %
    RecPath=('D:\Rachel\RDM06\RDM06_2022-05-19_11-09-24\Record Node 103\Artefactfree');
    bestChan=3;
    %
end

% first get the best channel to analyze
% this requires that the recordings are ordered properly 1:64


cd (RecPath)


Name=dir('*.dat');
if ~isempty (Name)
      opath=pwd;
    %     cd .. % là aussi les dossiers où on trouve quoi c'est le bordel
    %     cd ..
    % tres pratique cette structure oebin
    fname = 'structure.oebin';
    fid = fopen(fname);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    val = jsondecode(str);
    header=val.continuous(1);
    
    % ya surement une meilleure façon mais je la trouve pas
    
    for ch=1:header.num_channels
        
        if (strcmp(header.channels(ch).channel_name,'ADC1'))
            chfound=ch;
        end
        
    end
    
    cd (opath)
     [lfp2, t2]=Read_OEP_Binary (Name.name,bestChan, 0,-1, 20);
     lfp2=lfp2';
     SF2=floor(1/(t2(2)-t2(1)));
     
elseif isfile('100_1.continuous')
      prefix='100_';
elseif isfile('100_CH1.continuous') 
      prefix='100_CH';
elseif isfile('101_1.continuous')
      prefix='101_';
end
    
  

  

if isempty (Name)% change this because rachel changed the filename 
Nam=[prefix,num2str(bestChan,'%d'),'.continuous'];
[lfp,timestamps, header] = load_open_ephys_data_faster(Nam);
SF=header.header.sampleRate;
%
% params.thresholds=[2 5];
% params.frequency=SF;
decFactor=SF/1250;
SF2=1250;
lfp2=downsample(lfp,decFactor);
t2=timestamps(1:decFactor:end);

% clear 'lfp', 'timestamps'
end

[ripples] = bz_FindRipples(lfp2,t2,'EMGThresh',[], 'durations' , [10 500],'passband',[90 300], 'show','off');
figure

plotind=1;
for rr=1:length(ripples.peaks)
    
    if mod(rr,25)==0
        figure;
        plotind=1;
    end
    subplot(5,5,plotind)
    index=find(t2>ripples.timestamps(rr,1)-0.1 & t2<ripples.timestamps(rr,2)+0.1);
    plot(t2(index), lfp2(index),'k')
    plotind=plotind+1;
end


wn2=[100, 500]./(SF2/2);
[b,a]=butter(3,wn2,'bandpass');

filtered=filtfilt(b,a,lfp2);

[maps,data,stats] = bz_RippleStats(filtered,t2,ripples);

[Mx,I]=max(data.peakAmplitude);
[Mn,Im]=min(data.peakAmplitude);
X=((1:size(maps.amplitude,2))./SF2)*1000;

timeforRipples=X;

save('RippleResults.mat','ripples', 'maps','data','stats','timeforRipples')


% figure
% subplot(4,2,1)
% plot(X, maps.ripples(I,:),'k')
% ylabel ('Ripple')
% title('Biggest Ripple')
% ylim([Mx*-1 Mx])
% subplot(4,2,3)
% plot(X, maps.frequency(I,:),'k')
% ylabel ('Frequency (Hz)')
% subplot(4,2,5)
% plot(X,maps.phase(I,:),'k')
% ylabel ('Phase (rads)')
% subplot(4,2,7)
% plot(X, maps.amplitude(I,:),'k')
% ylabel ('Aplitude (mV)')
% ylim([0 Mx+(Mx*0.1)])
% xlabel('Time (ms)')
%
% subplot(4,2,2)
% plot(X, maps.ripples(Im,:),'k')
% ylim([Mx*-1 Mx])
% ylabel ('Ripple')
% title('Smallest Ripple')
% subplot(4,2,4)
% plot(X, maps.frequency(Im,:),'k')
% ylabel ('Frequency (Hz)')
% subplot(4,2,6)
% plot(maps.phase(Im,:),'k')
% ylabel ('Phase (rads)')
% subplot(4,2,8)
% plot(X, maps.amplitude(Im,:),'k')
% xlabel('Time (ms)')
% ylabel ('Aplitude (\muV)')
% ylim([0 Mx+(Mx*0.1)])
%
%
%
% figure
% subplot(3,1,1)
% hist(data.peakFrequency,100)
% xlabel('Peak frequency (Hz)')
%
% subplot(3,1,2)
% hist(data.peakAmplitude,100)
% xlabel('Peak Amplitude (mV)')
% subplot(3,1,3)
% hist(data.duration*1000,100)
% xlabel('Duration (ms)')