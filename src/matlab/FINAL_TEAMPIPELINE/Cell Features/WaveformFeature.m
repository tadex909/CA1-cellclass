%function [SpkDuration,SpkAsymmetry] = WaveformFeature(waveform,smplrate)
%
%waveform should be structure waveform.average waveform.samples (output of
%SpikeTimeWave function)
%
%spike duration aand asymmetry: See Royer et al. 2012 Nat Neurosci, suppl
%fig 7A

%%%detect the biggest waveform (channel with the biggest spike = peak to
%%%trough) on the average ones

%%%calculate the spike parameters on the samples ones

%%I add peak to trough SpkPeakTrough

function [SpkDuration,SpkAsymmetry,SpkPeakTrough] = WaveformFeature(samples,average,smplrate)

waveform.average = average;
waveform.samples = samples;


for ii = 1:size(waveform.average,3)%%%select cells
    
    for jj = 1:size(waveform.average,1)  %% select channels
        
        if isnan(sum(waveform.average(jj,:,ii)))
             p1(jj)=0;
             p2(jj)=0;
             tr(jj)=0;
             p1t(jj)=0;
             p2t(jj)=0;
             trt(jj)=0;
        else
             [p1(jj),p2(jj),tr(jj),p1t(jj),p2t(jj),trt(jj)]=spikeinfo(waveform.average(jj,:,ii),smplrate);
        end
    
    end
    
    [crap,Maxindex]=max([p1-tr]);%detect channel in which the spike with the highest amplitude
    
    for jj = 1:size(waveform.samples,3)
        
        [p1x(jj),p2x(jj),trx(jj),p1tx(jj),p2tx(jj),trtx(jj)]=spikeinfo(waveform.samples(Maxindex,:,jj,ii),smplrate);%calculate feaures on shanks with the biggest spike
        
    end
    
    %compute mean excluding the 5 edge values
    dd=sort(p2tx-trtx);
    SpkDuration(ii)=mean(dd(5:end-5));
    
    asym=sort((p1x-p2x)./(p1x+p2x));
    SpkAsymmetry(ii)=mean(asym(5:end-5));
    
    pt=sort(trtx-p1tx);
    SpkPeakTrough(ii)=mean(pt(5:end-5));
    
    clear p1 p2 tr p1t p2t trt p1x p2x trx p1tx p2tx trtx
    
end
