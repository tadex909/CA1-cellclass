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


%%%%ADAPTED for spyking circus%%%%%%%%%%%

function [SpkDuration,SpkAsymmetry,SpkPeakTrough] = WaveformFeature_SC(samples,smplrate)


waveform.samples = samples;

NBcells = length(waveform.samples(:,1))

for ii=1:NBcells

        [p1x,p2x,trx,p1tx,p2tx,trtx]=spikeinfo(waveform.samples(ii,:),smplrate);%calculate feaures on shanks with the biggest spike
        
    
    %compute mean excluding the 5 edge values
    dd=sort(p2tx-trtx);
    SpkDuration(ii)=dd; %mean(dd(5:end-5));
    
    asym=sort((p1x-p2x)./(p1x+p2x));
    SpkAsymmetry(ii)= asym;%mean(asym(5:end-5));
    
    pt=sort(trtx-p1tx);
    SpkPeakTrough(ii)= pt;%=mean(pt(5:end-5));
end
    
    clear p1 p2 tr p1t p2t trt p1x p2x trx p1tx p2tx trtx
    
end
