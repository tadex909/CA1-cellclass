function [SpkDuration, SpkAsymmetry, SpkPeakTrough] = WaveformFeatureVinca(samples, average, best, bests,smplrate)

waveform.average = average;
waveform.samples = samples;
waveform.best = best;
waveform.bests=bests;
for ii = 1:size(waveform.average,3)  % boucle sur les cellules

    % on récupère directement le canal optimal
%     Maxindex = waveform.best(ii);

     for jj = 1:size(waveform.samples,3)  % boucle sur les 50 spikes
        
        [p1x(jj), p2x(jj), trx(jj), p1tx(jj), p2tx(jj), trtx(jj)] = ...
            spikeinfo(waveform.bests(:,jj,ii), smplrate);

     end

%%ICI JE SUIS PAS DACCORD pour duration :
    dd=sort(p2tx-trtx);
    SpkDuration(ii)=mean(dd(5:end-5));
%    plot(best(:,ii))
%    hold on
%    plot(mean(p2tx(ii)*25000),mean(p2x(ii)),'o')
%    plot([mean(p2tx*25000) mean(trtx*25000)],[mean(p2x),mean(trx)])
%    hold off
    asym=sort((p1x-p2x)./(p1x+p2x));
    SpkAsymmetry(ii)=mean(asym(5:end-5));
    
    pt=sort(trtx-p1tx);
    SpkPeakTrough(ii)=mean(pt(5:end-5));
    
%     HalfWidth(ii)=mean(halfwidth);
    
    clear p1x p2x trx p1tx p2tx trtx  % nettoyage

end
