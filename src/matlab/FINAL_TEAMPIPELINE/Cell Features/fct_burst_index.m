


function [burstIndex]= fct_burst_index(refT, acg)

%refT=refractory period

%for this index, it's ike Seb ones, except that I use the refractory period
%to detect the maximum proportion of spikes (burst)

%if the cell is very theta modulated the index is high...need to be
%modified = newburstindex_corrected

%for the background level I used the autocorrelogram with 10 ms bins (max
%between bin 1 and edge - edge is the 20 ms)

%%use the new refractory period (ACGrefractoryT_new)

%%%%THIS ONE IS THE GOOD ONE, OCTOBER 2014%%%%%



if sum(acg)> 100 && isnan(refT) == 0
    
    %[maxacg,Imaxacg]=max(acg.acg1(51-(refT(ii)+5):50,ii));
    %[maxacg,Imaxacg]=max(acg.acg1(41:50,ii));
    maxacg=sum(acg(51-(refT+5):50));
    
    
    %edge=mean(acg.acg1(1:acg.refractoryT_new(ii)+5,ii));
    %edge=mean(acg.acg1(1:10,ii));
    %edge=(max(acg.acg10(1:71-69,ii)))/10;
    edge=sum(acg(1:10));
    
    
    amp = maxacg-edge;
    
    if amp>0
        burstIndex = amp/maxacg;
    else
        burstIndex = amp/edge;
    end
else
    burstIndex = NaN;
end








