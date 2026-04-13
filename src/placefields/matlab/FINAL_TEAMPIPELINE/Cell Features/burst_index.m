function [burstIndex]=burst_index(refT,G,acg)

%refT=refractory period

%for this index, it's ike Seb ones, except that I use the refractory period
%to detect the maximum proportion of spikes (burst)

%if the cell is very theta modulated the index is high...need to be
%modified = newburstindex_corrected

%for the background level I used the autocorrelogram with 10 ms bins (max
%between bin 1 and edge - edge is the 20 ms)

%%use the new refractory period (ACGrefractoryT_new)

%%%%THIS ONE IS THE GOOD ONE, OCTOBER 2014%%%%%



    %%%calculate new burst index
    burstIndex=NaN(length(G),1);
    for ii=1:length(G)
            
        if sum(acg.acg1(:,ii))>100 & isnan(refT(ii))==0
        
            %[maxacg,Imaxacg]=max(acg.acg1(51-(refT(ii)+5):50,ii));
            %[maxacg,Imaxacg]=max(acg.acg1(41:50,ii));
             maxacg=sum(acg.acg1(51-(refT(ii)+5):50,ii));
            
            
            %edge=mean(acg.acg1(1:acg.refractoryT_new(ii)+5,ii));
            %edge=mean(acg.acg1(1:10,ii));
            %edge=(max(acg.acg10(1:71-69,ii)))/10;
            edge=sum(acg.acg1(1:10,ii));
            
		
            amp=maxacg-edge;
		
            if amp>0
            	burstIndex(ii)=amp/maxacg;
            else
            	burstIndex(ii)=amp/edge;
            end
		
        end
	
    end




end

