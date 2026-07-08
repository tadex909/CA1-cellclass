%function [acg,acgt]=ACG(T,G,gsubset,smplrate,binsize,timewidth)
%
%smplrate in Hz
%binsize and timewidth in ms

function [acg,acgt]=ACG(T,G,gsubset,smplrate,binsize,timewidth)

bin=smplrate*binsize/1000;
halfbin=timewidth/binsize;

for ii=1:length(gsubset)

	
    T = double(T);
    G = double(G);
    
    [ccg,acgt]=CCG(T,G,bin,halfbin,smplrate,gsubset(ii),'count');
    bar(ccg)
    %[ccg,acgt]=CCG_buz(T,G);

	acg(:,ii)=ccg(:,1,1);

end


if nargout==0
	ncell=length(gsubset);
	for i=1:ncell
		subplot(ncell,1,i); bar(acgt,acg(:,i))
		xlim([-timewidth timewidth]);
		%axis off
	end
end