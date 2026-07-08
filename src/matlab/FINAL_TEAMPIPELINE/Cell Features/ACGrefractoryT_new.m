%function refractoryT=ACGrefractoryT(acg)
%
%acg bins should be in ms.
%%I change the snoothing. I went from 5 to 8



function refractoryT=ACGrefractoryT_new(acg)

if ~(sum(acg == 0) == length(acg))

startpoint=ceil(length(acg)/2)+1;

halfacg=acg(startpoint:end);
%figure(1)
%plot(halfacg)

%Shalfacg=smooth1D(halfacg,5);
%figure(2)
%plot(Shalfacg)

Shalfacg=smooth1D(halfacg,8);
%figure(3)
%plot(Shalfacg)

ndx=LocalMaxima(Shalfacg);

firstmode=ndx(1);

Dhalfacg=diff(halfacg(1:firstmode));

tmp=find(Dhalfacg>std(Dhalfacg));

if length(tmp)>0
	refractoryT=tmp(1)+1;
else
	refractoryT=nan; %if it's nan it's probably not a neuron
end

else
    refractoryT=nan;
end