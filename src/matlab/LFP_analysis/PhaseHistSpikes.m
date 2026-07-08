function [PhaseBin, PhaseHisto,theta, rbar, delta, sygma]=PhaseHistSpikes (SpikeTimes, FieldTS, FieldPhase, normalize);
%[PhaseBin, PhaseHisto,theta, rbar, delta, sygma]=PhaseHistSpikes (SpikeTimes, FieldTS, FieldPhase, normalize);
% uses my circular stats package
% inputs:
%
% SpikeTimes: times od spikes
% FieldTS: should have the same size as Filed phase and units should fit
% those of SpikeTimes... they all should have the same start too...
% FieldPhase: a vector with the phase of the signal (use getPhasefromField
% normalize: put 1 if you want to normalize the histogram so that its max=1
%otherwise put 0
% outputs:
% PhaseBin : list of phases for the histogram (in radians)
% PhaseHisto : the histogram
% theta: mean phase
% rbar: railegh vector/resultant length
% delta: dispersion
% sygma: standard dev

 unitsPhase=interp1(FieldTS,FieldPhase,SpikeTimes,'nearest');
 
PhaseBin=[-pi:pi/18:pi]; %one degree resolution
 PhaseHisto=histc((unitsPhase),PhaseBin);
 
 
 [theta, rbar, delta, sygma] = circmeanPP(unitsPhase);

 % filter maybe?
%  [b a] = cheby2(FilterOrd, Ripple, 10/round(30/2));
%  PhaseHisto2=filtfilt(b,a, PhaseHisto);
%normalize? 
if normalize
    PhaseHisto=PhaseHisto./max(PhaseHisto);
end;