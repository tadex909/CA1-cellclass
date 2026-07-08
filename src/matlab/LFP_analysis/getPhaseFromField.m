function [Phase]=getPhaseFromField (Field, SF)
%



     FreqRange = [0.1 5];
        FilterOrd = 2;
        Ripple = 20; 

        

        [b a] = cheby2(FilterOrd, Ripple, FreqRange/round(SF/2));
        Eegf = filtfilt(b,a,Field);

        %%%%% remove constant term to avoid bias
        Eegf = Eegf - mean(Eegf);
        Hilb = hilbert(Eegf);
        Phase = angle(Hilb);
