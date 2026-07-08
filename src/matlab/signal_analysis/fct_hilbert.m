% -----------------------------
% Written by MARTI Geoffrey 
% 02/16 
% 10/17 % Add Theta Module as third output
% -----------------------------


function [vfreq, vphase, vmod] = fct_hilbert(v, freq, varargin)


if any(strcmp(varargin, 'smooth'))
    sl = varargin{strcmp(varargin, 'smooth') + 1};
else
    sl = 0;    
end


vhilb = hilbert(v);
vphase = angle(vhilb);
vphase_uwp = unwrap(vphase);
if sl > 0
    vphase_uwp = fct_smoothgauss(vphase_uwp, sl);
end

vfreq = abs(diff(vphase_uwp)) / (2*pi*(1/freq));
vfreq = [vfreq(1) vfreq];

vmod = abs(vhilb);




