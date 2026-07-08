% -----------------------------
% Written by MARTI Geoffrey
% 15/04/15
% ----------------------------- 
% Mettre freq_down = NaN pour ne pas downsampler.
% ATTENTION : SI LA FREQ INITIALE NEST PAS DIVISIBLE PAR LA FREQ
% DOWNSAMPLEE, LA FREQUENCE OBTENUE NE SERA PAS EGALE A LA FREQ DOWNSAMPLEE
% CHOISIE

function [vec_d, time_d] = fct_dsampling(vec, freq, freq_d)

if isnan(freq_d)
    vec_d = vec;
    return
end

if mod(int64(freq), freq_d) == 0
    pas = floor(freq / freq_d);
    
    vec_d = vec(1:pas:end);
    time_d = 0:(length(vec_d) - 1);
    time_d = time_d / freq_d;
else
    warning('The freq is not a multiple of the down freq parameter. The signal will be interpolated.')
    time = 0:(length(vec)-1);
    time = time / freq;
    
    time_d = 0:(1/freq_d):time(end);
    vec_d = interp1(time, vec, time_d);    
end



    
    
    
    
       