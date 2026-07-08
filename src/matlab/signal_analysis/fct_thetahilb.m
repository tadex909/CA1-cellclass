% -----------------------------
% Written by MARTI Geoffrey 
% 02/16 
% 10/17 % Add Theta Module as third output
% -----------------------------

% Find Theta Phase & Theta Frequence


function [theta_freq, theta_phase, theta_mod] = fct_thetahilb(theta_sig, freq, varargin)


if any(strcmp(varargin, 'smooth'))
    smooth_level = varargin{strcmp(varargin, 'smooth') + 1};
else
    smooth_level = 0;    
end


theta_cmp = hilbert(theta_sig);
theta_phase = angle(theta_cmp);
theta_phase_uwp = unwrap(theta_phase);
if smooth_level > 0
    theta_phase_uwp = fct_smoothgauss(theta_phase_uwp, smooth_level);
end

theta_freq = abs(diff(theta_phase_uwp)) / (2*pi*(1/freq));
theta_freq = [theta_freq(1) theta_freq];

theta_mod = abs(theta_cmp);




