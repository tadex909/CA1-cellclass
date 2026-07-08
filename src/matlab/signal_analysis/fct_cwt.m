% Author(s): Marti Geoffrey
% Epsztein Lab 2017

% Commented by Michon François-Xavier 2017/11
% [wt, ffreq, time, prm, wlt, coi] = fct_cwt(y, varargin)
%Inputs:
%   y --> signal vec
%   varargin:
%       'verbose'--> plot
%       'fs' --> signal fqz
%       'frange'--> [a : b] interest all fqz 
%       'method' --> computing method
%           'fprod' --> fqz space computing
%           'conv'--> temporal space computing
%       'norm' --> normalization
%Outputs:
%   wt --> continuous wavelet transform matrix (~spectrogram)
%   ffreq --> freq
%   time --> time
%   prm --> user parameters(struct) 
%   wlt --> Wavelet matrix
%   coi --> cone of influence 
% -----------------------------



function [wt, ffreq, time, prm, wltf, coi] = fct_cwt(y, varargin)

[prm, verbose] = debug_input(varargin);

% frange is actually the scale denominator, not exactly the final
% frequencies

% Check signal dimension
if size(y, 1) > 1
    y = y';
end

% Sample Length/Interval
N = length(y);
prm.dt = 1/ prm.fs;

% Choose wavelets length
switch prm.samplesize
    case 'pow2'
        N2 = 2^nextpow2(N);
        ispad = true;
    case 'classic'
        N2 = N;
        ispad = false;
        % Make Odd Sample Size
        %     if mod(N, 2) == 0
        %         y(end+1) = y(end);
        %         N = N + 1;
        %         N2 = N;
        %     end
end


% Daughter Wavelets in Time Space (convolutions) or Freq Space (Fourier rules)
switch prm.method
    case 'fprod'
        [wltf, ~, ffreq, coi] = fct_wavelet_mat(N2, prm.fs, prm.frange, 'freq');
    case 'conv'
        [wlt, ~, ffreq, coi] = fct_wavelet_mat(N2, prm.fs, prm.frange, 'time');
        wltf = fft(wlt, N2, 2);
        wltf = abs(wltf);
end


% Pad Data
if ispad
    padlen = floor((N2-N)/2);
    ypad = padarray(y, [0 padlen], 'circular', 'both');
    if length(ypad) < N2
        ypad(length(ypad)+1:N2) = 0;
    end
    ind = padlen+1:padlen+N;
else
    ypad = y;
    ind = 1:N;
end



% Fourier Transform of the Data
yf = fft(ypad, N2);


% Frequency Space Product between Data and Wavelett (FFT(Data) * FFT(Wavelett))
% 1st method with bsxfun
fprod = bsxfun(@times, yf, wltf);

% 2nd method with repmat
% yfrep = repmat(yf, nb_scales, 1);
% fprod = yfrep.*wltf;


wt = ifft(fprod,[],2);
wt = wt(:, ind);


time = 0:prm.dt:((N-prm.dt)/prm.fs);

if verbose
    if  prm.isnorm
        fct_plot_cwt(time, ffreq, abs(wt), 'norm')
    else
        fct_plot_cwt(time, ffreq, abs(wt))
    end
end


end
function [prm, verbose] = debug_input(X)

if any(strcmp(X, 'fs'))
    prm.fs = X{find(strcmp(X, 'fs')) + 1};
else
    prm.fs = 1000;
end

if any(strcmp(X, 'frange'))
    prm.frange = X{find(strcmp(X, 'frange')) + 1};
else
    prm.frange = 1:100;
end

if any(strcmp(X, 'method'))
    prm.method = X{find(strcmp(X, 'method')) + 1};
else
    prm.method = 'fprod';
end

if any(strcmp(X, 'pow2'))
    prm.samplesize = 'pow2';
else
    prm.samplesize = 'classic';
end

if any(strcmp(X, 'norm'))
    prm.isnorm = true;
else
    prm.isnorm = false;
end

if any(strcmp(X, 'verbose'))
    verbose = true;
else
    verbose = false;
end



end


