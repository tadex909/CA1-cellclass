% Author(s): Marti Geoffrey
% Epsztein Lab 2017


function [ww, xx, ffreq, coi] = fct_wavelet_mat(N2, fs, frange, space)

% Sampling Interval
dt = 1 / fs;

% Wavelet Scales
scales = 1 ./ flip(frange);
nb_scales = length(scales);

% Non-dimensional frequency
w0 = 6;

% Fourier wavelength (mother)
fwavelen =(4*pi)/(w0 + sqrt(2 + w0^2));

% Fourier frequencies
ffreq = 1 ./ (fwavelen.*scales);

% Cone of influence
coi = fwavelen / sqrt(2);
coi = coi*dt*[1E-5,1:((N2+1)/2-1),fliplr((1:(N2/2-1))),1E-5];

isnorm = false; 
% When isnorm was true, the wavelet matrix was normalized  to compensate
% the high frequencies spread. It will no longer be set as true. Only the
% final matrix (CWT) can be normalized for vizualisation purpose.

switch space
    case 'freq'
        % Angular Frequencies (wk) array (from 0 to pi)
        wk = 1:fix(N2/2);
        wk = wk.*((2*pi)/(N2*dt));
        wk = [0, wk, -wk(fix((N2-1)/2):-1:1)];
        
        
        % Daughter wavelets in Fourier Space
        % 1st method with loop
        fdaughter = zeros(nb_scales, N2);
        
        for s = 1:nb_scales
            fwlt = fct_wavelet(scales(s)*wk, 'morlet', 'freq');
            fnorm = sqrt(2*pi*scales(s) / dt)/ (isnorm*sqrt(scales(s)) + (1 - isnorm));
            fdaughter(s,:) = fnorm*fwlt;
        end
        
        % 2nd method with bsxfun
        %        wks = bsxfun(@times, wk, scales');
        %        fwlt = fct_wavelet(wks, 'morlet', 'freq');
        %        fnorm = sqrt(2*pi*scales / dt) ./ (isnorm*sqrt(scales) + (1 - isnorm));
        %        fdaughter = bsxfun(@times, fwlt, fnorm');
        
        xx = wk;
        ww = fdaughter;
    case 'time'
        x = linspace(-(N2-1)*dt/2, (N2-1)*dt/2, N2);
        
        % Daughter wavelets in Time Space
        % 1st method with loop
        daughter = zeros(nb_scales, N2);
        for s = 1:nb_scales
            twlt = fct_wavelet(x / scales(s), 'morlet', 'time');
            fnorm = ((dt./scales(s)).^(0.5)) ./ (isnorm*sqrt(scales(s)) + (1 - isnorm));
            daughter(s,:) = fnorm*twlt;
        end
        
        
        % 2nd method with bsxfun
        %         xs = bsxfun(@rdivide, x, scales');
        %         fnorm = ((dt./scales).^(0.5)) ./ (isnorm*sqrt(scales) + (1 - isnorm));
        %         twlt = fct_wavelet(xs, 'morlet', 'time');
        %         daughter = bsxfun(@times, twlt, fnorm');
        
        xx = x;
        ww = daughter;       
end









