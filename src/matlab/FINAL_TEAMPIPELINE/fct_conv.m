% Author(s): Marti Geoffrey, Bourboulou Romain
% Epsztein Lab 2019

% Convolution between the row vectors of two matrix 'x' and 'y' using Fourier
% transforms. Not that Strong but Fast.

function out = fct_conv(x, y)

N = size(x, 2) + size(y, 2) - 1;
N2 = pow2(nextpow2(N));
 
x(isnan(x)) = 0;
y(isnan(y)) = 0;

x_f = fft(x, N2 , 2); 
y_f = fft(y, N2 , 2);

out = ifft(bsxfun(@times , x_f, y_f), N2, 2);
out = out(:, 1:N);




