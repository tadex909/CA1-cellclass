% Author(s): Marti Geoffrey, Bourboulou Romain
% Epsztein Lab 2019

% This function smooths each row of the matrix 'data' using a convolution
% between the matrix and a Gaussian with a half-window size 'hwin' expressed in
% samples.


function data_s = fct_smoothgauss(data, hwin)

if iscolumn(data)
    data = data';
end
  
idxnan = isnan(data);
data_p = padarray(data, [0 hwin], 'replicate' , 'both');

win = (-hwin:hwin);
kernel = exp(-win.^2/(hwin/2)^2);
kernel = kernel / sum(kernel);

data_s = fct_conv(data_p, kernel);
data_s = data_s(:, (2*hwin) + 1:end-(2*hwin));
data_s(idxnan) = NaN;




