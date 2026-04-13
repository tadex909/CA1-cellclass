% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% Based on a firing rate vector/matrix 'fr' and 'dwell' (row: lap, col: spatial bin), 
% this function will determine the spatial information (Skaggs et al.,
% 1996) for each lap.

function si = fct_spatial_info(fr, dwell)

p = bsxfun(@rdivide, dwell,  nansum(dwell, 2));
x = fr;
r = nanmean(fr, 2);
xr = bsxfun(@rdivide, x, r);

si = nansum(p .* xr .* log2(xr), 2);