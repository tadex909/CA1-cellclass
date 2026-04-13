% Author(s): Marti Geoffrey
% Epsztein Lab 2019

function si = fct_spatial_info(fr, dwell)

p = bsxfun(@rdivide, dwell,  nansum(dwell, 2));
x = fr;
r = nanmean(fr, 2);
xr = bsxfun(@rdivide, x, r);

si = nansum(p .* xr .* log2(xr), 2);