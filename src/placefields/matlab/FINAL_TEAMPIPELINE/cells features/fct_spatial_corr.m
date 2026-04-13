% -----------------------------
% Written by MARTI Geoffrey
% 10/15
% -----------------------------

function [spatial_corr_dist, spatial_corr_pos] = fct_spatial_corr(firing_rate_way1, firing_rate_way2)


vec1 = firing_rate_way1;
vec2 = firing_rate_way2;

spatial_corr_pos = fct_pearson(vec1, vec2);
spatial_corr_dist = fct_pearson(vec1, vec2(end:-1:1));





