% -----------------------------
% Written by MARTI Geoffrey
% 07/15
% -----------------------------

function [spatial_info] = fct_spatial_info_old(firing_rate, time_spent)

% The spatial_info is the same if we delete the nanvalues or not!!
% ind_del = find(isnan(firing_rate));
% firing_rate(ind_del) = [];
% time_spent(ind_del) = [];


p = time_spent / nansum(time_spent);
x = firing_rate;
r = nanmean(firing_rate);


spatial_info = nansum(p .* (x / r) .* log2(x / r));