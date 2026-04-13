% Author(s): Marti Geoffrey   
% Epsztein Lab 2019

% Given a matrix B with 0 and 1, this function will find for each row the
% index of the first 1 (istart) and the last 1 (istop) of the sequence of 1.
% if several sequences of 1 are present in a given row, the function will
% return the index of the first sequence

function [istart, istop] = fct_find_startstop(B)

% B = [0     1     0     1
%      1     1     1     0
%      0     0     1     1
%      0     0     0     1
%      0     0     0     0
%      1     0     0     0
%      1     1     0     0
%      1     1     1     1
%      0     1     1     0
%      0     1     1     1];
 
 
[nb_tbin, nb_xbin] = size(B);
A = [false(nb_tbin, 1) diff(B')'];

[~, istart] = nanmax(A, [], 2);
[~, istop] = nanmin(A, [], 2);
idx = istart == istop;
istart(idx) = NaN;
istop(idx) = NaN;
istop = istop - 1;
istop(istop == 0) = nb_xbin;
idx = nansum(B, 2) == nb_xbin;
istart(idx) = 1;
istop(idx) = nb_xbin;