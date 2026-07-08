% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% Based on animal's 1D position and continuous vector 'vm' of the same length (e.g., velocity, membrane potential), 
% this function returns a median/std map structure with median and std of the vector 'vm' over spatial
% and temporal bins.

% INPUT
% xpos    = animal's position over time
% vm      = continuous vector (same length as xpos)
% prm     = rate map parameters
%       -> xbin        = set of spatials bins (e.g., xbin = linspace(0, 200, 100))
%       -> tbin        = set of times intervals in samples (first col = start index, second col = stop index, row = interval number)
%       -> freq        = sample frequency
%       -> ismooth     = half-width of the Gaussian filter in spatial bins
%       -> idx_rem     = index of spikes to remove (e.g., animal's velocity below threshold)
%       -> xbin_rem    = number of xbin to remove on both xbin edges
% missing parameters will be automatically filled with default values


% OUTPUT
% rmap = map structure with a set of matrix:
%    -> vmed   = vector median in each spatial bin
%    -> vstd   = vector std in each spatial bin
% matrix subscripts give information about matrix format:
% _tx : row: time bins (t), col: position bins (x)
% _s subscript refers to a smoothed matrix


function rmap = fct_rmapc(xpos, vm, prm)


xpos(prm.idx_rem) = NaN;
prm.xbin([1:prm.xbin_rem (end -(prm.xbin_rem - 1)):end]) = [];
vm(prm.idx_rem) = NaN;

nb_tbin = size(prm.tbin, 1);
nb_xbin = length(prm.xbin) - 1;

rmap.vmed_tx = zeros(nb_tbin, nb_xbin);

for t = 1:nb_tbin
    idx = prm.tbin(t, 1):prm.tbin(t, 2);
    xposlap = xpos(idx);
    vlap = vm(idx);
    ibinpos = fct_discretize(xposlap, prm.xbin, 'bin');
    isvnan = isnan(ibinpos);
    ibinpos(isvnan) = [];
    vlap(isvnan) = [];
    rmap.vmed_tx(t, :) = accumarray(ibinpos', vlap, [nb_xbin 1], @nanmedian, NaN);
    rmap.vstd_tx(t, :) = accumarray(ibinpos', vlap, [nb_xbin 1], @nanstd, NaN);
end

rmap.vmed_s_tx = fct_smoothgauss(rmap.vmed_tx, prm.ismooth);
rmap.vstd_s_tx = fct_smoothgauss(rmap.vstd_tx, prm.ismooth);
