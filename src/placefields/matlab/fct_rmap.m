% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% Based on animal's 1D position and spikes trains, this function returns a
% rate map structure with dwell, spike counts and firing rate over spatial
% and temporal bins.

% INPUT
% xpos    = animal's position over time
% idx_spk = index of spikes times expressed in xpos samples
% prm     = rate map parameters
%       -> xbin        = set of spatials bins (e.g., xbin = linspace(0, 200, 100))
%       -> tbin        = set of times intervals in samples (first col = start index, second col = stop index, row = interval number)
%       -> freq        = sample frequency
%       -> ismooth     = half-width of the Gaussian filter in spatial bins
%       -> idx_rem     = index of spikes to remove (e.g., animal's velocity below threshold)
%       -> xbin_rem    = number of xbin to remove on both xbin edges
% missing parameters will be automatically filled with default values


% OUTPUT
% rmap = ratemap structure with a set of matrix:
%    -> nbspk   = spikes count in each spatial bin 
%    -> dwell   = time spent by the animal in each spatial bin 
%    -> fr      = firing rate in each spatial bin 
% matrix subscripts give information about matrix format:
% _tx : row: time bins (t), col: position bins (x)
% _s subscript refers to a smoothed matrix


function rmap = fct_rmap(xpos, idx_spk, prm)


xpos(prm.idx_rem) = NaN;
prm.xbin([1:prm.xbin_rem (end -(prm.xbin_rem - 1)):end]) = [];
idx_spk(ismember(idx_spk, prm.idx_rem)) = [];

nb_tbin = size(prm.tbin, 1);
nb_xbin = length(prm.xbin) - 1;

rmap.nbspk_tx = zeros(nb_tbin, nb_xbin);
rmap.dwell_tx = zeros(nb_tbin, nb_xbin);

for t = 1:nb_tbin
    aa=find(idx_spk >= prm.tbin(t, 1) & idx_spk <= prm.tbin(t, 2));
    a=xpos(idx_spk(idx_spk >= prm.tbin(t, 1) & idx_spk <= prm.tbin(t, 2)));
    %b=fct_hist(a',prm.xbin);
    rmap.nbspk_tx(t, :) = fct_hist(xpos(idx_spk(idx_spk >= prm.tbin(t, 1) & idx_spk <= prm.tbin(t, 2)))', prm.xbin);
    rmap.dwell_tx(t, :) = (fct_hist(xpos(prm.tbin(t, 1):prm.tbin(t, 2))', prm.xbin) / prm.freq) + eps;
    a=sum(rmap.dwell_tx(t, :));
end

rmap.fr_tx = rmap.nbspk_tx ./ rmap.dwell_tx;


rmap.nbspk_s_tx = fct_smoothgauss(rmap.nbspk_tx, prm.ismooth);
rmap.nbspk_s_tx(rmap.nbspk_s_tx < 0) = 0;
rmap.dwell_s_tx = fct_smoothgauss(rmap.dwell_tx, prm.ismooth) + eps;
rmap.dwell_s_tx(rmap.dwell_s_tx < 0) = 0;
rmap.fr_s_tx = rmap.nbspk_s_tx ./ (rmap.dwell_s_tx);