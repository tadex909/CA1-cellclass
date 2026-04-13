% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% Based on the firing rate map matrix 'fr_tx' (row: time bin /lap number,
% col: spatial bin) and a half-window 'hwin', this function will return a vector 
% 'slidstab' consisting of the local stability of the rate map at each
% spatial bin.
% The function will first determine two mean firing rate vectors for even and odd trials 
% in the neighborhood of each spatial bin (hwin spatial bins on the left and on the right and after). 
% Local stability index is then defined as the spatial correlation between 
% these two vectors for a given spatial bin.


function slidstab = fct_slidstab(fr_tx, hwin)

fr_tx_ode = [nanmean(fr_tx(1:2:end, :), 1) ; nanmean(fr_tx(2:2:end, :), 1)];
fr_tx_ode_p = padarray(fr_tx_ode, [0 hwin], 'replicate' , 'both');
[nb_tbin, nb_xbin] = size(fr_tx_ode_p);
slididx = bsxfun(@plus, 0:2*hwin, (1:(nb_xbin - (2*hwin)))');
fr_tx_slid = reshape(fr_tx_ode_p(:, slididx), nb_tbin* (nb_xbin - 2*hwin), 2*hwin + 1);
slidstab = fct_pearson(fr_tx_slid(1:2:end, :), fr_tx_slid(2:2:end, :))';


