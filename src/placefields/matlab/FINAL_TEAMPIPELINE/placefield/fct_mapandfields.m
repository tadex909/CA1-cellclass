% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% Based on animal's 1D position and spikes trains, this function returns
% rate maps and place fields detected with their properties.

% INPUT
% xpos    = animal's position over time
% idx_spk = index of spikes times expressed in xpos samples
% prm     = parameters structure with 2 sub-structures
%   -> prm.rmap = rate map parameters
%       -> xbin        = set of spatials bins (e.g., xbin = linspace(0, 200, 100))
%       -> tbin        = set of times intervals in samples (first col = start index, second col = stop index, row = interval number)
%       -> freq        = sample frequency
%       -> ismooth     = half-width of the Gaussian filter in spatial bins
%       -> idcond_t    = vector to affect a condition to each time interval (tbin)
%       -> nb_cond     = user can define a number of conditions to get standardized matrix format
%       -> idx_rem     = index of spikes to remove (e.g., animal's velocity below threshold)
%       -> xbin_rem    = number of xbin to remove on both xbin edges
%   -> prm.pf = parameters for place fields detection
%       -> mtd         = method used to simulate spikes trains ('random', 'poisson' or 'circular_shift') 
%       -> nb_rep      = number of rate maps generated with simulated spikes trains
%       -> min_len     = place field min width (in spatial bins)
%       -> max_len     = place field max width (in spatial bins)
%       -> max_btw_pf  = maximum inter place field distance to merge 2 place fields (in spatial bins) 
%       -> max_ext_pf  = maximum edge extension of place field based on a second pval threshold ('pval_edge_pf')
%       -> pval_crit   = pval threshold for place field detection
%       -> pval_edge_pf= pval threshold to extend place field borders
%       -> pctstab     = percentage of stable laps with mean place field
%       -> stab        = stability threshold to consider stable laps 
%       -> frmax       = active cell criteria: peak firing rate has to be greater than this value
%       -> frmean      = active cell criteria: mean firing rate has to be greater than this value
%       -> pctspklap   = active cell criteria: percentage of laps in which the cell has to fire at least one spike 
% missing parameters will be automatically filled with default values


% OUTPUT
% rmap = ratemap structure with a set of matrix:
%    -> nbspk   = spikes count in each spatial bin 
%    -> dwell   = time spent by the animal in each spatial bin 
%    -> fr      = firing rate in each spatial bin 
% matrix subscripts give information about matrix format:
% _tx : row: time bins (t), col: position bins (x)
% _tx matrix format are averaged over a set of laps in a given condition (using prm.rmap.idcond_t) to get _cx matrix format
% _cx : row: condition (c), col: position bins (x)
% _s subscript refers to a smoothed matrix
%    -> feat = ratemap features with the following fields:
%        -> isactive     = returns 1 for active, 0 otherwise
%        -> stb_oddeven  = stability index consisting of spatial correlation between the mean of the odd/even laps 
%        -> stb_meancorr = stability index defined by the mean of the spatial correlations between the mean firing rate and firing rates of each lap
%        -> stb_allpairs = stability index determined with the mean spatial correlation between all combinations of laps
%        -> sparsity     = spartsity index
%        -> si           = spatial information determined on the mean vectors
%        -> si_meanlap   = mean of spatial informations determined on all trials
%        -> localstab    = local stability index assessed at each spatial bin

% pf = place field structure with a set of matrix:
%    -> pval_pf = place field p-value for each spatial bin
%    -> ispf    =  boolean vector of placefield presence/absence for each bin (0 for no placefield, 1 otherwise)
% matrix subscripts give information about matrix format:
% _tx matrix refers to place fields per lap
% _tx : row: time bins (t), col: position bins (x)
% _cx matrix refers to mean place fields per condition
% _cx : row: condition (c), col: position bins (x)
%    -> feat: structure of place fieds features for all conditions, each feature
%    vector is calculated for all place fields detected, and ordered from the place field with the highest
%    peak (first component) to the lowest (last component)
%       ->  ratein               = mean firing rate inside the place field
%       ->  rateout              = mean firing rate outside the place field
%       ->  outovin              = rateout/ratein ratio
%       ->  outovinc             = corrected rateout/ratein ratio
%       ->  size                 = place field width in spatial bins
%       ->  size2                = another measure of place field width in spatial bins
%       ->  height               = place field height
%       ->  ifrmax               = index of spatial bin with hightest firing rate for the place field of interest
%       ->  devcenter_index      = index to quantify how much a place field peak is shifted from its center
%       ->  lap_size_t           = vector of place fields width per lap
%       ->  lap_ifrmax_t         = index of spatial bin with place field peak for each lap
%       ->  lap_dev              = standard deviation of the place field peaks per lap
%       ->  lap_maxdev           = max deviation of the place field peaks per lap
%       ->  lap_devwithmean      = deviation of the place field peaks per lap compared to the mean place field
%       ->  lap_size             = mean width of place fields per lap
%       ->  lap_varsize          = variability of place field width per lap  
%       ->  lap_diffsizewithmean = difference between mean place field width and mean width per lap
%       ->  lap_pct              = percentage of laps with detected fields
% parameters are included in rmap and pf structures 
    

function [rmap, pf] = fct_mapandfields(xpos, idx_spk, prm)


rmap = fct_rmap_wrap(xpos, idx_spk, prm.rmap);
pf = fct_find_placefield(xpos, rmap, prm.pf);
