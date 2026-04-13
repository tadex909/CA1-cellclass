% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% This function returns several features (especially the stability, spatial information, sparsity and local stability) of a rate map.

% INPUT
% fr_tx     = firing rate map matrix (row: time bin /lap number, col: spatial bin)
% dwell_tx  = time spent matrix (row: time bin /lap number, col: spatial bin)
% prm       = parameters structure to consider active cells
%       -> frmean       = mean firing rate threshold
%       -> frmax        = peak firing rate threshold
%       -> pctspklap    = % of trials in which the cells has to fire at least one spike

% OUTPUT
% feat = features structure with the following fields:
% -> isactive     = returns 1 for active, 0 otherwise
% -> stb_oddeven  = stability index consisting of spatial correlation between the mean of the odd/even laps 
% -> stb_meancorr = stability index defined by the mean of the spatial correlations between the mean firing rate and firing rates of each lap
% -> stb_allpairs = stability index determined with the mean spatial correlation between all combinations of laps
% -> sparsity     = spartsity index
% -> si           = spatial information determined on the mean vectors
% -> si_meanlap   = mean of spatial informations determined on all trials
% -> localstab    = local stability index assessed at each spatial bin



function feat = fct_rmap_features(fr_tx, dwell_tx, prm)


fr_t = nanmean(fr_tx, 1);
dwell_t = nanmean(dwell_tx, 1);

feat.fr = nanmean(fr_t);
feat.isactive = fct_isactive(fr_tx, prm);
feat.stb_oddeven = fct_stability_index(fr_tx, 'oddeven');
feat.stb_meancorr = fct_stability_index(fr_tx, 'meancorr');
feat.stb_allpairs = fct_stability_index(fr_tx, 'allpairs');
feat.sparsity = fct_sparsity_index(fr_t);
feat.si = fct_spatial_info(fr_t, dwell_t);
feat.si_meanlap = nanmean(fct_spatial_info(fr_tx, dwell_tx));
feat.localstab = fct_slidstab(fr_tx, 2);

