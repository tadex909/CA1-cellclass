% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% Based on animal's 1D position and rate map structure from 'fct_rmap_wrap' function, this function returns
% place fields detected with their properties.

% INPUT
% xpos    = animal's position over time
% rmap    = rate map structure from 'fct_rmap_wrap' function
% prm     = parameters for place fields detection
%       -> mtd         = method used to simulate spike trains ('random', 'poisson' or 'circular_shift') 
%       -> nb_rep      = number of rate maps generated with simulated spike trains
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

function pf = fct_find_placefield(xpos, rmap, prm)

nb_tbin = size(rmap.prm.tbin, 1);
nb_xbin = length(rmap.prm.xbin) - 2*rmap.prm.xbin_rem - 1;
nb_cond = rmap.prm.nb_cond;


if nargin == 2
    prm = struct;
end

prm = dealin(prm);


rmap.prm.nb_rep = prm.nb_rep;
rmap.prm.mtd = prm.mtd;
fr_s_txrep = fct_placefield_simtrain(xpos, rmap, rmap.prm);


pf.pval_pf_tx = ones(nb_tbin, nb_xbin);
pf.ispf_tx = zeros(nb_tbin, nb_xbin);
for t = 1:nb_tbin
    pf.pval_pf_tx(t, :) = sum(bsxfun(@ge, squeeze(fr_s_txrep(t, :, :)), rmap.fr_s_tx(t, :)'), 2) / prm.nb_rep;
    pf.ispf_tx(t, :) = fct_placefield_threshold(pf.pval_pf_tx(t, :), prm);
end


pf.pval_pf_cx = ones(nb_cond, nb_xbin);
pf.ispf_cx = zeros(nb_cond, nb_xbin);
for c = 1:nb_cond
    idx = rmap.prm.idcond_t == c;
    
    A = squeeze(nanmean(fr_s_txrep(idx, :, :), 1));
    pf.pval_pf_cx(c, :) =  sum(bsxfun(@ge, A, rmap.fr_s_cx(c, :)'), 2) / prm.nb_rep;
    
    pf.ispf_cx(c, :) = fct_placefield_threshold(pf.pval_pf_cx(c, :), prm);
    pf.ispf_cx(c, :) = fct_placefield_stability(rmap.fr_s_tx(idx, :), pf.ispf_cx(c, :), prm);  
    pf.ispf_cx(c, :) = fct_placefield_peakorder(rmap.fr_s_cx(c, :), pf.ispf_cx(c, :));
    pf.ispf_tx(idx, :) = fct_placefield_keepoverlap(pf.ispf_tx(idx, :), pf.ispf_cx(c, :));
    pf.feat(c) = fct_placefield_features(rmap.fr_s_tx(idx, :), pf.ispf_tx(idx, :), pf.ispf_cx(c, :));
end
pf.prm = prm;

end


function prm = dealin(prm)

defprm = set_defprm;

defprm_names = fieldnames(defprm);
prm_names = fieldnames(prm);

idx = find(~ismember(defprm_names, prm_names));

for f = idx'
       warning(['Parameter ''' defprm_names{f} ''' is missing. Default value affected [' num2str(defprm.(defprm_names{f})) '].'])
       prm.(defprm_names{f}) = defprm.(defprm_names{f});   
end

end

% OLD
% function prm = set_defprm
% 
% prm.mtd = 'random';
% prm.nb_rep = 1000;
% 
% prm.min_len = 3;
% prm.max_len = 45;
% prm.max_btw_pf = 5;
% prm.max_ext_pf = 5;
% prm.pval_crit = 0.01;
% prm.pval_edge_pf = 0.30;
% 
% prm.pctstab = 0.40;
% prm.stab = 0.60;
% 
% prm.frmax = 1;%1.5
% prm.frmean = 0.3;%0.3
% prm.pctspklap = 0.50;
% 
% end

function prm = set_defprm

prm.mtd = 'random';
prm.nb_rep = 1000;

prm.min_len = 2;
prm.max_len = 45;
prm.max_btw_pf = 3; %old 5 
prm.max_ext_pf = 3; %old 5 
prm.pval_crit = 0.05; %old 0.01
prm.pval_edge_pf = 0.30;

prm.pctstab = 0.40;
prm.stab = 0.60;

prm.frmax = 1.5;%1.5
prm.frmean = 0.3;%0.3
prm.pctspklap = 0.45;

end
