% Author(s): Marti Geoffrey   
% Epsztein Lab 2019

% Based on a p-value for place field detection for each spatial bin, this 
% function will first detect a set of contiguous values below a critical
% p-value ('prm.pval_crit') with a minimum length ('prm.min_len' in spatial bins), then merge
% these values if the distance between the closest edges is below 'prm.max_btw_pf'. 
% Place fields’ edges are extended by at most ('prm.max_ext_pf') spatial bins (for each edge) 
% when the p-value was below a more permissive threshold ('prm.pval_edge_pf'). 
% A field with a size greater than 'prm.max_len' spatial bins is excluded.

% INPUT
% pval_x  = p-value vector of place field bootstrap for each spatial bin
% prm     = parameters structure 
%       -> min_len     = place field min width (in spatial bins)
%       -> max_len     = place field max width (in spatial bins)
%       -> max_btw_pf  = maximum inter place field distance to merge 2 place fields (in spatial bins) 
%       -> max_ext_pf  = maximum edge extension of place field based on a second pval threshold ('pval_edge_pf')
%       -> pval_crit   = pval threshold for place field detection
%       -> pval_edge_pf= pval threshold to extend place field borders


% OUTPUT
% ispf_x = boolean vector for presence (0) / absence (1) of place field in
% each spacial bin
% iseq_pf = matrix of start/stop of detected place field(s)

function [ispf_x, iseq_pf] = fct_placefield_threshold(pval_x, prm)

iseq_pf = fct_find_seq(pval_x <= prm.pval_crit, 'minlen', prm.min_len);
iseq_pf = fct_merge_seq(iseq_pf, prm.max_btw_pf);
iseq_pf = fct_extent_edges_seq(iseq_pf, fct_find_seq(pval_x <= prm.pval_edge_pf), prm.max_ext_pf);
iseq_pf = fct_remove_seq(iseq_pf, 'maxlen', prm.max_len); 

nb_xbin = length(pval_x);
ispf_x = zeros(1, nb_xbin);
ispf_x(fct_itv2idx(iseq_pf)) = 1;