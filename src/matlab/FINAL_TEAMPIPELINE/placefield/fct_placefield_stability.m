% Author(s): Marti Geoffrey  
% Epsztein Lab 2019

% Based on a firing rate map matrix 'fr_s_tx' (row: lap number, col: spatial bin) and a boolean vector of
% detected mean place fields 'ispf_x' for each spatial bin, this function will validate a given place field when 2 criteria are verified:
% 1) the cell is active for the considered rate map 
% 2) the place field is stable
% To assess the stability of the field, spatial correlations are calculated between the firing
% rate vector of each trial and the mean firing rate vector. The place field is validated if
% the spatial correlations are greater than 'prm.stab' for at least 'prm.pctstab' % of trials.

% INPUT
% fr_s_tx  = smoothed firing rate matrix (row = lap / temporal bin, col = spatial bin)
% ispf_x   = vector with 0 (no placefield) and 1 (placefield) for each spatial bin, several place fields can be present
% prm      = parameters structure 
%       -> pctstab     = placefield stability criterion: percentage of stable laps with mean place field
%       -> stab        = placefield stability criterion: threshold to consider stable laps 
%       -> frmax       = active cell criterion: peak firing rate has to be greater than this value
%       -> frmean      = active cell criterion: mean firing rate has to be greater than this value
%       -> pctspklap   = active cell criterion: percentage of laps in which the cell has to fire at least one spike 


% OUTPUT
% ispf_x   = boolean vector for presence (0) / absence (1) of place field in each spacial bin
% isgoodpf = vector of stability for each lap
% pctstb   = percentage of stable laps


function [ispf_x, isgoodpf, pctstb] = fct_placefield_stability(fr_s_tx, ispf_x, prm)



[nb_tbin, nb_xbin] = size(fr_s_tx);
if ~fct_isactive(fr_s_tx, prm)
    isgoodpf = [];
    pctstb = [];
    ispf_x = zeros(1, nb_xbin);
    return
end



mtd = 0;


switch mtd
    case 0 % faster
        ispflab_x = bwlabel(ispf_x);
        nb_pf = max(ispflab_x(:));
        
        
        fr_s_tx_rep = repmat(fr_s_tx, nb_pf, 1);
        
        meanfr_s_x = nanmean(fr_s_tx, 1);
        ispf_x_rep = repmat(ispflab_x, nb_pf, 1);
        pf_id = repmat(1:nb_pf, nb_xbin, 1)';
        
        meanfr_s_x_rep = repmat(meanfr_s_x, nb_pf, 1);
        meanfr_s_x_rep(ispf_x_rep ~= pf_id & ispf_x_rep ~= 0) = NaN;
        meanfr_s_x_rep = kron(meanfr_s_x_rep, ones(nb_tbin, 1));
        
        stab_pft = reshape(fct_pearson(fr_s_tx_rep, meanfr_s_x_rep), nb_tbin, nb_pf);
        
        pctstb = (sum(stab_pft > prm.stab, 1) / nb_tbin);
        isgoodpf =  pctstb > prm.pctstab;
        
        ispf_x(ismember(ispflab_x, find(~isgoodpf))) = 0;
    case 1
        
        iseq_pf = fct_find_seq(ispf_x);
        nb_pf = size(iseq_pf, 1);
        
        pctstb = zeros(nb_pf, 1);
        
        
        for p = 1:nb_pf
            vkeeppf = true(1, nb_xbin);
            idx = setxor(p, 1:nb_pf);
            vkeeppf(fct_itv2idx(iseq_pf(idx, :))) = false;
            vkeeppfmat = repmat(vkeeppf, nb_tbin, 1);
            
            fr_tmp = fr_s_tx;
            fr_tmp(~vkeeppfmat) = NaN;
            [~, stb_t] = fct_stability_index(fr_tmp, 'meancorr');
            
            pctstb(p) = (sum(stb_t > prm.stab) / nb_tbin);
            
        end
        
        isgoodpf = pctstb >= prm.pctstab;
end








