% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% Based on a firing rate matrix 'rate_tr_bin' (row: time, col: position), 
% the function returns the state 1 (active) when:  
% - the mean firing rate is greater than 'prm.frmean' Hz, 
% - the peak firing rate is greater than 'prm.frmax' Hz, 
% - and the cell fired at least one spike in 'prm.pctspklap' % of the trials.
% and 0 (not active) otherwise

function state = fct_isactive(rate_tr_bin, prm)

mrate_bin = nanmean(rate_tr_bin);
nb_tr = size(rate_tr_bin, 1);

crit1 = nanmean(mrate_bin) > prm.frmean;
crit2 = nanmax(mrate_bin) > prm.frmax;
crit3 = (length(find(nanmean(rate_tr_bin, 2) > 0)) / nb_tr) > prm.pctspklap; 

state = crit1 & crit2 & crit3;