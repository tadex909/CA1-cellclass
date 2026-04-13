% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% Based on animal's 1D position and spikes trains, this function returns a
% rate map structure with dwell, spike counts and firing rate over spatial
% and temporal bins.

% INPUT
% xpos    = animal's position over time
% idx_spk = index of spikes times expressed in xpos samples
% prm     = rate map parameters
%       -> xbin        = set of spatial bins (e.g., xbin = linspace(0, 200, 100))
%       -> tbin        = set of times intervals in samples (first col = start index, second col = stop index, row = interval number)
%       -> freq        = sample frequency
%       -> ismooth     = half-width of the Gaussian filter in spatial bins
%       -> idcond_t    = vector to affect a condition to each time interval (tbin)
%       -> nb_cond     = user can define a number of conditions to get standardized matrix format (default is max(idcond_t))
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



function rmap = fct_rmap_wrap(xpos, idx_spk, prm)


if nargin == 2
    prm = struct;
end

prm = dealin(xpos, prm);
rmap = fct_rmap(xpos, idx_spk, prm);

% très important : prm.nb_cond peut être fixé pour avoir un vecteur de taille standard
nb_xbin = length(prm.xbin) - 2*prm.xbin_rem - 1;
nb_tbin = size(prm.tbin, 1);
nb_cond = prm.nb_cond;


mtd = 0;

switch mtd
    case 0 % Loop Version (faster)
        rmap.nbspk_cx = NaN(nb_cond, nb_xbin);
        rmap.dwell_cx = NaN(nb_cond, nb_xbin);
        rmap.fr_cx = NaN(nb_cond, nb_xbin);
        
        rmap.nbspk_s_cx = NaN(nb_cond, nb_xbin);
        rmap.dwell_s_cx = NaN(nb_cond, nb_xbin);
        rmap.fr_s_cx = NaN(nb_cond, nb_xbin);
        
        for c = 1:nb_cond
            idx = prm.idcond_t == c;
            
            rmap.nbspk_cx(c, :) = nanmean(rmap.nbspk_tx(idx, :), 1);
            rmap.dwell_cx(c, :) = nanmean(rmap.dwell_tx(idx, :), 1);
            rmap.fr_cx(c, :) = nanmean(rmap.fr_tx(idx, :), 1);
            
            rmap.nbspk_s_cx(c, :) = nanmean(rmap.nbspk_s_tx(idx, :), 1);
            rmap.dwell_s_cx(c, :) = nanmean(rmap.dwell_s_tx(idx, :), 1);
            rmap.fr_s_cx(c, :) = nanmean(rmap.fr_s_tx(idx, :), 1);
            rmap.feat(c) = fct_rmap_features(rmap.fr_s_tx(idx, :), rmap.dwell_s_tx(idx, :), prm);
        end
        
    case 1 % Accumarray Version
        prm.idcond_t(isnan(prm.idcond_t)) = nb_cond + 1;
        idx_c_rep = repmat(prm.idcond_t, 1, nb_xbin);
        idx_x_rep = repmat(1:nb_xbin, nb_tbin, 1);
        meanmapf = @(X) accumarray([idx_c_rep(:) idx_x_rep(:)], X(:), [nb_cond + 1 nb_xbin], @nanmean, NaN);
        remlast = @(X) X(1:end-1, :);
        meanmap = @(X) remlast(meanmapf(X));
        
        
        rmap.nbspk_cx = meanmap(rmap.nbspk_tx);
        rmap.dwell_cx = meanmap(rmap.dwell_tx);
        rmap.fr_cx = meanmap(rmap.fr_tx);
        
        rmap.nbspk_s_cx = meanmap(rmap.nbspk_s_tx);
        rmap.dwell_s_cx = meanmap(rmap.dwell_s_tx);
        rmap.fr_s_cx = meanmap(rmap.fr_s_tx);
        
        
        for c = 1:nb_cond
            idx = prm.idcond_t == c;
            rmap.feat(c) = fct_rmap_features(rmap.fr_s_tx(idx, :), rmap.dwell_s_tx(idx, :), prm);
        end
end

rmap.prm = prm;

end


function prm = dealin(xpos, prm)

defprm = set_defprm(xpos);

defprm_names = fieldnames(defprm);
prm_names = fieldnames(prm);

idx = find(~ismember(defprm_names, prm_names));

for f = idx'
       warning(['Parameter ''' defprm_names{f} ''' is missing. Default value affected [' num2str(defprm.(defprm_names{f})) '].'])
       prm.(defprm_names{f}) = defprm.(defprm_names{f});   
end

if ~isfield(prm, 'idcond_t')
    prm.idcond_t = ones(size(prm.tbin, 1), 1);
end

if ~isfield(prm, 'nb_cond')
    prm.nb_cond = nanmax(prm.idcond_t);
end

end




function prm = set_defprm(xpos)
prm.freq = 1000;
prm.ismooth = 7;

prm.xbin = linspace(floor(nanmin(xpos)), ceil(nanmax(xpos)), 100);
prm.tbin = [1 length(xpos)];

prm.xbin_rem = 10;
prm.idx_rem = [];
end





