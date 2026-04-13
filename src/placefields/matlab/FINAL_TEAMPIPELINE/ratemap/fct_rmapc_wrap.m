% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% Based on animal's 1D position and continuous vector 'vm' of the same length (e.g., velocity, membrane potential), 
% this function returns a median/std map structure with median and std of the vector 'vm' over spatial
% and temporal bins.

% INPUT
% xpos    = animal's position over time
% vm      = continuous vector (same length as xpos)
% prm     = rate map parameters
%       -> xbin        = set of spatial bins (e.g., xbin = linspace(0, 200, 100))
%       -> tbin        = set of times intervals in samples (first col = start index, second col = stop index, row = interval number)
%       -> freq        = sample frequency
%       -> ismooth     = half-width of the Gaussian filter in spatial bins
%       -> idcond_t    = vector to affect a condition to each time interval (tbin)
%       -> nb_cond     = user can define a number of conditions to get standardized matrix format (default is max(idcond_t))
%       -> idx_rem     = index to remove (e.g., animal's velocity below threshold)
%       -> xbin_rem    = number of xbin to remove on both xbin edges
% missing parameters will be automatically filled with default values


% OUTPUT
% rmap = map structure with a set of matrix:
%    -> vmed   = vector median in each spatial bin
%    -> vstd   = vector std in each spatial bin
% matrix subscripts give information about matrix format:
% _tx : row: time bins (t), col: position bins (x)
% _tx matrix format are averaged over a set of laps in a given condition (using prm.rmap.idcond_t) to get _cx matrix format
% _cx : row: condition (c), col: position bins (x)
% _s subscript refers to a smoothed matrix


function rmap = fct_rmapc_wrap(xpos, vm, prm)


if nargin == 2
    prm = struct;
end

prm = dealin(xpos, prm);
rmap = fct_rmapc(xpos, vm, prm);

% très important : prm.nb_cond peut être fixé pour avoir un vecteur de taille standard
nb_xbin = length(prm.xbin) - 2*prm.xbin_rem - 1;
nb_cond = prm.nb_cond;


rmap.vmed_cx = NaN(nb_cond, nb_xbin);
rmap.vmed_s_cx = NaN(nb_cond, nb_xbin);


rmap.vstd_cx = NaN(nb_cond, nb_xbin);
rmap.vstd_s_cx = NaN(nb_cond, nb_xbin);


for c = 1:nb_cond
    idx = prm.idcond_t == c;
    
    rmap.vmed_cx(c, :) = nanmean(rmap.vmed_tx(idx, :), 1);
    rmap.vmed_s_cx(c, :) = nanmean(rmap.vmed_s_tx(idx, :), 1);
    
    rmap.vstd_cx(c, :) = nanmean(rmap.vstd_tx(idx, :), 1);
    rmap.vstd_s_cx(c, :) = nanmean(rmap.vstd_s_tx(idx, :), 1);
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
prm.freq = 100;
prm.ismooth = 10;

prm.xbin = linspace(floor(nanmin(xpos)), ceil(nanmax(xpos)), 100);
prm.tbin = [1 length(xpos)];

prm.xbin_rem = 0;
prm.idx_rem = [];
end





