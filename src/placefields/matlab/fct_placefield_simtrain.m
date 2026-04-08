% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% This function generates a fixed number of rate maps ('prm.nb_rep') with simulated spikes trains.

% INPUT
% xpos    = animal's position over time
% rmap    = rate map structure from "fct_rmap_wrap" function
% prm     = must be the same parameters used for the "fct_rmap_wrap" function with two other parameters:
%    -> mtd       = method used to simulate spike trains ('random' or 'poisson') 
%    -> nb_rep    = number of rate maps generated with simulated spikes trains

% OUTPUT
%   fr_s_txrep = 3D matrix corresponding to a set of simulated rate maps (row: time bins (t), col: position bins (x), depth: ratemap number)


function fr_s_txrep = fct_placefield_simtrain(xpos, rmap, prm)


xpos(prm.idx_rem) = NaN;
prm.xbin([1:prm.xbin_rem (end -(prm.xbin_rem - 1)):end]) = [];

nb_tbin = size(prm.tbin, 1);
nb_xbin = length(prm.xbin) - 1;
nb_rep = prm.nb_rep;


dwell_s_tx = rmap.dwell_s_tx;
nbspk_tx = rmap.nbspk_tx;
tbin = prm.tbin;
xbin = prm.xbin;
ismooth = prm.ismooth;
mtd = prm.mtd;

fr_s_txrep = zeros(nb_tbin, nb_xbin, nb_rep);


parfor t = 1:nb_tbin
%for t = 1:nb_tbin
    
    idx = tbin(t, 1):tbin(t, 2);
    idxnonan = idx(1) + find(~isnan(xpos(idx)) & xpos(idx) >= xbin(1) & xpos(idx) <= xbin(end)) - 1;
    N = length(idxnonan);

    nbspk = sum(nbspk_tx(t, :));
    
    if N == 0 % very special case where there is no spatial candidate for spikes
        N = 1; % to make it works
        nbspk = 0;
    end
    
    switch mtd
        case 'random'
            itmp_spk = randi(N, nb_rep, nbspk);
            xpos_randspk = xpos(idxnonan(itmp_spk));
            
            if nbspk == 1 % very special case where the argument of xpos is not a matrix
                %xpos_randspk = xpos_randspk';
            end
            
            nbspk_repx = fct_hist(xpos_randspk, xbin);
            
        case 'poisson'
            nbspk_trep = poissrnd(nbspk / N, [nb_rep N]);
            
            idx_xbin = fct_discretize(xpos(idxnonan), xbin, 'bin');
            
            cc=find(idx_xbin==146); %%%bug quand y'a un 146 dans le bordel
            if isempty(cc) ~= 1
            idx_xbin(cc)=145;    
            end
   
            idx_xbin = repmat(idx_xbin, nb_rep, 1);
            idx_rep = repmat((1:nb_rep)', 1, N);
            
            %a = idx_rep(:);
            %b = idx_xbin(:);
            %cc=find(b==146);
            %c = nbspk_trep(:);
            %n = max([a b]);
            
            nbspk_repx = accumarray([idx_rep(:) idx_xbin(:)], nbspk_trep(:), [nb_rep nb_xbin]);
    end
    
    nbspk_s_repx = fct_smoothgauss(nbspk_repx, ismooth);
    
    fr_repx = bsxfun(@rdivide, nbspk_s_repx, dwell_s_tx(t, :));
    fr_s_repx = fct_smoothgauss(fr_repx, ismooth);
    
    fr_s_txrep(t, :, :) = fr_s_repx';
end