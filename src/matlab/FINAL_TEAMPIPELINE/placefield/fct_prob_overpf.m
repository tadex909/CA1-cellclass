% -----------------------------
% Written by MARTI Geoffrey
% 03/16
% -----------------------------
% This probability density function is only defined when an overlap occured
% between the two placefield. Otherwise, the probability is null.


function [proba_noover] =  fct_prob_overpf(size_pf1, size_pf2, nb_bin, nb_overlap)

% size_pf1 = 10;
% size_pf2 = 10;
% nb_bin = 100;
nb_comb1 = nb_bin - size_pf1 + 1;
nb_comb2 = nb_bin - size_pf2 + 1;
% nb_overlap = 10;


gen_pf_way1 = zeros(nb_comb1, nb_bin);
for k  = 1:nb_comb1
    gen_pf_way1(k, k:(k+size_pf1-1)) = 1;  
end

gen_pf_way2 = zeros(nb_comb2, nb_bin);
for k  = 1:nb_comb2
    gen_pf_way2(k, k:(k+size_pf2-1)) = 1;  
end

overlap = zeros(nb_comb1*nb_comb2, nb_bin);
for i =1:nb_comb1
    for j = 1:nb_comb2
        overlap(nb_comb2*(i-1) + j, :) = gen_pf_way1(i, :) + gen_pf_way2(j, :) - ones(1, nb_bin);  
    end
end


overlap(overlap < 0) = 0;
overlap_ok = find(sum(overlap') >= nb_overlap);

proba = length(overlap_ok) / (nb_comb1*nb_comb2); % at least "nb_overlap" overlaps

proba_noover = 1 - proba; % at most "nb_overlap" - 1 overlap(s)
% ind = (nb_overlap / (min(size_pf1, size_pf2)))*proba_noover;