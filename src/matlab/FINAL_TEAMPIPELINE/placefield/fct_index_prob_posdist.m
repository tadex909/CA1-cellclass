% -----------------------------
% Written by MARTI Geoffrey
% 03/16
% -----------------------------

% IL s'agit ici de comparer la quantité d'overlap entre les PF des deux sens 
% en distinguant deux overlaps :  si 2ème sens conservé, overlap de pos si 2ème sens inversé, overlpa de dist
% avec les overlaps de toutes les combinaisons de PF de même taille
% Si notre overlap est supérieur à l'overlap d'une combinaison, alors la
% probabilité que cet overlap ne soit pas arrivé par hasard augmente, sinon
% elle diminue. A noter que cette probabilité a une précision relative à la
% taille du PF sur l'échantillonage de l'espace, plus le PF est grand,
% moins les combinaisons sont possibles et donc plus l'overlap est
% important et plus la probabilité tend vers 0.



function [index_ppos, index_pdist] = fct_index_prob_posdist(ispf_bin_way1, ispf_bin_way2)


% V1 (with discontinuities)
% ispf_bin_way2_pos = ispf_bin_way2;
% ispf_bin_way2_dist = ispf_bin_way2(end:-1:1);
% 
% nb_bin = length(ispf_bin_way1);
% 
% 
% overlap_pos = ispf_bin_way1 + ispf_bin_way2_pos - ones(1, nb_bin);
% overlap_pos(overlap_pos < 0) = 0;
% overlap_dist = ispf_bin_way1 + ispf_bin_way2_dist - ones(1, nb_bin);
% overlap_dist(overlap_dist < 0) = 0;
% 
% nb_overlap_pos = sum(overlap_pos);
% nb_overlap_dist = sum(overlap_dist);
% 
% size_pf1 = length(find(ispf_bin_way1 == 1));
% size_pf2 = length(find(ispf_bin_way2_pos == 1));
% 
% index_ppos = fct_prob_overpf(size_pf1, size_pf2, nb_bin, nb_overlap_pos);
% index_pdist = fct_prob_overpf(size_pf1, size_pf2, nb_bin, nb_overlap_dist);
% 
% % Before: (proportions of overlapped field times the proba)
% % if any(overlap_pos)
% %     over_ratio = length(find(overlap_pos == 1)) / min(size_pf1, size_pf2);
% %     prob = fct_prob_overpf(size_pf1, size_pf2, nb_bin, nb_overlap_pos);
% %     index_ppos = over_ratio*prob;
% % else
% %     index_ppos = 0;
% % end

% V2 (proba continu)

ispf_bin_way2_pos = ispf_bin_way2;
ispf_bin_way2_dist = ispf_bin_way2(end:-1:1);

nb_bin = length(ispf_bin_way1);


nb_overlap_pos = fct_overlap_quantum(ispf_bin_way1, ispf_bin_way2_pos);
nb_overlap_dist = fct_overlap_quantum(ispf_bin_way1, ispf_bin_way2_dist);

size_pf1 = length(find(ispf_bin_way1 == 1));
size_pf2 = length(find(ispf_bin_way2_pos == 1));

index_ppos = fct_prob_overpf_cont(size_pf1, size_pf2, nb_bin, nb_overlap_pos);
index_pdist = fct_prob_overpf_cont(size_pf1, size_pf2, nb_bin, nb_overlap_dist);


