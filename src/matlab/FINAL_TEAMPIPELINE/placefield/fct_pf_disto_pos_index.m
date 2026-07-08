% -----------------------------
% Written by MARTI Geoffrey
% 07/15
% -----------------------------

function [dist_index, pos_index] = fct_pf_disto_pos_index(firing_rate_way1, pf_way1, firing_rate_way2, pf_way2)



% L'idée des disto et position index est qu'il faut un terme qui quantifie
% si une cellule code une position ou une vitesse. Pour cela, on définit
% pour un sens donné un place field en termes de position et en termes de
% distance. Ainsi définie, on regarde dans l'autre sens si ce PL a lieu ou
% non. S'il a lieu selon la distance, l'indice de distance augmente. S'il a
% lieu selon la position, l'indice de position augmente. On fait la même
% chose en définissant un PL distance et position à partir de l'autre sens.
% Exemple : 
% 3 bins
% 0cm   20cm    60cm
% firing_rate_way1 = [50 0 20]; sens : de 0cm vers 60cm
% firing_rate_way2 = [0 0 50]; sens : de 60cm vers 0cm
% pf_way1 = [1 0 0];
% pf_way2 = [0 0 1];

% Si on prend le sens 1, le PL est défini en termes de positions comme
% 0cm et en termes de distances comme 0cm après le début du track. 
% (vecteur [1 0 0])
% On regarde ainsi la décharge dans le deuxième sens. En termes de distance, on voit qu'à 0cm du
% début, c'est-à-dire à 60cm, il y a 50Hz qui viennent contribuer à
% l'augmentation de l'indice de distance. Il n'y a rien ailleurs, POUR CE SENS, l'indice
% de distance est donc maximal.
% En termes de position, on voit qu'à 0cm, il n'y a aucun spike, POUR CE SENS, l'indice
% de position est donc minimal.
% Si on se tient à cette formalisation, on voit que les indices ne
% quantifient pas exactement la position/distance. En effet, ils ne
% prennent pas en compte le 20 Hz dans le sens 1 qui devrait faire tendre
% l'indice de distance à diminuer et l'indice de position à augmenter. Il
% faut donc faire la même démarche dans le sens inverse. On définit un PL à
% partir du sens 2 comme étant 60 cm en termes de positions et 0cm du début
% du maze en termes de distance.  On regarde désormais dans le sens 1 et on
% voit qu'à 0cm du début, nous avons le terme 50 Hz qui vient augmenter
% l'indice de distance et le terme 20 qui vient le diminuer puisqu'il est
% hors du PL. Dans le sens 1, à la position 60cm, nous avons le terme 20 Hz
% qui vient augmenter l'indice de position.


pf_way1_pos = pf_way2;
pf_way1_dist = pf_way2(end:-1:1);

pf_way2_pos = pf_way1;
pf_way2_dist = pf_way1(end:-1:1);

mean_pos_in_field = nanmean(firing_rate_way2(pf_way2_pos == 1)) + nanmean(firing_rate_way1(pf_way1_pos == 1));
mean_pos_out_field = nanmean(firing_rate_way2(pf_way2_pos == 0)) + nanmean(firing_rate_way1(pf_way1_pos == 0));

pos_index = (mean_pos_in_field - mean_pos_out_field) ./ (mean_pos_in_field + mean_pos_out_field);

mean_dist_in_field = nanmean(firing_rate_way2(pf_way2_dist == 1)) + nanmean(firing_rate_way1(pf_way1_dist == 1));
mean_dist_out_field = nanmean(firing_rate_way2(pf_way2_dist == 0)) + nanmean(firing_rate_way1(pf_way1_dist == 0));


dist_index = (mean_dist_in_field - mean_dist_out_field) ./ (mean_dist_in_field + mean_dist_out_field);

