

function [D, P] = fct_distposindexpeaktopeak(A, B)


nb_fields = 5;
nb_bin = 80;

P = 1 - (nanmin(nanmin(abs(repmat(A, nb_fields, 1) - repmat(B, nb_fields, 1)')))) / nb_bin;
D = 1 - (nanmin(nanmin(abs(repmat(A, nb_fields, 1) - repmat(nb_bin - B, nb_fields, 1)')))) / nb_bin;

P = (P - 0.5)*2;
D = (D - 0.5)*2;