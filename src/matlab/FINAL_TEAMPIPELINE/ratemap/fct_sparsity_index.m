% Author(s): Marti Geoffrey
% Epsztein Lab 2019

function sparsity_index = fct_sparsity_index(fr)

N = length(fr);
sparsity_index = (1 - (1/N)*(((nansum(fr))^2) / (nansum(fr .* fr))))*(N / (N - 1));
