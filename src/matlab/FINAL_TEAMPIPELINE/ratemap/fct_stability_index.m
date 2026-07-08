% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% This function determines a stability index of the rate map matrix 'fr_tx' 
% (row: time bin (lap), col: position bin). 
% In second argument, you can put:
% 'oddeven' -> spatial correlation between the mean of the odd laps and the
% mean of the even laps
% 'meancorr'-> mean of the spatial correlations between the mean firing
% rate and all laps
% 'allpairs'-> mean spatial correlation between all combinations of laps
% (option by default)

function [stabindex, stabvec] = fct_stability_index(fr_tx, varargin)


if nargin == 1
    opt = 'allpairs';
else
    opt = varargin{1};
end

nb_tbin = size(fr_tx, 1);
stabvec = [];
stabindex = NaN;

if nb_tbin <= 1
    return
end

switch opt
    case 'oddeven'
        mat1 = nanmean(fr_tx(1:2:end, :), 1);
        mat2 = nanmean(fr_tx(2:2:end, :), 1);
        stabindex = fct_pearson(mat1, mat2);
    case 'meancorr'
        mat1 = fr_tx;
        mat2 = repmat(nanmean(fr_tx, 1), nb_tbin, 1);
        stabvec = fct_pearson(mat1, mat2);
        stabindex = nanmean(stabvec);
    case 'allpairs'
        allpairs = nchoosek(1:nb_tbin, 2);   
        mat1 = fr_tx(allpairs(:, 1), :);
        mat2 = fr_tx(allpairs(:, 2), :);
        stabvec = fct_pearson(mat1, mat2);
        stabindex = nanmean(stabvec);
end











