% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% This function determines the Pearson correlation coefficient between each
% row of two input matrix 'mat1' and 'mat2'.

function r = fct_pearson(mat1, mat2)


if iscolumn(mat1)
    mat1 = mat1';
end

if iscolumn(mat2)
    mat2 = mat2';
end

idx_rem = isnan(mat1) | isnan(mat2);


mat1(idx_rem) = NaN;
mat2(idx_rem) = NaN;

[~, vmat] = version;
vmat = str2double(vmat(end-3:end));

if vmat >= 2018
    mean1 = nanmean(mat1, 2);
    mean2 = nanmean(mat2, 2);
else
    % Difference between Matrix and Vector not allowed in previous Matlab
    % version
    nrep = size(mat1, 2);
    mean1 = repmat(nanmean(mat1, 2), 1, nrep);
    mean2 = repmat(nanmean(mat2, 2), 1, nrep);
end


cova = nansum(((mat1 - mean1) .* (mat2 - mean2)), 2);
sd1 = sqrt(nansum((mat1 - mean1).^2, 2));
sd2 = sqrt(nansum((mat2 - mean2).^2, 2));


r = cova ./ (sd1 .* sd2);