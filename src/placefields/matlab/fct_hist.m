% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% Given a matrix 'mat', this function will count how many elements of each
% row are located in the 'bins' intervals and does not take into account values outside the edges.
% Bins must be equally spaced
% As third argument, you can put 'n' to get a normalized count (%) or 'nc' to
% get a cumulative normalized count (%)


function count = fct_hist(mat, bins, varargin)


if abs(std(diff(bins))) > (1000*eps)
    error('Bins are not equally spaced.');
end

pas = bins(2) - bins(1);

bins = bins + 0.5*pas;
bins = [(bins(1) - pas) bins];
bins = bins - 2*eps;


if isempty(mat)
    sz = size(mat);
    if sz(1) ~= 0
        count = zeros(size(mat, 1), length(bins) - 2);
    else
        count = zeros(1, length(bins) - 2);
    end
    return
end

if iscolumn(mat)
    mat = [mat NaN(length(mat), 1)];
end

count = hist(mat', bins);
 
if ~isvector(mat)
    count = count';
end


count(:, 1) = [];
count(:, end) = [];

if any(strcmp(varargin, 'n'))
    count = (count ./ sum(count, 2))*100;
end

if any(strcmp(varargin, 'nc'))
    countc = cumsum(count, 2);
    count = bsxfun(@rdivide, countc , countc(:, end))*100; 
end










