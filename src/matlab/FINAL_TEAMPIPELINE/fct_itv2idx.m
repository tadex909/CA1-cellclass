% Author(s): Marti Geoffrey
% Epsztein Lab 2019

% INPUT
% iseq = matrix with several lines which consist of start (first col) and stop (second col) 

% OUTPUT
% idx_full = vector of all indexes between each start/stop intervals


function idx_full = fct_itv2idx(iseq, mtd)


if nargin == 1
    mtd = 1;
end

switch mtd
    case 1 % Good Tradeoff
        seqlen = diff(iseq') + 1;
        idx_full = zeros(1, sum(seqlen));
        idx = [0 cumsum(seqlen)];
        for k = 1:(length(idx) - 1)
            idx_full((idx(k) + 1):idx(k + 1)) = iseq(k, 1):iseq(k, 2);
        end
    case 2 % Good for Many Intervals, Weak for Big Index 
        tmp(iseq(:, 2) + 1) = -1;
        tmp(iseq(:, 1)) = 1;
        idx_full = find(cumsum(tmp));
    case 3 % Weak for Many Intervals, Good for Big Index
        idx_full = cell2mat(arrayfun(@(n) iseq(n, 1):iseq(n, 2), 1:size(iseq, 1), 'uni', 0));
    case 4
        iseq(:, 2) = iseq(:, 2) + 1;
        idx = fct_discretize(1:max(iseq(:)), iseq, 'bin');
        idx_full = find(~isnan(idx));
end




