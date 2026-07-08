% Author(s): Marti Geoffrey
% Epsztein Lab 2019


% This function will find for each element of the vector 'vec' its closest
% value/bin in the vector 'bins' and return the index of that value/bin.
% The algorithm will look for closest values by default. 
% You can put  'bin'  as third argument to get the closest bin 
% (e.g., if bins = [0 2 4 6], then [0 2[, [2 4[ and [4 6[  bins will be considered).
% bins can also be a matrix of intervals (col1 = start, col2 = stop, row = interval number) 
% Finally, you can add 'outedgebincount' as 4th argument to take into account values 
% outside the edges (to be used with 'bin' option)



function ivec = fct_discretize(vec, bins, varargin)


[vec, bins, prm] = dealin(vec, bins, varargin);


% vec = [-10 -15 -14 0.1 0.2 0.3 2.1 2.2 2.3 2.56 3.1 4 10 15 20 30 12];
% bins = [0 1 2 3 4];

% vec = [-5 0 0.1 0.8 1.5 1.8 3.1 3.5 4.5 10 15 16 18 25 59 70 78];
% bins = [0 1 ; 3 9 ; 20 60];


switch prm.option
    case 'closest' % look for the closest values from the vector vec in the vector bins
        switch prm.mtd 
            case 1 % bins values have to be equally spaced
                ts = mean(diff(bins));
                bins = bins - 0.5*ts;
                bins = [bins (bins(end) + ts)];
                
                bins(1) = -inf;
                bins(end) = +inf;
                
                [~, ivec] = histc(vec, bins);
            case 2 % find the closest values by error minimization but needs much memory  
                tmp = abs(bsxfun(@minus, vec, bins'));
                [~, ivec] = min(tmp);   
        end
             
    case 'bin'
        if prm.outedgebincount
            bins(1) = -inf;
            bins(end) = +inf;
        end
        
        [~, ivec] = histc(vec, bins);
        ivec(ivec == 0) = NaN;   
        
        % very special case when there is only one bin which does not cover
        % all the vector (e.g., fct_discretize(1:10, [0 9], 'bin'))
        if numel(bins) == 2
            ivec(ivec == 2) = 1;
        end
        
        if prm.isbinmatrix
            ivec(mod(ivec, 2) == 0) = NaN;
            ivec(mod(ivec, 2) == 1) = floor(ivec(mod(ivec, 2) == 1) / 2) + 1;
        end
end


end


function [vec, bins, prm] = dealin(vec, bins, X)

prm.isbinmatrix = false;
prm.option = 'closest';
prm.mtd = 2;
prm.outedgebincount = false;

if iscolumn(vec)
    vec = vec';
end
 
if ~isvector(bins)
    prm.isbinmatrix = true;
    if size(bins, 2) ~= 2
        warning('Bins matrix bad conditionned, it should be col = start/stop index, row = #interval')
    end
    bins = bins';
    bins = bins(:)';
end


if any(strcmp(X, 'bin'))
    prm.option = 'bin'; 
end

if strcmp(prm.option, 'closest')
    if (abs(std(diff(bins))) <= 1000*eps)
        prm.mtd = 1;    
    end
end

if any(strcmp(X, 'outedgebincount'))
    if strcmp(prm.option, 'bin')
        prm.outedgebincount = true;
    else
        prm.outedgebincount = false;
        warning('''outedgebincount'' option only allowed when ''bin'' option enabled')
    end
end
end


