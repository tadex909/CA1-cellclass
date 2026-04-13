

function slidstab = fct_slidstab_old(rate_tr_bin, shwin)

% shwin = 2; % size half-window in points

nb_tr = size(rate_tr_bin, 1);

nb_bin = size(rate_tr_bin, 2);

%% Create an Index Matrix to perform the mean
init = 1:(2*shwin+1);

mat = repmat(init, nb_bin - (2*shwin + 1), 1);
mat = bsxfun(@plus, mat, (1:(nb_bin - (2*shwin + 1)))');

edgeup = repmat(init, shwin + 1, 1); 
edgedn = repmat(mat(end, :), shwin, 1); 

mat = [edgeup ; mat ; edgedn];

%%
slidstab = NaN(1, nb_bin);


for k  = 1:nb_bin
    
    % Correlation / mean
%     meanvec = nanmean(rate_tr_bin(:, mat(k, :)), 2);
%     curvec = rate_tr_bin(:, k);
%     slidstab(k) = fct_det_r(meanvec, curvec);
    
    % Pairwises
    
%     mattmp = rate_tr_bin(:, mat(k, :));
%     
%         allpairs = nchoosek(1:(2*shwin + 1), 2);
%         nb_pairs = size(allpairs, 1);
%         vec_stab = zeros(1, nb_pairs);
%         
%         for kk = 1:nb_pairs
%             vec1 = mattmp(:, allpairs(kk, 1));
%             vec2 = mattmp(:, allpairs(kk, 2));
% 
%             vec_stab(kk) = fct_det_r(vec1, vec2);
%         end
%         slidstab(k) = nanmean(vec_stab);
        
        
        
           % Pairwises
%             mattmp = rate_tr_bin(:, mat(k, :));
%     
%         allpairs = nchoosek(1:nb_tr, 2);
%         nb_pairs = size(allpairs, 1);
%         vec_stab = zeros(1, nb_pairs);
%         
%         for kk = 1:nb_pairs
%             vec1 = mattmp(allpairs(kk, 1), :);
%             vec2 = mattmp(allpairs(kk, 2), :);
% 
%             vec_stab(kk) = fct_det_r(vec1, vec2);
%         end
%         slidstab(k) = nanmean(vec_stab);
        
        
        
        
          mattmp = rate_tr_bin(:, mat(k, :));
    iimp = nanmean(mattmp(1:2:end, :), 1);
     ipp = nanmean(mattmp(2:2:end, :), 1);
   
        slidstab(k) = fct_pearson(iimp, ipp);
end
















