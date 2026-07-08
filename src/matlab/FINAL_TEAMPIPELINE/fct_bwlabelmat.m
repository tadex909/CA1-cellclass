% Author(s): Bourboulou Romain, Marti Geoffrey   
% Epsztein Lab 2019

% Given a matrix of zeros and ones, this function will label each
% contiguous sequence of ones with index from 1 (first sequence) to N (N-th
% sequence). This operation is performed for each matrix row independently.


function mat_lab = fct_bwlabelmat(mat_in)

sz = size(mat_in);

% Label field per laps 
% pad zero to split field that are at the beginning and end.
mat_lab = [mat_in zeros(sz(1), 1)]';

% vectorise vector and label connected components. Then reshape in the
% original form
mat_lab = bwlabel(mat_lab(:));
mat_lab = reshape(mat_lab ,  sz(2)+1 , sz(1));
mat_lab = mat_lab(1:end-1,:)';

% remove minimum per row to obtain the number of field per row and not in
% total
idx = (mat_lab == 0);
mat_lab(idx) = NaN;

%%%ORIGINAL CODE%%%%
%{
crap = nanmin(mat_lab , [] , 2);
mat_lab = (mat_lab - nanmin(mat_lab , [] , 2)) + 1;
mat_lab(idx) = 0;
%}

%bug in my code (julie): 
%I replace by : 
crap = nanmin(mat_lab , [] , 2);
crap2 = nan(length(mat_lab(:,1)),length(mat_lab(1,:)));
for tt=1:length(crap)
 crap2(tt,:) = crap(tt);    
end
mat_lab = crap2 +1;
mat_lab(idx) = 0;





