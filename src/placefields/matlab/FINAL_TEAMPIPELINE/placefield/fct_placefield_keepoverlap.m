% Author(s): Bourboulou Romain, Marti Geoffrey   
% Epsztein Lab 2019

% This is a super cool function that will take a matix of place field per
% lap (mat_in) and a vector of principal median field (majField). It will return a matrix of
% same size than mat_in with lapfields labeled with the id of the major
% field they overlap the most.
% 
% INPUT:
%       - mat_in: 'logical' matrix of fields per lap (lap * space bins) 
%       - majField : logical or labeled vector of major fields
% 
% OUTPUT:
%       - mat_out: labeled matrix of fields per lap (lap * space bins) with
%       label corresponding with the major field they overlap the most
%


% mat_in = [0 1 0 0 0 0 0; 0 1 0 1 1 0 0 ; 0 0 0 0 0 0 1 ; 0 0 0 0 0 1 0 ; 0 0 0 0 0 0 0 ; 0 1 0 1 0 0 0 ; 0 1 0 1 0 0 1 ; 0 0 0 0 1 1 1 ; 0 0 0 1 1 1 1];
% majField = [0 2 2 2 2 0 1];
% mat_out = fct_placefield_keepoverlap(mat_in , majField)
% figure
% subplot(311)
% imagesc(mat_in)
% subplot(312)
% imagesc(mat_out)
% subplot(313)
% imagesc(majField)


function mat_out = fct_placefield_keepoverlap(mat_in , majField)



%% Initialise values
sz = size(mat_in);

% check input
assert(sz(2) == size(majField, 2), 'size of the input do not match')

% add the possiblity to affect the sub field per lap to the major field (id
% = 1) default behavior will affect to the left side field
if any(majField > 1)
    majField_l = majField;
else
    majField_l = bwlabel(majField);
end
nFieldId = max(majField_l);

if nFieldId == 0
    mat_out = zeros(sz);
    return
end

%% Label field per laps 
mat = fct_bwlabelmat(mat_in);

% find max number of fields per lap 
nmaxfield = nanmax(mat(:));
if nmaxfield == 0
    mat_out = zeros(sz);
    return
end

% Make bigMat of dim [lap * spacebin * number of field per lap]
repMat_in = repmat(mat, 1, 1, nmaxfield);
fieldIdMat = repmat(reshape(1:nmaxfield , 1 , 1 , []) , [sz(1) sz(2)  1]);
bigMat = repMat_in == fieldIdMat;

%%
sumForField = zeros(sz(1), nmaxfield, nFieldId);

for f = 1:nFieldId
    tmpField = majField_l == f; 
    tmpField3D = repmat(tmpField, sz(1), 1, nmaxfield);
    tmpField3D = bigMat & tmpField3D; 
    overlapF = reshape(nansum(tmpField3D, 2), sz(1) , nmaxfield);

    % Put Nan when there is no overlap
    % overlapF: the number of overlapped bins of each pf
    % per lap with the mean pf
    overlapF(overlapF == 0 ) = NaN; 

    sumForField(: , : , f) = overlapF;
end
% find the max overlap 
[val , mIDx] = nanmax(sumForField , [] , 3);
mIDx(isnan(val)) = NaN ;
mIDxRep = repmat(reshape( mIDx , sz(1) , [] ,  nmaxfield ) , 1 ,sz(2), 1 );

% make output
mat_out = nansum(mIDxRep.*bigMat, 3);

