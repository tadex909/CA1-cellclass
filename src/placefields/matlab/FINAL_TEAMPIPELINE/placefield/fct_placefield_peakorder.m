% Author(s): Marti Geoffrey   
% Epsztein Lab 2019

% Given a firing rate ('fr_s_x') and a boolean vector of place fields
% ('ispf_x') per spatial bin, this function will assign a label to each field, from the one
% with the highest peak firing rate (label 1, principal field) to the one
% with the lowest peak (label N, N being the total number of fields).

function ispflab_x = fct_placefield_peakorder(fr_s_x, ispf_x)


ispflab_x = bwlabel(ispf_x);
nb_pf = max(ispflab_x);

if nb_pf == 0
    ispflab_x = ispf_x;
    return
end
    
[~, idx] = sort(accumarray((ispflab_x + 1)', fr_s_x, [], @nanmax), 'descend');
idx(idx == 1) = [];
ispflab_x_rep = repmat(ispflab_x, nb_pf, 1);
ispflab_x = sum(bsxfun(@times, ispflab_x_rep == (1:nb_pf)',  (idx - 1)), 1);





