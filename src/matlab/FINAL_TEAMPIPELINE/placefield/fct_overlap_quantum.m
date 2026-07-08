% -----------------------------
% Written by MARTI Geoffrey
% 03/16
% -----------------------------

% This function is intended to quantify the overlap between two vectors ispf_way1 and ispf_way2 
% including one placefield (only one)

function overlap_quantum = fct_overlap_quantum(ispf_way1, ispf_way2)


if isequal(ispf_way1, ispf_way2)
    overlap_quantum = length(find(ispf_way1)); 
else
    pf1_bound = [find(ispf_way1 == 1, 1, 'first') find(ispf_way1 == 1, 1, 'last')];
    pf2_bound = [find(ispf_way2 == 1, 1, 'first') find(ispf_way2 == 1, 1, 'last')];
    if pf1_bound(2) == pf2_bound(2) % e.g., ispf_way1 = [0 0 0 1 1 1] and ispf_way2 = [0 0 0 0 1 1]
        [~, imax_track] = min([pf1_bound(1) pf2_bound(1)]);
        imax_track = fct_mod(imax_track + 1, 2); 
    else
        [~, imax_track] = max([pf1_bound(2) pf2_bound(2)]);
    end
    if imax_track == 1
        overlap_quantum = pf2_bound(2) - pf1_bound(1) + 1;
    else
        overlap_quantum = pf1_bound(2) - pf2_bound(1) + 1;
    end
end


