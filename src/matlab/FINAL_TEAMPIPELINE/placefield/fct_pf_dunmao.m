% -----------------------------
% Written by MARTI Geoffrey
% 09/17
% -----------------------------

function ispf = fct_pf_dunmao(rate_tr, prm)

%"Sparse orthogonal population representation of spatial context in the retrosplenial cortex Dun Mao1,2, Steffen Kandler1,3, Bruce L. McNaughton1,2 & Vincent Bonin1,3,4,5
% Nature com 2017
% First, criteria for place cell selection were adapted from previous literature30, 33. Briefly, each place cell had
%  to satisfy the following criteria:
% (1) Initial threshold was set at 30% of the difference between highest and lowest activity of the position tuning curve. Place fields must
%  be a continuous region with minimum 15 cm width and maximum 120 cm width.
% (2) The mean in-field activity must be at least three times larger than the mean out-
%  of-field activity.
% (3) The peaks of the position activity map across trials must be
%  within the potential field for at least one-third of all trials. Neurons that met these
%  criteria were selected as potential place cells.

% prm.th = 0.30;
% prm.nosth = 3;
% prm.nb_tr = 1/3;

rate = nanmean(rate_tr, 1);
nb_bin = size(rate_tr, 2);
nb_tr = size(rate_tr, 1);
[~, peak_tr] = nanmax(rate_tr, [], 2);

minr = nanmin(rate);
maxr = nanmax(rate);

ispf = rate > (minr + prm.th*(maxr - minr));

[iseq, ~, seql] = fct_find_seq(ispf);
critl = (seql < prm.min_len) | (seql > prm.max_len);
iseq(critl, :) = [];

nb_pf = size(iseq, 1);
ispf = zeros(1, nb_bin);
if fct_isactive(rate_tr, prm) && (nb_pf > 0)
    for k = 1:nb_pf
        ipf = iseq(k, 1):iseq(k, 2);
        ind_in = ipf;
        ind_out = setxor(1:nb_bin, ind_in);
        
        ratein = nanmean(rate(ind_in));
        rateout = nanmean(rate(ind_out));
        nos = ratein / rateout;
        
        if nos >= prm.nosth
            iint = (peak_tr >= iseq(k, 1)) & (peak_tr <= iseq(k, 2));
            if (sum(iint) / nb_tr) > prm.nb_tr
                ispf(ipf) = 1;
            end
        end
    end
    
end
