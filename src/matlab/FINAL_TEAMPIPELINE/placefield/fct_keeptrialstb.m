function idx = fct_keeptrialstb(ratem, prm)

if nargin == 2
    prm.th_stb = 0.60;
    prm.th_nbtr = 0.40;
end


[~, stb_tr] = fct_stability_index(ratem, 'meancorr');
nb_tr = size(ratem, 1);
if (sum(stb_tr >= prm.th_stb) / nb_tr) >= prm.th_nbtr
    idx = find(stb_tr >= prm.th_stb);
else
    idx = NaN;  
end




