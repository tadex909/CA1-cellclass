

function stab_bin = fct_stabin(rate_tr_bin, nb_bin, mtd)
% nb_bin = size(rate_tr_bin, 2);
% nb_tr = size(rate_tr_bin, 1);
%
%
% stab_bin = NaN(1, nb_bin);
%
%
% if isempty(rate_tr_bin)
% end


% for k = 1:nb_bin
%     mean_bin = nanmean(rate_tr_bin(:, k))*ones(nb_tr, 1);
%     stab_bin(k) = fct_det_r(rate_tr_bin(:, k), mean_bin);
% end


if ~isempty(rate_tr_bin)
    switch mtd
        case 'meannorm' % moyenne normalisée
            nrate_tr_bin = fct_matnorm(rate_tr_bin, 2);
            stab_bin = nanmean(nrate_tr_bin, 1);
        case 'meannormcorr' % moyenne normalisée avec 0 remplacé par le min de chaque colonne pour corriger l'indice en cas de non activité
            nrate_tr_bin = fct_matnorm(rate_tr_bin, 2);
            for k = 1:size(rate_tr_bin, 2)
                minval = min(nrate_tr_bin(nrate_tr_bin(:, k) > 0, k));
                if ~isempty(minval)
                nrate_tr_bin(nrate_tr_bin(:, k) == 0, k) = minval;
                end
            end     
            stab_bin = nanmean(nrate_tr_bin, 1);
            
        case 'stdnorm' % std normalisée
            nrate_tr_bin = fct_matnorm(rate_tr_bin, 2);
            stab_bin = 1 - nanstd(nrate_tr_bin, 1);
            
        case 'varcoef' % inverse du coefficient de variation
            stab_bin = nanmean(rate_tr_bin, 1) ./ nanstd(rate_tr_bin, 1) ;
            
        case 'idisp' % inverse de l'indice de dispersion
            stab_bin = nanmean(rate_tr_bin, 1) ./ (nanstd(rate_tr_bin, 1).^2);
            
        case 'idispnorm' % inverse de l'indice de dispersion normalisé
            nrate_tr_bin = fct_matnorm(rate_tr_bin, 2);
            stab_bin = nanmean(nrate_tr_bin, 1) ./ (nanstd(nrate_tr_bin, 1).^2);
            
            
    end
else
    stab_bin = NaN(1, nb_bin);
end