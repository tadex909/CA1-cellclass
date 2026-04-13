function [si_corr, si_suro, si] = fct_si_suro(frmap, tsmap)


% isgood_tr= fdt_find_waycond(1, eprm.icond_tr, bhv.way_tr);
% 
% 
% frmap = cel(1).frmap_nr.rate_s_tr_bin(isgood_tr, :);
% tsmap = cel(1).frmap_nr.ts_s_tr_bin(isgood_tr, :);

%% Surrogate Distribution
nb_rep = 100;
frmap_suro = frmap;
tsmap_suro = tsmap;


idx = find(~isnan(frmap_suro) & ~isnan(tsmap_suro));
si_suro = NaN(1, nb_rep);

for p = 1:nb_rep
    ip = randperm(length(idx));
    
    idx_perm = idx(ip);
    
    frmap_suro(idx) = frmap(idx_perm);
    tsmap_suro(idx) = tsmap(idx_perm);
    
    
  
            si_suro(p) = fct_spatial_info(nanmean(frmap_suro), nanmean(tsmap_suro));
   
    
    %     figure
    %     subplot(221)
    %     imagesc(frmap)
    %
    %     subplot(222)
    %     imagesc(frmap_suro)
    %
    %     subplot(223)
    %     imagesc(tsmap)
    %
    %     subplot(224)
    %     imagesc(tsmap_suro)
end

%% Spatial Info
si = fct_spatial_info(nanmean(frmap), nanmean(tsmap));

%% Normalized Spatial Info
si_corr = (si - nanmean(si_suro)) / nanstd(si_suro);











