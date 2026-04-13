% Author(s): Marti Geoffrey   
% Epsztein Lab 2019

function [allrmap, allpf] = fct_mapandfields_reshape(rmap, pf)


fnames = fieldnames(rmap);
idx = find(~ismember(fnames, {'prm', 'feat'}))';
for k = idx
    allrmap.([fnames{k} 'u']) = cat(3, rmap.(fnames{k}));
end

fnames = fieldnames(rmap(1).feat);
idx = find(~ismember(fnames, {'localstab'}))';
for k = idx
    allrmap.feat.([fnames{k} '_uc']) = cell2mat(arrayfun(@(X) [X.feat.(fnames{k})], rmap(:), 'uni', 0));
end


fnames = fieldnames(pf);
idx = find(~ismember(fnames, {'prm', 'feat'}))';
for k = idx
    allpf.([fnames{k} 'u']) = cat(3, pf.(fnames{k}));
end

fnames = fieldnames(pf(1).feat);
idx = find(~ismember(fnames, {'lap_size_t', 'lap_ifrmax_t'}))';
nb_cond = rmap(1).prm.nb_cond;
for k = idx
    tmp1 = arrayfun(@(X) [X.feat.(fnames{k})], pf(:), 'uni', 0);
    tmp2 = cellfun(@(X) reshape(X, 5, nb_cond), tmp1, 'uni', 0);
    allpf.feat.([fnames{k} '_ucf']) = permute(cat(3, tmp2{:}), [3 2 1]);
end
