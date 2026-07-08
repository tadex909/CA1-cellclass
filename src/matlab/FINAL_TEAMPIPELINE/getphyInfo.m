function [Spikes,phyData,goodClusListIdx,datacell,headercell,outputt]=getphyInfo(CircusPath)
% gets all info from phy

% cd (CircusPath)
[Spikes,phyData,goodClusListIdx]=Load_phyResults(CircusPath);
cd (CircusPath)
[datacell, headercell,outputt]=tsvread2('cluster_info.tsv');
end

