function [Spikes,phyData,goodClusListIdx]=Load_phyResults(phyOutputFolder)
curDir=cd; %keep current directory in memory
% phyOutputFolder=cd;
cd(phyOutputFolder);

%list files in phy GUI outputs
dirListing=dir(phyOutputFolder);

%keep npy files
npyListing=dirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'.npy'),...
{dirListing.name},'UniformOutput',false)));

%read each file into a data structure
phyData=struct;
for fileNum=1:size(npyListing,1)
%     disp(['reading', (npyListing(fileNum).name)])
phyData.(npyListing(fileNum).name(1:end-4))=readNPY(npyListing(fileNum).name);
end

%get cluster_id group from csv file
tsvFile=dirListing(~cellfun('isempty',cellfun(@(x) strfind(x,'cluster_group.tsv'),...
{dirListing.name},'UniformOutput',false)));
formatSpec = '%s%[^\n\r]';
startRow = 2;
fileID = fopen(tsvFile.name,'r');
phyData.cluster_ids = textscan(fileID, formatSpec, 'HeaderLines' ,startRow-1);
fclose(fileID);
phyData.cluster_ids =[phyData.cluster_ids{:}];

% get the cluster_info file
% tsvFile2=fopen('cluster_info.tsv','r');
% [DD]=textscan(tsvFile2, '%d %5.2f %d %3.2f %s %s %d %d','HeaderLines' ,1);


%go back to original directory
cd(curDir);
clusgroup= phyData.cluster_ids(:,2);
%find good ones
goodClusListIdx=cellfun(@(clusgroup) strcmp(clusgroup,'good'), phyData.cluster_ids(:,2));
goodClus=phyData.cluster_ids(goodClusListIdx,1);
goodClusIdx=ismember(phyData.spike_clusters,str2double(goodClus));
goodTemplates=unique(phyData.spike_clusters(goodClusIdx)+1); % FIX ME for some reason, need to add +1 to template number
% goodTemplates=unique(phyData.spike_templates(goodClusIdx)+1);
% electrodeIds=phyData.template_ind(goodTemplates)+1; %+1 to convert electrode number back to proper base index
electrodeIds=zeros(size(goodTemplates));
for id=1:length(goodTemplates)
electrodeIds(id)= str2double(cell2mat(phyData.cluster_ids(goodClusListIdx(id),1)))+1;
end
% electrodeIds=str2double(cell2mat(phyData.cluster_ids(goodTemplates,1)))+1; %+1 to convert electrode number back to proper base index
Spikes.Electrodes=nan(length(goodClusIdx),1);
for tempNum=1:length(goodTemplates)
Spikes.Electrodes(phyData.spike_templates==(goodTemplates(tempNum)))=electrodeIds(tempNum);
end
%store data about good clusters
Spikes.Units=phyData.spike_clusters(goodClusIdx);
Spikes.SpikeTimes=phyData.spike_times(goodClusIdx);
Spikes.Electrodes=uint16(Spikes.Electrodes(goodClusIdx));
Spikes.PC=phyData.pc_features(goodClusIdx,:,:); % it's bugging here, it's not the good size. j'ai remis vinca 08/08/24


% %plots
% figure;
% for clusNum=1:length(goodClus)
% goodClusIdx=phyData.spike_clusters==str2double(goodClus{clusNum});
% goodTemplates=unique(phyData.spike_templates(goodClusIdx));
% electrodeIds=phyData.template_ind(goodTemplates+1)+1;
% subplot(1,length(goodClus),clusNum); hold on;
% for tempNum=1:length(goodTemplates)
% plot(squeeze(phyData.templates(goodTemplates(tempNum)+1,:,:)))
% legend(['Electrode ' num2str(electrodeIds) ...
% ', template ' num2str(goodTemplates(tempNum))]);
% end
% end

 