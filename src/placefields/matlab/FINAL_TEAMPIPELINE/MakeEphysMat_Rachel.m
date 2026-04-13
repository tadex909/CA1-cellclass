function []=MakeEphysMat_Rachel (AnlPath)


if nargin<1
    AnlPath='D:\Rachel\RAM12_2024-04-09_11-12-39';
end

% ____________________Where are the data?_____________________

% CircusPath=('X:\MATLAB\Data\Data_raw\VS15\VS15_2022-03-02_18-18-27\Continuous_Data\Continuous_Data-merged.GUI');%where circus data is .GUI !!!

slashes=strfind(AnlPath,'\');
Namesession=[AnlPath(slashes(end)+1:end)];


CircusPath=[AnlPath, '\', Namesession, '\', Namesession ,'.GUI'];
OutputFolder='D:\Rachel\Results' ;
cd (AnlPath)

% ispheno=0; %pour le moment  laisser comme ça c'est pour les nouvelles modifs 19/08/24 Vinca
%% basic data extraction
if exist(strcat(Namesession,'_Phenosys.mat')) ~=2 % if it exists in the folder as a file
    ispheno=0;
    disp('Cannot find Phenosys.mat')
    disp('Please run MakeTrajFile_Rachel first')
    
else
    ispheno=1;
    load (strcat(Namesession,'_Phenosys.mat'))
end

if exist(strcat(Namesession, '_TrajData.mat'))==2
    load(strcat(Namesession, '_TrajData.mat'));
    istraj=1;
else
    disp('Cannot find TrajData.mat')
    disp('Please run MakeTrajFile_Rachel first')
    istraj=0;
    
end
%TESTS PB CIRCUS REMETTRE APRES
if exist(strcat(Namesession,'_Ratemap.mat')) ==2
    load(strcat(Namesession,'_Ratemap.mat'));
    isplacecell =1;
else
    isplacecell=0;
end


% uncomment bellow for unit analysis so that you can run code from geoffrey and Romain
% unit analysis
     [Spikes,phyData,goodClusListIdx]=Load_phyResults(CircusPath);
        cd (CircusPath)
%         [datacell, headercell,outputt]=tsvread2('cluster_info.tsv');
  

     [Spikes,phyData,goodClusListIdx]=Load_phyResults(CircusPath);
    cd (AnlPath);
    config = parse_params_file([Namesession '.params']);
    %get info from the structure.oebin
    fname = 'structure.oebin';
    fid2 = fopen(fname);
    raw = fread(fid2,inf);
    str = char(raw');
    fclose(fid2);
    val = jsondecode(str);
    header=val.continuous(1);% plenty of info here
    
    cd (CircusPath);

     data=readtable('cluster_info.tsv','FileType','text','Delimiter', '\t');
     CID=data.cluster_id(goodClusListIdx);% remove base 0?
%     data=readtable('cluster_info.tsv','FileType','text','Delimiter', '\t');
 SF=config.sampling_rate;
allcel.id_cel = CID;
allcel.itime_spk = double(Spikes.SpikeTimes);
allcel.time_spk = double(Spikes.SpikeTimes)./SF;
allcel.time_spk2 = double(Spikes.SpikeTimes)./(SF/1000);
allcel.id_spk = Spikes.Units+1; %NEW +1 pour quitter 0base et faire1base pour matlab (22/08/24 Vinca)
allcel.Spikes = Spikes;

tcid=CID(1);% start with first cell
indexes=find(Spikes.Units==tcid);
SpikeTimes=double(Spikes.SpikeTimes(indexes)); %if I use the dat cele, in samples
        
        Channel= data.ch(data.cluster_id==tcid)+1; %truc de fou... faut trouver lequel est le canal que data donne dans channel map
        shank=data.sh(data.cluster_id==tcid)+1;
        [WVMat]=DatWV (shank, AnlPath,SpikeTimes, config, header );

 allcel.waveform=nan(size(WVMat,1),size(WVMat,2),size(WVMat,3),length(allcel.id_cel));
%     this guy should be: 150 by 16 by 50 by ncells
%    the 50 is an issue since I can get less spikes
    
 for cel=1:length(CID)
        
    disp (['Cell ', num2str(cel), ' id ', num2str(CID(cel),'%d')]) ;
    cellID=CID(cel);
        indexes=find(Spikes.Units==cellID);
        if ~isempty (indexes)
        %     there is an issue when the cell fires less than 50APs
        
        SpikeTimes=double(Spikes.SpikeTimes(indexes)); %if I use the dat cele, in samples
            if ~isempty (SpikeTimes)
                Channel= data.ch(data.cluster_id==cellID)+1; %truc de fou... faut trouver lequel est le canal que data donne dans channel map
                shank=data.sh(data.cluster_id==cellID)+1;
                
                [WVMat]=DatWV (shank, AnlPath,SpikeTimes, config,header );
            else
                disp('WTF!')
            end
        end
allcel.waveform(:,:,:,cel)= WVMat;  
% allcel.meanwaveform () = squeeze(mean(WVMat
end
allcel.meanwaveform=squeeze(mean(allcel.waveform,3));     

%% Save outputs:

 save([Namesession '_ePhy.mat'],'allcel')

% save([OutputFolder '\' Namesession '_TrajData.mat'],'Traj','newW')
% save([OutputFolder '\' Namesession '_Phenosys.mat'],'t','X','transition','Sess', 'Whl','starts','sessS','SessStarts','Speed','XSpeed','transitionThresh')

%


