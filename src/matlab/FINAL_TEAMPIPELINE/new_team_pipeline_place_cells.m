    % ____________________Where are the data?_____________________
AnlPath='X:\MATLAB\Data\object_def _good_2025-05-07_16-19-42\Record Node 103\experiment1\recording1\continuous\Rhythm_FPGA-100.0';
CircusPath=('X:\MATLAB\Data\Data_raw\VS109\VS109_2024-12-20_17-20-40\Record_Node_104\experiment1\recording1\continuous\Rhythm_FPGA-101.0\continuous\continuous-merged.GUI');%where circus data is .GUI !!!
Namesession= 'VS109_2024-12-20_17-20-40';
whatprob=77; %can only be numbers matching cases in probechannels function 5,64 or 77 for now

OutputFolder='N:\Vinca\MATLAB\Data\All_data\process'; % where you want your results

    % ____________________What are the data?_____________________
% please CondNam can be 'PO' 'PNO' 'PO_NOTREE' ' POnM' 'POnM9' 'POMb' 'POM'
CondNam={'PO' 'PO' 'PNO' 'PO' 'PO'};
Circu=2*pi*15; %ball circumference  %%REINCLURE DAnS CODE
%  Circu=43; %old wheel circonf
cd (AnlPath);
nospike=0; %put 1 if you want to do bhv only CHANGE AS WANTED
isoe=0; %create a variable to classify as open ephys format (=1) or not (=0, means .dat binary format) DO NOT MODIFY
A=dir;
oebpath=pwd;

%% basic data extraction
cd(strcat(AnlPath,'\process'));
if exist('Phenosys.mat') ~=2 % if it exists in the folder as a file
    ispheno=0;    
else
    ispheno=1;
    load ('Phenosys.mat')
end

if exist('TrajData.mat') ==2
    load('TrajData.mat');
    istraj=1;
else
   plotting=1;
   istraj=0;
end
%TESTS PB CIRCUS REMETTRE APRES
if exist(strcat(Namesession,'_Ratemap.mat')) ==2
    load(strcat(Namesession,'_Ratemap.mat'));
    isplacecell =1;
else
        isplacecell=0;
end
cd (AnlPath);
%% Data extraction :  Continuous Open Ephys file format 
        %____________________preparation_____________________
    if startsWith(A(4).name,'10');  %NEW !! haha trop contente de cette idée : les 3 formats commence par 100 donc valide tous les 100_trucmuch mais comme j'ai des sessions avec 101_ADC et pas 100 au moins ca marchera avec elles aussi ! 
        isoe=1; 
        if strcmp( A(4).name(4:5),'_1') %NEW j'ai changé pour que ça marche avec tous les processor id voir premier if de cette section. 
            prefix =  A(4).name(1:4);
            fff(1)=65; %XposVR
            fff(2)=66;% transitions in Y -back and fourth trip 
            fff(3)=67;% Xpos real
            fff(4)=68;%Sess??
        elseif  strcmp(A(4).name(4:7),'_CH1') %NEW j'ai changé pour que ça marche avec tous les processor id voir premier if de cette section. 
            prefix =  A(4).name(1:6);
            fff(1)=65;
            fff(2)=66;
            fff(3)=67;% Xpos real
            fff(4)=68;
        elseif strcmp(A(4).name(4:8),'_ADC1')
            prefix = A(4).name(1:7);
            fff(1)=1;
            fff(2)=2;
            fff(3)=3;% Xpos real
            fff(4)=4;

        end
    end


    %____________________file opening_____________________
    if isoe 
        tic
        %_____________________way/back data_____________________
        Nam=[prefix,num2str(fff(2)),'.continuous'];
        [Y, ty, info] = fct_read_continuous_file(fullfile(AnlPath, Nam)); 
        sf1=info.sampleRate; 
        sf2 =sf1/1000; 
        Y_ds=downsample(Y,sf2);
        t_y_ds=ty(1:sf2:end); 

        %_____________________position data_____________________

        Nam=[prefix,num2str(fff(1)),'.continuous'];
        [X, tx, info] = fct_read_continuous_file(fullfile(AnlPath, Nam)); %NEW
        FirstTime=tx(1);
        X_ds=downsample(X,sf2); %was 30, now auto
        tx_ds=tx(1:sf2:end);

        %_____________________wheel data_____________________

        Nam=[prefix,num2str(fff(3)),'.continuous'];
        [Whl_u, twl, ~] = fct_read_continuous_file(fullfile(AnlPath, Nam)); %NEW
        %
        Whl_ds=downsample(Whl_u,25);
        t_wl_ds=twl(1:25:end);
        
        %_____________________wheel data_____________________

        Nam=[prefix,num2str(fff(4)),'.continuous'];
        [Y_r, tses, info] = fct_read_continuous_file(fullfile(AnlPath, Nam)); %NEW
        Y_r_ds=downsample(Y_r,25);
        

        sf=sf1/sf2; %so 1000Hz...
        si=(1/sf1)*1E6; % in microseconds
        
        t_ds=tx_ds; %arbitraire pourrait être ty ou twl
 toc
    end
   
    clear Nam 

    
%% Data extraction : continuous.dat Binary file format  

%_____________________initialize thanks to oebin_____________________

    if  ~isoe

            oebpath=pwd; 


            fname = 'structure.oebin';
            fid = fopen(fname);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            val = jsondecode(str);
            header=val.continuous(1);
            sr = val.continuous.sample_rate;
            idADC=  strfind(arrayfun(@(x) startsWith(val.continuous.channels(x).channel_name,'ADC') , 1:val.continuous.num_channels),1);  
    end
    

  if  ~isoe
%_____________________localise .dat file and read bhv_____________________
    tic
    if isfile('continuous.dat')
        [data, t_ds]=Read_OEP_Binary ('continuous.dat',[idADC(1:4)], 0,-1, sr/1000);
    elseif ~isfile('continuous.dat')  
        path_dat = strcat(AnlPath, filesep, 'continuous\Rhythm_FPGA-', num2str(val.continuous.source_processor_id), '.', num2str(val.continuous.source_processor_sub_idx), filesep); 
        [data, t_ds]=Read_OEP_Binary_altered(strcat(path_dat,'continuous.dat'),idADC(1:4), 0,-1, sr/1000);
            if ~isfile(strcat(path_dat, 'continuous.dat'))
                error('file continuous.dat not find, please check location')
            end
    end

    if length(t_ds)~=(size(data,2))
        disp('OUPS TIMESTAMPS ARE FALSE !!! REPAIRING')
        t_ds=linspace(t_ds(1), t_ds(1)+ size(data,2)/1000,size(data,2));
    end

    toc

%_____________________make bhv variables_____________________

    Y_ds=data(2,:);
    X_ds=data(1,:);
    Whl_ds=data(3,:);
    Y_r_ds=data(4,:);

   
    si=(t_ds(2)-t_ds(1)); % in microseconds
    sf=1/si;
    si=round(si*1E6);
    FirstTime=t_ds(1); 
 end

%%  Vérification des variables visuellement 
figure
subplot(3,1,1)
plot(t_ds,X_ds)
hold on

subplot(3,1,2)
plot(t_ds,Y_ds)
hold on
subplot(3,1,3)
plot(t_ds,  Y_r_ds)


%% Compute behavior info 
%%%%Start at 0

    if min(X_ds)<0
        X_ds_n=X_ds+abs(min(X_ds));  %X downsamplé normalisé
    elseif min(X_ds)>=0
        X_ds_n=X_ds-(min(X_ds));
    end

    %
    X_ds_n=(X_ds_n./max(X_ds_n))*145; % to fit what the animal runs for real, like if gain=1


    %%Start at 0
    if min(Whl_ds)<0
        Whl_ds_n=Whl_ds+(abs(min(Whl_ds)));
    elseif min(Whl_ds)>=0
        Whl_ds_n=Whl_ds-((min(Whl_ds)));
    end
    %

    if max(Whl_ds_n)>4
        Turn=2*pi*15;%in cm
        VoltTurn = 1725.75;%in whl units estimated by 10 real turns of wheel by julie on April 2024
        Whl_u=(Whl_ds_n.*Turn)./VoltTurn; %Wheel downsamplé normalisé unwrappé
        transitionThresh=2600;%%

    %     highYlevel=2900; % for info, to get rid of middle y values at the start of a session
    else
        Turn=2*pi*15;%in cm
        VoltTurn=1.3435;%in whl units estimated by 50 real turns of wheel by julie April 2024
        Whl_u=(Whl_ds_n.*Turn)./VoltTurn;
        %transitionThresh=2;%OLD pose pb avec détection resets (sur-détecté avec tT=2)
        transitionThresh=2;
    %     highYlevel=3;% voir avec un enregistrement sans rec?
    end

    % detect when wheel data jumps
    % so that I patch it back together after

    thresh=5;% can't decently run 5cm in 1ms

    jumpsup=find(diff(Whl_u)>thresh);
    jumpsdown=find(diff(Whl_u)<-thresh);
    Whl2=Whl_u;
    if ~isempty(jumpsup)
        for jj=1:length(jumpsup)
            dist2remove=abs(Whl2(jumpsup(jj)+1)-Whl2(jumpsup(jj)));
            Whl2(jumpsup(jj)+1:end)=Whl2(jumpsup(jj)+1:end)-dist2remove;
            Whl2=Whl2+abs(min(Whl2));
        end

        for jj=1:length(jumpsdown)
            dist2add= abs(Whl2(jumpsdown(jj)+1)-Whl2(jumpsdown(jj)));
            Whl2(jumpsdown(jj)+1:end)=Whl2(jumpsdown(jj)+1:end)+dist2add;
        end

        Whl_u=Whl2.*-1;%it was going dozn I need it to go up
    else
        disp('BAAD ISSUE WITH Wheel acquisition!!!')

    end

    Whl_u=Whl_u+abs(min(Whl_u));% put it to zero

    Whl_u=Whl_u(:);
    t_ds=t_ds(:);
    X_ds_n=X_ds_n(:);
    Y_r_ds=Y_r_ds(:);

%% get separation between sessions S1-S2-S3

hf=figure;
ax1=subplot(2,1,1);
plot(Y_ds,'k')
hold on
set(hf,'Units','normalized');
set(hf, 'Position',[0.05 0.05 0.8 0.5])

NSess=str2double(cell2mat(inputdlg('How Many Sessions?','Sessions',1,{'1'})));
prompts=cell(1,NSess);

prompts(1:end)=compose('Name for Sess %d', 1:NSess);
Defaults=cell(1,NSess);
Defaults(1:end)={'PO'};
CondNam=inputdlg(prompts,'Names',1,Defaults);

close(hf)

hf=figure;
set(hf,'Units','normalized');
set(hf, 'Position',[0.05 0.05 0.8 0.5])

times=zeros(NSess,2);
sessS=zeros(NSess,1);

t2=t_ds-t_ds(1);

for c=1:NSess%tryin to plot a piece at a time, for precision
    idx=find(t2>times(c,1)& t2<times(c,1)+(max(t2)));
    plot(t2(idx),Y_ds(idx))
%     xlim ([t2(idx(1)),t2(idx(end))])

    hold on
    title(['Click on Session',num2str(c,'%d'),' start'])
    [times(c,1),y,button]=ginput(1);
    plot(times(c,1),y,'or')
     
    title(['Click on Session',num2str(c,'%d'),' stop'])
    [times(c,2),y,button]=ginput(1);
    if c<NSess
    times(c+1,1)=times(c,2);
    end
    plot(times(c,2),y,'or')
    pause(1)
    hold off
    sessS(c)=find(t2>=times(c,1),1,'first');
end

SessStarts=times(:,1);

    %I use this to remove first session transition if it is too close from  0
    if ~isempty(sessS)
        if sessS(1)<(sf*30) %if it's less than 30 seconds from start... I could use more
            SessStarts(1)=[]; 
            SessStarts=[t_ds(sessS(1));SessStarts; t_ds(end)];
             sessS(1)=[];
    %         
        else
            SessStarts=[t_ds(1);SessStarts; t_ds(end)];
        end
    else
        %one condition
        SessStarts=[t_ds(1); t_ds(end)];
        disp('Session detection ISSUE Check your cables!!!!!!!')
    end

    trialstarts=find(abs(diff(X_ds_n))>50); %OLD (transition2))>0.03
    remove=find(diff(trialstarts)==1);%
    trialstarts(remove+1)=[];

    %% extract speed from wheel distance

    Whl2=smooth(Whl_u,300);
    Whl2=resample(Whl2,1,300);
    RealSpeed=(diff(Whl2).*1000)./300; % estimate speed in 300ms windows but in cm/sec
    RealSpeed(end-10:end)=0;
    RealSpeed(find(RealSpeed<=0))=0;
    Speed=interp1(t_ds(1:300:end-300),RealSpeed,t_ds); %put back at 1kHz...


    %% extract speed from VR distance
    jumpsup=find(abs(diff(X_ds_n))>50);
    if ~isempty(jumpsup)
        X2=X_ds_n;
        for jj=1:length(jumpsup)
            dist2add=abs(X_ds_n(jumpsup(jj)+2)-X_ds_n(jumpsup(jj)));
            X2(jumpsup(jj)+1:end)=X2(jumpsup(jj)+1:end)+dist2add;% here sometimes you get 2 samples that go through transition... annoying
            X2=X2-abs(min(X2));
        end
    else
        disp('BAAD ISSUE WITH VR DETECTION')

    end
     X3=X2;
    jumpsup=find(abs(diff(X2))>50);
    if ~isempty(jumpsup)
        for jj=1:length(jumpsup)
            dist2remove=abs(X2(jumpsup(jj)+2)-X2(jumpsup(jj)));
            X3(jumpsup(jj)+1:end)=X3(jumpsup(jj)+1:end)-dist2remove;
        end

    end



     X2=X3(1:300:end);
    X2Speed=(diff(X2).*1000)./300; % estimate speed in 300ms windows but in cm/sec
    X2Speed(find(X2Speed<=0))=0;
    XSpeed=interp1(t_ds(1:300:end-300),X2Speed,t_ds); %put back at 1kHz...
    starts=trialstarts;

    
     save( [AnlPath '\process\' 'Phenosys.mat'],'t_ds','X_ds','X_ds_n','Y_ds','Y_r_ds', 'Whl_u','starts','sessS','SessStarts','Speed','XSpeed','transitionThresh')
    save([  OutputFolder '\' Namesession '_Phenosys.mat'],'t_ds','X_ds','X_ds_n','Y_ds','Y_r_ds', 'Whl_u','starts','sessS','SessStarts','Speed','XSpeed','transitionThresh')

 
%% preparing a structure for behavior outputs
  sf=1000; % should be for rachel


cd (AnlPath)
if exist( 'TrajData.mat','file')
    
    answer = questdlg('Ovewrite TrajData file?','TrajData analysis',...
        'Yes','No','No');
    switch answer
        case 'Yes'
            
            [Traj, newW]= runTraj_vinca (t_ds, X_ds_n, Y_ds, Whl_u, starts, sessS, Speed, XSpeed,transitionThresh,AnlPath, OutputFolder, Namesession,CondNam)  ;
            
        case 'No'
            cd (AnlPath)
            disp('Using previously saved Phenosys.mat file')
            load ('TrajData.mat')
    end
else
    
    
    [Traj, newW]= runTraj_vinca (t_ds, X_ds_n, Y_ds, Whl_u, starts, sessS, Speed, XSpeed,transitionThresh,AnlPath, OutputFolder, Namesession,CondNam);
    %newW is reset at zero for each trajectory
end
% 
% if plotting
%     edgesX=[0:2.5:max(X)];
%     edgesW=[0:5:max(newW)];
%     plotmybehav (Traj,edgesX,edgesW, OutputFolder,Namesession,CondNam);
%     cd (AnlPath)
% end


%%

%Number of condition
nbsess=length(unique([Traj.Cond])); %max should be 5 based on the CondNam but things can happen and make the recording shorter

%get information for cells and clusters
[Spikes,phyData,goodClusListIdx,datacell,headercell,clusters_info_new]=getphyInfo(CircusPath); %J'ai modifié pour bien récupérer les infos du tsv car ça buguait chez moi.

if exist(strcat(AnlPath,'\',Namesession,'_ePhy.mat'))~=2
cluster_infos_tsv=[];
for i = 1:9:length(clusters_info_new);
    cluster_infos_tsv=[cluster_infos_tsv;[clusters_info_new(i:i+8)].'];
end


shkcol=find(strcmp(headercell,'sh'));%looking in the header what column corresponds to shank
idcol=find(strcmp(headercell,'cluster_id'));% and clusterID
datacell(1,:)=[]; %%%datacell is the table with all the units caracteristics
clusters=phyData.cluster_ids(:,idcol);

%extracting the name of every single cell
CID=zeros(size(clusters));
for cl=1:length(clusters)
    CID(cl)=str2double(cell2mat(clusters(cl,1)))+1; %NEW +1 pour quitter 0base et faire1base pour matlab (22/08/24 Vinca)
end


CID=CID(goodClusListIdx); %keeping only the good ones (as saved in Phy -'good')
%%
allcel.id_cel = CID;
allcel.itime_spk = Spikes.SpikeTimes;
allcel.time_spk = double(Spikes.SpikeTimes)/25000;
allcel.time_spk2 = double(Spikes.SpikeTimes)/25;
allcel.id_spk = Spikes.Units+1; %NEW +1 pour quitter 0base et faire1base pour matlab (22/08/24 Vinca)
allcel.Spikes = Spikes;

%OLD : 
%waveform = phyData.templates(CID,:,:); 
%waveform = phyData.templates(goodClusListIdx,:,:); 

%NEW : 
bestchannels=str2double([cluster_infos_tsv(2:end,3)])+1; %NEW +1 
bestchannelsunits=str2double([cluster_infos_tsv(2:end,1)])+1; %NEW +1 
bestchannelsunitsgood=[bestchannelsunits(strcmp(cluster_infos_tsv(2:end,6),'good'))]; %NEW +1 
bestchannelsgood=[bestchannels(strcmp(cluster_infos_tsv(2:end,6),'good'))]; 
cd(AnlPath)
shanksize=probechannels(whatprob,1);
allcel.waveform=nan(23,length(shanksize),50,length(allcel.id_cel));

for i=1:length(unique(Spikes.Units))
    disp(strcat('cell n°',string(i),'/',string(length(allcel.id_cel))))
    channeldata=[];
    idx =bestchannelsunitsgood(i);
    channel=bestchannelsgood(i);% logical indexing 
    mask=(double(Spikes.Units+1)==idx); %NEW -1 pour des questions de base zéro
    selectedtimespk= allcel.itime_spk(mask,1);
    channelshank=probechannels(whatprob,channel);

 %openephys format 
    for k=1:length(channelshank);    
        if isoe
             Nam=[prefix(1:end-3),'CH',num2str(channelshank(k)),'.continuous']; % OLD
%             Nam=[prefix,num2str(channelshank(k)),'.continuous']; 
%channell=find(channelshank(k)==chorder);
  %Nam=[prefix, num2str(channell),'.continuous'];
            %Nam=[prefix(1:end-3),'CH',num2str(channelshank(k)),'.continuous']; 
            [dischanneldata, ~, ~] = fct_read_continuous_file(fullfile(AnlPath, Nam)); 
            channeldata=[channeldata;dischanneldata];   
        elseif isfile('continuous.dat')
%             [dischanneldata, ~]=Read_OEP_Binary ('continuous.dat',channelshank, 0,-1, sr); %19/08/24 il faudrait remonter time window sans tout casser pour extraire que ces timestamps je l'ai pas fait désolée VINCA
%              channeldata=[channeldata;dischanneldata];   
            [channeldata, ~]=Read_OEP_Binary ('continuous.dat',channelshank, 0,-1, 1); %19/08/24 il faudrait remonter time window sans tout casser pour extraire que ces timestamps je l'ai pas fait désolée VINCA
        elseif ~isfile('continuous.dat')  %equivalent à si continuous.dat est ailleurs (pas sorti de rythmFGPA) ; ne marchera surement pas car val n'exisste pas si on a déja fait traj etc. j'ai une soluce en tête si besoin :) 
            path_dat = strcat(AnlPath, filesep, 'continuous\Rhythm_FPGA-', num2str(val.continuous.source_processor_id), '.', num2str(val.continuous.source_processor_sub_idx), filesep); %va chercher la localisation du .dat dans les sous-dossiers
            [dischanneldata, ~]=Read_OEP_Binary_altered(strcat(path_dat,'continuous.dat'),channelshank, 0,-1, sr); %la version "altérée consiste en la modification d'une ligne (D.timestamps) où j'ai ajouté le path_dat (comme au dessus) avant continuous.dat dans les inputs de ReadNPY pour qu'il trouve le fichier quin'est pas au même endroit que oebin ! c'est tout!
             channeldata=[channeldata;dischanneldata];   
            end
    end

    selectedamptimespk=selectedtimespk(1:50);
    for j=1:50 %j = nb de spike sélectionnés pour le sampling
        timewindow=[0:22]; %fenetre pour le spike 
        timewindow=timewindow+double(selectedamptimespk(j));
           if timewindow(1) <= 0
               timewindow=timewindow+abs(timewindow(1))+1;
           end
           subplot(25,25,j)
       allcel.waveform(:,:,j,i)=channeldata(:,timewindow).';
       plot(channeldata(:,timewindow))
    end
end
allcel.meanwaveform=squeeze(mean(allcel.waveform,3));
else 
    load([OutputFolder '\' Namesession '_ePhy.mat'])
end


save([AnlPath  '\process\ePhy.mat'],'allcel')
save([OutputFolder '\' Namesession '_ePhy.mat'],'allcel')
%save([OutputFolder '\' Namesession '_TrajData.mat'],'Traj','newW')
% save([OutputFolder '\' Namesession '_Phenosys.mat'],'t','X','transition','Sess', 'Whl','starts','sessS','SessStarts','Speed','XSpeed','transitionThresh')

%%
if ~isplacecell 
    Rmap_G
end