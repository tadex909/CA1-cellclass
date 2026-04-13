%% Initialisation : Ajustez à vos data

    % ____________________Where are the data?_____________________
 %AnlPath='N:\Vinca\MATLAB\Data\Data_raw\VS57\VS57_2022-12-18_18-51-04\Record Node 104\experiment1\recording1'; %binary test pour PP section mettre doss avec continuous.dat pour vinca mettre le doss de struct.oebin
%AnlPath='N:\Vinca\MATLAB\Data\Data_raw\VS57\VS57_2022-12-19_16-30-47'; %EXEMPLE CONTINUOUS O.E.
%AnlPath='X:\MATLAB\Data\Data_raw\VS64\VS64_2023-02-14_18-04-10';
%AnlPath='X:\MATLAB\Data\Data_raw\VS88\VS88_2023-10-30_17-02-23\Record Node 104\experiment1\recording1'
AnlPath='N:\Vinca\MATLAB\Data\Data_raw\VS15\VS15_2022-03-01_16-02-15'
CircusPath=('N:\Vinca\MATLAB\Data\Data_raw\VS15\VS15_2022-03-01_16-02-15\Continuous_Data\Continuous_Data.GUI');%where circus data is .GUI !!!

OutputFolder=AnlPath; % where you want your results
Namesession= 'VS15_2022-03-01_16-02-15';


    % ____________________What are the data?_____________________
% ordre de rachel et vinca systematique:
% 'Obj' 'Obj' 'NO' 'NO' 'Obj'
% Sessions as they are organized by Rachel:
CondNam={'Obj' 'Obj' 'NO' 'NO' 'Obj'};
Circu=2*pi*15; %wheel circumference
%Circu=94.25; %from rachel's tape measure
cd (AnlPath);
nospike=0; %put 1 if you want to do bhv only CHANGE AS WANTED
isoe=0; %create a variable to classify as open ephys format (=1) or not (=0, means .dat binary format) DO NOT MODIFY
A=dir;
oebpath=pwd;

%% basic data extraction
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

%% Data extraction :  Continuous Open Ephys file format 
if ~ispheno
        %____________________preparation_____________________
    if startsWith(A(3).name,'10');  %NEW !! haha trop contente de cette idée : les 3 formats commence par 100 donc valide tous les 100_trucmuch mais comme j'ai des sessions avec 101_ADC et pas 100 au moins ca marchera avec elles aussi ! 
        isoe=1; 
        if strcmp( A(3).name(4:5),'_1') %NEW j'ai changé pour que ça marche avec tous les processor id voir premier if de cette section. 
            prefix =  A(3).name(1:4);
            fff(1)=65; %XposVR
            fff(2)=66;% transitions in Y -back and fourth trip 
            fff(3)=67;% Xpos real
            fff(4)=68;%Sess??
        elseif  strcmp(A(3).name(4:7),'_CH1') %NEW j'ai changé pour que ça marche avec tous les processor id voir premier if de cette section. 
            prefix =  A(3).name(1:6);
            fff(1)=65;
            fff(2)=66;
            fff(3)=67;% Xpos real
            fff(4)=68;
        elseif strcmp(A(3).name(4:8),'_ADC1')
            prefix = A(3).name(1:7);
            fff(1)=1;
            fff(2)=2;
            fff(3)=3;% Xpos real
            fff(4)=4;

        end
    end
    %____________________file opening_____________________
    if isoe %encore plus simple et rigide je pense mais je peux me tromper
        tic
        Nam=[prefix,num2str(fff(2)),'.continuous'];
        [transition, t, info] = fct_read_continuous_file(fullfile(AnlPath, Nam)); %NEW
        sf1=info.sampleRate; %NEW1
        sf2 =sf1/1000; 
        transition=downsample(transition,sf2);%NEW2 (automatiser au cas où on ai pas le meme SR pour pas tout casser)
        t=t(1:sf2:end); %NEW 
        Y=transition;
        %  plot(t,transition)

        %_____________________position data_____________________

        Nam=[prefix,num2str(fff(1)),'.continuous'];
        [X, tx, info] = fct_read_continuous_file(fullfile(AnlPath, Nam)); %NEW
        FirstTime=tx(1);
        X=downsample(X,sf2); %was 30, now auto


        % plot(tx, X)

        %_____________________wheel data_____________________

        Nam=[prefix,num2str(fff(3)),'.continuous'];
        [Whl, twl, ~] = fct_read_continuous_file(fullfile(AnlPath, Nam)); %NEW
        %
        Whl=downsample(Whl,25);
        twl=twl(1:25:end);


        Nam=[prefix,num2str(fff(4)),'.continuous'];
        [Sess, tses, info] = fct_read_continuous_file(fullfile(AnlPath, Nam)); %NEW
        Sess=downsample(Sess,25);
        
        
        tx=tx(1:sf2:end);
        sf=sf1/sf2; %so 1000Hz...
        si=(1/sf1)*1E6; % in microseconds

    end
    toc
    clear Nam 

%% Data extraction : Continuous.dat Binary file format  

%_____________________initialize thanks to oebin_____________________

    if  ~isoe

            oebpath=pwd; %moved


            fname = 'structure.oebin';
            fid = fopen(fname);
            raw = fread(fid,inf);
            str = char(raw');
            fclose(fid);
            val = jsondecode(str);
            header=val.continuous(1);
            sr = val.continuous.sample_rate;
            idADC=  strfind(arrayfun(@(x) startsWith(val.continuous.channels(x).channel_name,'ADC') , 1:val.continuous.num_channels),1);  %ça parait compliqué mais c'est l'équivalent d'une simple compréhension de liste de python sur matlab je pense;


    %_____________________localise .dat file and read bhv_____________________
        tic
        if isfile('continuous.dat')
            %path_dat=pwd    not used in this case for now :)
            [data, t]=Read_OEP_Binary ('continuous.dat',[idADC(1:4)], 0,-1, sr/1000);
            top=t;
        elseif ~isfile('continuous.dat')  %equivalent à si continuous.dat est ailleurs (pas sorti de rythmFGPA) 
            path_dat = strcat(AnlPath, filesep, 'continuous\Rhythm_FPGA-', num2str(val.continuous.source_processor_id), '.', num2str(val.continuous.source_processor_sub_idx), filesep); %va chercher la localisation du .dat dans les sous-dossiers
            [data, t]=Read_OEP_Binary_altered(strcat(path_dat,'continuous.dat'),idADC(1:4), 0,-1, sr/1000); %la version "altérée consiste en la modification d'une ligne (D.timestamps) où j'ai ajouté le path_dat (comme au dessus) avant continuous.dat dans les inputs de ReadNPY pour qu'il trouve le fichier quin'est pas au même endroit que oebin ! c'est tout!

            if ~isfile(strcat(path_dat, 'continuous.dat'))
                error('file continuous.dat not find, please check location')
            end

        end

        if length(t)~=(size(data,2))
            disp('OUPS TIMESTAMPS ARE FXXXD!!! REPAIRING')
            t=linspace(t(1), t(1)+ size(data,2)/1000,size(data,2));
        end

        toc

    %_____________________make bhv variables_____________________

        transition=data(2,:);
        X=data(1,:);
        Y=data(2,:);
        Whl=data(3,:);
        Sess=data(4,:);

        %!!!! attention différent de Julie ici, voir BasicBehavFromRec ligne 93 to 106 !!!!!

        si=(t(2)-t(1)); % in microseconds
        sf=1/si;
        si=round(si*1E6);
        FirstTime=t(1); %SI PROBLEME PREMIERE SESSION PERSISTE CREER BOUCLE IF POUR FIRST TIME OU SI INF A x ALORS PRENDRE SECONDE TIME
    end

%%  Vérification des variables visuellement 
figure
subplot(3,1,1)
plot(t- FirstTime,X)
hold on

subplot(3,1,2)
plot(t- FirstTime,Y)
hold on
subplot(3,1,3)
plot(t- FirstTime,  Sess)

%% je rajoute une demande d'input pour simplifier car trouver les bons paramètres est un casse tête sinon
% isfirstsess = 'Is the first detected session change a real change (Y) or the onset of the env (N)? [Y]';
% txt = input(isfirstsess,'s');
% if isempty(txt)
%     isfirstsess=1;
% elseif strcmp(txt,'Y')
%     isfirstsess=1;
% else
%     isfirstsess = 0;
% end
%y = x*10
%% Compute behavior info 
%%%%Start at 0

    if min(X)<0
        X=X+abs(min(X));
    elseif min(X)>0
        X=X-(min(X));
    end

    X1=X;
    %
    X=(X./max(X))*145; % to fit what the animal runs for real, like if gain=1


    %%Start at 0
    if min(Whl)<0
        Whl=Whl+(abs(min(Whl)));
    elseif min(Whl)>0
        Whl=Whl-((min(Whl)));
    end
    %

    if max(Whl)>4
        Turn=2*pi*15;%in cm
        VoltTurn = 1725.75;%in whl units estimated by 10 real turns of wheel by julie on April 2024
        Whl=(Whl.*Turn)./VoltTurn;
        transitionThresh=2600;%%

    %     highYlevel=2900; % for info, to get rid of middle y values at the start of a session
    else
        Turn=2*pi*15;%in cm
        VoltTurn=1.3435;%in whl units estimated by 50 real turns of wheel by julie April 2024
        Whl=(Whl.*Turn)./VoltTurn;
        %transitionThresh=2;%OLD pose pb avec détection resets (sur-détecté avec tT=2)
        transitionThresh=2;
    %     highYlevel=3;% voir avec un enregistrement sans rec?
    end

    % detect when wheel data jumps
    % so that I patch it back together after

    thresh=5;% can't decently run 5cm in 1ms

    jumpsup=find(diff(Whl)>thresh);
    jumpsdown=find(diff(Whl)<-thresh);
    Whl2=Whl;
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

        Whl=Whl2.*-1;%it was going dozn I need it to go up
    else
        disp('BAAD ISSUE WITH Wheel acquisition!!!')

    end

    Whl=Whl+abs(min(Whl));% put it to zero

    Whl=Whl(:);
    t=t(:);
    X=X(:);
    Sess=Sess(:);


    %% get separation between sessions S1-S2-S3

    %%%%%%%%%ATTENTION CHANGEMEN%%%%%%%%%%%
    %SessStarts=t(LocalMinimaPP(Sess, 1000,500)); OLD JULIE !!
    %[SessStarts,b] = max(Sess); %%%%pour le nouvel ordi 
    if isoe 
        SessStarts=t(LocalMinimaPP(Sess, 1000,0.05));  %POUR CONTINUOUS
        sessS=LocalMinimaPP(Sess, 1000,0.05);
    else
        SessStarts=t(LocalMinimaPP(Sess, 1000,500)); %OLD JULIE !!
        sessS=LocalMinimaPP(Sess, 1000,500);
    end

    % sessS=[1; sessS];
    %I use this to remove first session transition if it is too close from  0
    if ~isempty(sessS)
        if sessS(1)<(sf*30) %if it's less than 30 seconds from start... I could use more
            SessStarts(1)=[]; 
            SessStarts=[t(sessS(1));SessStarts; t(end)];
             sessS(1)=[];
    %         
        else
            SessStarts=[t(1);SessStarts; t(end)];
        end
    else
        %one condition
        SessStarts=[t(1); t(end)];
        disp('Session detection ISSUE Check your cables!!!!!!!')
    end



    transition2=transition./(max(transition));% here is a potential issue 

    trialstarts=find(abs(diff(transition2))>0.03);
    remove=find(diff(trialstarts)==1);%
    trialstarts(remove+1)=[];
    % keep=ones(size(trialstarts));
    % 
    % 
    % firstTrajSess=find(transition(trialstarts)>lowYlevel & transition(trialstarts)<highYlevel);
    % keep(firstTrajSess(2:end)-1)==0;% remove the trialstart before the jump






    %% extract speed from wheel distance

    Whl2=smooth(Whl,300);
    Whl2=resample(Whl2,1,300);
    RealSpeed=(diff(Whl2).*1000)./300; % estimate speed in 300ms windows but in cm/sec
    RealSpeed(end-10:end)=0;
    RealSpeed(find(RealSpeed<=0))=0;
    Speed=interp1(t(1:300:end-300),RealSpeed,t); %put back at 1kHz...


    %% extract speed from VR distance
    % le souci c'est qu'il est pas teleporté à 0 mais à 6 ou 5
    %here I am unwrapping VR distance
    jumpsup=find(abs(diff(X))>50);
    if ~isempty(jumpsup)
        X2=X;
        for jj=1:length(jumpsup)
            dist2add=abs(X(jumpsup(jj)+2)-X(jumpsup(jj)));
            X2(jumpsup(jj)+1:end)=X2(jumpsup(jj)+1:end)+dist2add;% here sometimes you get 2 samples that go through transition... annoying
            X2=X2-abs(min(X2));
        end
    else
        disp('BAAD ISSUE WITH VR DETECTION')

    end
     X3=X2;
    jumpsup=find(abs(diff(X2))>50);
    if ~isempty(jumpsup)
    %     X3=X2;
        for jj=1:length(jumpsup)
            dist2remove=abs(X2(jumpsup(jj)+2)-X2(jumpsup(jj)));
            X3(jumpsup(jj)+1:end)=X3(jumpsup(jj)+1:end)-dist2remove;% here sometimes you get 2 samples that go through transition... annoying
    %         X3=X3-abs(min(X2));
        end

    end



     X2=X3(1:300:end);
    % X2=resample(X3,1,300);
    X2Speed=(diff(X2).*1000)./300; % estimate speed in 300ms windows but in cm/sec
    % X2Speed(end-300:end)=0;
    X2Speed(find(X2Speed<=0))=0;% this is questionnable (when ball goes backwards
    XSpeed=interp1(t(1:300:end-300),X2Speed,t); %put back at 1kHz...
    % XSpeed=X2Speed;
    starts=trialstarts;

     save('Phenosys.mat','t','X','transition','Sess', 'Whl','starts','sessS','SessStarts','Speed','XSpeed','transitionThresh','t')
    %save([savepath '\' Namesession '_Phenosys.mat'],'t','X','transition','Sess', 'Whl','starts','sessS','SessStarts','Speed','XSpeed','transitionThresh')

end

%% preparing a structure for behavior outputs
if ~istraj
  sf=1000; % should be for rachel
% W=Whl;

cd (AnlPath)
if exist( 'TrajData.mat','file')
    
    answer = questdlg('Ovewrite TrajData file?','TrajData analysis',...
        'Yes','No','No');
    switch answer
        case 'Yes'
            
            [Traj, newW]= runTraj (t, X, transition, Whl, starts, sessS, Speed, XSpeed,transitionThresh,AnlPath, OutputFolder, Namesession)  ;
            
        case 'No'
            cd (AnlPath)
            disp('Using previously saved Phenosys.mat file')
            load ('TrajData.mat')
    end
else
    
    
    [Traj, newW]= runTraj (t, X, transition, Whl, starts, sessS, Speed, XSpeed,transitionThresh,AnlPath, OutputFolder,Namesession);
    %newW is reset at zero for each trajectory
end

%% for figures
if plotting
    edgesX=[0:2.5:max(X)];
    edgesW=[0:5:max(newW)];
    plotmybehav (Traj,edgesX,edgesW, OutputFolder,Namesession,CondNam);
    cd (AnlPath)
end
end
%%

%Number of condition
nbsess=length(unique([Traj.Cond])); %max should be 5 based on the CondNam but things can happen and make the recording shorter

%Moment of sessions in sample
mysessS=[1; sessS; length(t)]; 

%get information for cells and clusters
[Spikes,phyData,goodClusListIdx,datacell,headercell,clusters_info_new]=getphyInfo(CircusPath); %J'ai modifié pour bien récupérer les infos du tsv car ça buguait chez moi.

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
    CID(cl)=str2double(cell2mat(clusters(cl,1)));
end


CID=CID(goodClusListIdx); %keeping only the good ones (as saved in Phy -'good')
%%
allcel.id_cel = CID;
allcel.itime_spk = Spikes.SpikeTimes;
allcel.time_spk = Spikes.SpikeTimes/25000;
allcel.time_spk2 = Spikes.SpikeTimes/25
allcel.id_spk = Spikes.Units;
allcel.Spikes = Spikes;

%OLD : 
%waveform = phyData.templates(CID,:,:); 
%waveform = phyData.templates(goodClusListIdx,:,:); 

%NEW : 
bestchannels=str2double([cluster_infos_tsv(2:end,3)]);
bestchannelsunits=str2double([cluster_infos_tsv(2:end,1)]);
bestchannelsunitsgood=[bestchannelsunits(strcmp(cluster_infos_tsv(2:end,6),'good'))];
bestchannelsgood=[bestchannels(strcmp(cluster_infos_tsv(2:end,6),'good'))];
waveformnew=struct();
cd(AnlPath)
    waveformnew.meantemp=[]
allcel.waveform=nan(50,21,length(allcel.id_cel));
allcel.waveform_edited=nan(50,10,length(allcel.id_cel));
for i=1:length(unique(Spikes.Units));
    idx =bestchannelsunitsgood(i);
    channel=bestchannelsgood(i);% logical indexing 
    mask=(double(Spikes.Units)==idx);
    selectedtimespk= allcel.time_spk2(mask,1);
    waveformnew(i).unit=idx;
    waveformnew(i).bestchan=channel;

    if isoe %openephys format %REGLER CA IL FAUT QUE ISOE SOIT JUSTE MEME SI PAS APPELE
        Nam=[prefix(1:end-3),'CH',num2str(channel),'.continuous'];
        [channeldata, ~, info] = fct_read_continuous_file(fullfile(AnlPath, Nam)); 
    else  %binary .dat
        if isfile('continuous.dat')
            [channeldata, t]=Read_OEP_Binary ('continuous.dat',channel, 0,-1, sr/1000);
            top=t;
        elseif ~isfile('continuous.dat')  %equivalent à si continuous.dat est ailleurs (pas sorti de rythmFGPA) ; ne marchera surement pas car val n'exisste pas si on a déja fait traj etc. j'ai une soluce en tête si besoin :) 
            path_dat = strcat(AnlPath, filesep, 'continuous\Rhythm_FPGA-', num2str(val.continuous.source_processor_id), '.', num2str(val.continuous.source_processor_sub_idx), filesep); %va chercher la localisation du .dat dans les sous-dossiers
            [channeldata, t]=Read_OEP_Binary_altered(strcat(path_dat,'continuous.dat'),channel, 0,-1, sr/1000); %la version "altérée consiste en la modification d'une ligne (D.timestamps) où j'ai ajouté le path_dat (comme au dessus) avant continuous.dat dans les inputs de ReadNPY pour qu'il trouve le fichier quin'est pas au même endroit que oebin ! c'est tout!
        end
    end
% %     disp(i)
    selectedamptimespk=selectedtimespk(1:50);
    for j=1:50 %j = nb de spike sélectionnés pour le sampling
        timewindow=[-10:10]; %fenetre pour le spike 
% %         disp(j)
       timewindow=timewindow+double(selectedamptimespk(j));
       if timewindow(1) <= 0
           timewindow=timewindow+abs(timewindow(1))+1;
       end
           waveformnew(i).waveform(:,j)=channeldata(timewindow);
%              waveformnew(i).meantemp=nan(1,length(timewindow)); 
             allcel.waveform(j,1:length(timewindow),i) =waveformnew(i).waveform(:,j);

    end
    end
    %% select spike in the time window: 
    for i = 1:length(unique(Spikes.Units))_
    figure
    
                 plot(allcel.waveform(j,1:length(timewindow),i) );
             prompt = char("Choisissez les 10 timestamps dans lesquels on a le spike et les écrire sous forme 2:11 svp");
                x = input(prompt)
              allcel.waveform_edited(:,:,i)=allcel.waveform(:,x,i)
              close all
    for y=1:length(timewindow);
    waveformmoy(i,1:10)=[mean(allcel.waveform_edited(:,:,i))]; 
    end
   
end
% 
% allcel.waveform=cat(3,nan(8),nan(8),nan(8));
% %new
% allcel.waveform(length(timewindow),j,i) =waveformnew(1:end).waveform(1:end);
% allcel.waveform_circus_mean =  struct('meantemp', {waveformnew(1:end).meantemp});

prompt = "Choisissez les 10 timestamps dans lesquels on a le spike et les écrire sous forme 2:12 svp";
x = input(prompt)
y = x*10
save([OutputFolder '\' Namesession '_ePhy.mat'],'allcel')
save([OutputFolder '\' Namesession '_TrajData.mat'],'Traj','newW')
save([OutputFolder '\' Namesession '_Phenosys.mat'],'t','X','transition','Sess', 'Whl','starts','sessS','SessStarts','Speed','XSpeed','transitionThresh')

%EdgesX=[0:1:max(X)];%preparing bins for X and ball
%EdgesBall=[0:2:max(newW)];


%% Do ratemap and place fields
if ~isplacecell 
    EssaiRmap_G
end
