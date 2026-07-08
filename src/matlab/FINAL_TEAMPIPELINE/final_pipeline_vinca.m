% === Paramètres de base à modifier à chaque lancement ===
session_path='X:\MATLAB\Data\Data_raw\VS137\VS137_2026-01-18_15-39-11\Record Node 104\experiment1\recording1\continuous\Rhythm_FPGA-101.0\';
whatprob =5;                           % Identifiant de la sonde 5, 64 ou 77 

%% === Paramètres de base à modifier pour adapter le code à vos data 

Rawfolder='X:\MATLAB\Data\Data_raw\';  %root
session_path=session_path(size(Rawfolder,2)+1:end);
OutputFolder = 'X:\MATLAB\Data\All_data_final\process';      % Dossier de sortie
AnlPath = strcat(Rawfolder,session_path);

if isstrprop(session_path(5), 'digit') 
    Namesession = session_path(7:31);       % VS100 et +
    Mouse=Namesession(3:5);
else
    Namesession = session_path(6:29);   %VS99 et -
    Mouse=Namesession(3:4);
end

if str2num(Mouse) <57 %NE MARCHE QUE ¨POUR VS !! 
    Circu=43; %old wheel circonf
else 
    Circu=2*pi*15; %ball circumference
end
%%
if contains(session_path, 'continuous', 'IgnoreCase', false)
    CircusPath = strcat(Rawfolder,session_path,'\continuous\continuous-merged.GUI');
    isoe=false;
    prompt='_DAT_';
    format='continuous';
else
    CircusPath = strcat(Rawfolder,session_path,'\Continuous_Data\Continuous_Data-merged.GUI');
    isoe=true;
    prompt='_OPENEPHYS_';
    format='Continuous_Data';
end

cd(AnlPath)

prompt=strcat('Analysing session is___ ', Namesession, '____with electrode____', num2str(whatprob),'____in the format_',prompt);
disp(prompt)

%% Phenosys 
exists=exist(fullfile(strcat(AnlPath,'\process'), 'Phenosys.mat'), 'file') == 2;
if exists 
    prompt='Pheno already exist, do you want to redo Phenosys?';
    if askConfirmation(prompt)
        Phenosys(AnlPath,OutputFolder,isoe,Circu,Namesession)
        load (fullfile(strcat(AnlPath,'\process'), 'Phenosys.mat'))
        disp('Phenosys.mat done')
    else
        load (fullfile(strcat(AnlPath,'\process'), 'Phenosys.mat'))
        disp('loading Phenosys.mat')
    end
else 
    disp('Doing Phenosys on session')
    Phenosys(AnlPath,OutputFolder,isoe,Circu,Namesession)
    load (fullfile(strcat(AnlPath,'\process'), 'Phenosys.mat'))
    disp('Phenosys.mat done')
end

%% TrajData
% if ~exist('WhlSpeed', 'var'), WhlSpeed= Speed; end %FOR CODING PURPOSE ONLY ENLEVER APRES DEBUGGAGE
% if ~exist('CondNam', 'var')  %FOR CODING PURPOSE ONLY ENLEVER APRES DEBUGGAGE
%     NSess=length(SessStarts)
%      NSess=length(sessStartsIdx)
%     prompts=cell(1,NSess);
%     ; %FOR CODING PURPOSE ONLY ENLEVER APRES DEBUGGAGE
%     prompts(1:end)=compose('Name for Sess %d', 1:NSess); %FOR CODING PURPOSE ONLY ENLEVER APRES DEBUGGAGE
%     Defaults=cell(1,NSess); %FOR CODING PURPOSE ONLY ENLEVER APRES DEBUGGAGE
%     Defaults(1:end)={'PO'}; %FOR CODING PURPOSE ONLY ENLEVER APRES DEBUGGAGE
%     CondNam=inputdlg(prompts,'Names',1,Defaults); %FOR CODING PURPOSE ONLY ENLEVER APRES DEBUGGAGE
% end %FOR CODING PURPOSE ONLY ENLEVER APRES DEBUGGAGE


exists=exist(fullfile(strcat(AnlPath,'\process'), 'TrajData.mat'), 'file') == 2;

if exists 
    prompt='Traj already exist, do you want to redo TrajData?';
    if askConfirmation(prompt)
         [Traj, newW]= runTraj_vinca(t_ds, X_ds_n, Y_ds, Whl_r, starts, stops, sessStartsIdx, WhlSpeed, XSpeed,transitionThresh,AnlPath, OutputFolder, Namesession,CondNam,Circu) ;
         load (fullfile(strcat(AnlPath,'\process'), 'TrajData.mat'))
        disp('TrajData.mat done')
    else
        load (fullfile(strcat(AnlPath,'\process'), 'TrajData.mat'))
        disp('loading TrajData.mat')
    end
else 
    disp('Doing Traj on session')
     [Traj, newW]= runTraj_vinca (t_ds, X_ds_n, Y_ds, Whl_r, starts, stops, sessStartsIdx, WhlSpeed, XSpeed,transitionThresh,AnlPath, OutputFolder, Namesession,CondNam,Circu)  ;
     load (fullfile(strcat(AnlPath,'\process'), 'TrajData.mat'))
    disp('TrajData.mat done')
end
%%
NSess=length(unique([Traj.Cond]));
exists=exist(fullfile(strcat(AnlPath,'\process'),'ePhy.mat'));
if exists
    prompt = 'ePhy already done, do you want to redo it?';
    if askConfirmation(prompt)
        [Spikes,phyData,goodClusListIdx,datacell,headercell,clusters_info_new]=getphyInfo(CircusPath); 
        if size(clusters_info_new,2)==1
            disp('reshaping csv')
            
            cluster_infos_tsv=[];
            for i = 1:9:length(clusters_info_new);
                cluster_infos_tsv=[cluster_infos_tsv;[clusters_info_new(i:i+8)].'];
               
            end
        xlswrite( 'clusters_info.xlsx',cluster_infos_tsv);    
        end
        [allcel]=make_ephy(cluster_infos_tsv,datacell,phyData,Spikes,goodClusListIdx,Namesession,AnlPath,whatprob,isoe,OutputFolder)
    end
else
    disp('ePhy doesnt exist, doing ePhy now')
    
        [Spikes,phyData,goodClusListIdx,datacell,headercell,clusters_info_new]=getphyInfo(CircusPath); 
        if size(clusters_info_new,2)==1
            disp('reshaping csv')
            
            cluster_infos_tsv=[];
            for i = 1:9:length(clusters_info_new);
                cluster_infos_tsv=[cluster_infos_tsv;[clusters_info_new(i:i+8)].'];
               
            end
        xlswrite( 'clusters_info.xlsx',cluster_infos_tsv);    
        end
        [allcel]=make_ephy(cluster_infos_tsv,datacell,phyData,Spikes,goodClusListIdx,Namesession,AnlPath,whatprob,isoe,OutputFolder)

end
  load (fullfile(strcat(AnlPath,'\process'), 'ePhy.mat'))
    disp('ePhy.mat loaded')
    
 
Rmap_G

%% === Fonctions Majeures === 
function Phenodone = Phenosys(AnlPath, OutputFolder, isoe,circu,Namesession)
          cd (AnlPath);
          A=dir;
          idx = find(startsWith({A.name}, '10'), 1);
         if isempty(idx), idx = 4; end
         
% OPENEPHYS FORMAT
    if startsWith(A(idx).name,'10'); 
        isoe=1; 
        if strcmp( A(idx).name(4:5),'_1') 
            prefix =  A(idx).name(1:4);
            fff(1)=65; %XposVR
            fff(2)=66;% transitions in Y -back and fourth trip 
            fff(3)=67;% Xpos real
            fff(4)=68;%Sess??
        elseif  strcmp(A(idx).name(4:7),'_CH1') 
            prefix =  A(idx).name(1:6);
            fff(1)=65;
            fff(2)=66;
            fff(3)=67;% Xpos real
            fff(4)=68;
        elseif strcmp(A(idx).name(4:8),'_ADC1')
            prefix = A(idx).name(1:7);
            fff(1)=1;
            fff(2)=2;
            fff(3)=3;% Xpos real
            fff(4)=4;
        end
    end


    %____________________file opening_____________________
    if isoe 
        tic
        Nam=[prefix,num2str(fff(2)),'.continuous'];
        [Y, ty, info] = fct_read_continuous_file(fullfile(AnlPath, Nam)); 
        Nam=[prefix,num2str(fff(1)),'.continuous'];
        [X, tx, info] = fct_read_continuous_file(fullfile(AnlPath, Nam)); 
        Nam=[prefix,num2str(fff(3)),'.continuous'];
        [Whl_r, twl, ~] = fct_read_continuous_file(fullfile(AnlPath, Nam)); 
        Nam=[prefix,num2str(fff(4)),'.continuous'];
        [Y_r, tses, info] = fct_read_continuous_file(fullfile(AnlPath, Nam)); 
        clear Nam 
        
        sf1=info.sampleRate; 
        sf2 =sf1/1000; 
        sf=sf1/sf2; %so 1000Hz...
        si=(1/sf1)*1E6; % in microsecond
        FirstTime=tx(1);
        
        Y_ds=downsample(Y,sf2);
        t_y_ds=ty(1:sf2:end); 
        X_ds=downsample(X,sf2); 
        tx_ds=tx(1:sf2:end);
        Whl_ds=downsample(Whl_r,sf2);
        t_wl_ds=twl(1:sf2:end);       
        Y_r_ds=downsample(Y_r,sf2);
        t_ds=tx_ds; 
        
        elapsedTime = toc;
        disp(strcat('Behavior loaded in_ ',num2str(elapsedTime),'seconds'))
    
    end
   
%CONTINUOUS.DAT FORMAT
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

           elapsedTime = toc;
           disp(strcat('Behavior loaded in_ ',num2str(elapsedTime),'seconds'))

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
    
    if circu == 43 %WHEEL
        endVR=200;
         X_ds_n=(X_ds_n./max(X_ds_n))*endVR; 
    elseif circu == 2*pi*15 %BALL
        endVR=145;
         X_ds_n=(X_ds_n./max(X_ds_n))*endVR; 
    end
        
   % to fit what the animal runs for real, like if gain=1

    %%Start at 0
    if min(Whl_ds)<0
        Whl_ds_n=Whl_ds+(abs(min(Whl_ds)));
    elseif min(Whl_ds)>=0
        Whl_ds_n=Whl_ds-((min(Whl_ds)));
    end
    %

        
%ESSAIE RECalcul du VoltTurn
jumpX=find(diff(X_ds_n)<-130); %Select teleportations
Xchosen=X_ds_n((jumpX(2)+1):(jumpX(3)-1)); %Prends 2eme trajectoire VR dont on connait longueur
endVRidx=find(Xchosen==max(Xchosen),1,'first')+(jumpX(2)+1);
trajVRidx=[jumpX(2)+1,endVRidx];

Whlchosen=Whl_ds_n(trajVRidx(1):trajVRidx(2))
jumpWhl=find(diff(Whlchosen)>2*std(Whlchosen),1,'first')
if ~isempty(jumpWhl)
    Whldiff(jumpWhl) = nan;
end

Whldiff=diff(Whlchosen);
Whldiff(jumpWhl)=nan;
Whldist=abs(nansum(Whldiff));

deltaW=Whldist;
distance_en_X=endVR;
   if circu == 43 %WHEEL
VoltTurnWheel= deltaW*circu / distance_en_X ;
    elseif circu == 2*pi*15 %BALL
VoltTurnWheel= deltaW*circu / distance_en_X ;
    end
        
% VoltTurnWheel= deltaW*100 / distance_en_X ;
Whl_r=(Whl_ds_n.*circu)./VoltTurnWheel;
    Whlverif=Whl_r(trajVRidx(1):trajVRidx(2));
plot(Whlverif)
Whlverif(jumpWhl)=nan;
% prompt=strcat('confirm that the total distance of ', num2str(abs(nansum(diff(Whlverif)))), 'is close to the expected one of_', num2str(endVR));
% disp(prompt)
% 

  transitionThresh=mean (Y_ds);
    % detect when wheel data jumps
    % so that I patch it back together after

    thresh=10;% can't decently run 5cm in 1ms

    jumpsup=find(diff(Whl_r)>thresh);
    jumpsdown=find(diff(Whl_r)<-thresh);
    Whl_u=Whl_r;
    if ~isempty(jumpsup)
        for jj=1:length(jumpsup)
            dist2remove=abs(Whl_u(jumpsup(jj)+1)-Whl_u(jumpsup(jj)));
            Whl_u(jumpsup(jj)+1:end)=Whl_u(jumpsup(jj)+1:end)-dist2remove;
            Whl_u=Whl_u+abs(min(Whl_u));
        end

        for jj=1:length(jumpsdown)
            dist2add= abs(Whl_u(jumpsdown(jj)+1)-Whl_u(jumpsdown(jj)));
            Whl_u(jumpsdown(jj)+1:end)=Whl_u(jumpsdown(jj)+1:end)+dist2add;
        end

        Whl_r=Whl_u.*-1;%it was going down I need it to go up
    else
        disp('BAAD ISSUE WITH Wheel acquisition!!!')

    end

    Whl_r=Whl_r+abs(min(Whl_r));% put it to zero
    Whl_r=Whl_r(:);
    t_ds=t_ds(:);
    X_ds_n=X_ds_n(:);
    Y_r_ds=Y_r_ds(:);

%% get separation between sessions S1-S2-S3

hf=figure;
ax2=subplot(2,1,2)

plot(X_ds,'k')
ax1=subplot(2,1,1);

plot(Y_ds,'k')
linkaxes([ax1 ax2], 'x');
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
sessStartsIdx=zeros(NSess,1);

t2=t_ds-t_ds(1);

  figure;
for c = 1:NSess  % trying to plot a piece at a time, for precision

    % Sélection des indices valides
    idx = find(t2 > times(c,1) & t2 < times(c,1) + max(t2));

    % --- Subplot pour X_ds (non interactif)
    subplot(2,1,1);
    plot(t2(idx), X_ds(idx), 'b');
    title(['X\_ds - Session ', num2str(c), 'pour comparaison']);
    xlabel('Temps');
    ylabel('X');
    grid on;

    % --- Subplot pour Y_ds (interactif)
    subplot(2,1,2);
    plot(t2(idx), Y_ds(idx), 'r');
    hold on;
    title(['Click on Session ', num2str(c,'%d'), ' precise start HERE']);
    xlabel('Temps');
    ylabel('Y');
    grid on;

    % --- GINPUT : Début session
    [times(c,1), y, ~] = ginput(1);
    plot(times(c,1), y, 'or');

    % --- GINPUT : Fin session
    title(['Click on Session ', num2str(c,'%d'), ' stop, before the drop for all sess']);
    [times(c,2), y, button] = ginput(1);
    plot(times(c,2), y, 'or');

    % Préparer la session suivante
    if c < NSess
        times(c+1,1) = times(c,2);
    end

    hold off;
    pause(1);

    % Enregistrer l'index de début
    sessStartsIdx(c) = find(t2 >= times(c,1), 1, 'first');
end


SessStartsTimes=times(:,1);

X_ds_n=X_ds_n(sessStartsIdx(1):end);
X_ds_n= (X_ds_n - min(X_ds_n)) / (max(X_ds_n)- min(X_ds_n)) * endVR;
t_ds=t_ds(sessStartsIdx(1):end);
Whl_r=Whl_r(sessStartsIdx(1):end);
Y_r_ds=Y_r_ds(sessStartsIdx(1):end);
    if ~isempty(sessStartsIdx)
    else
        %one condition
        SessStartsTimes=[t_ds(1); t_ds(end)];
        disp('Session detection ISSUE Check your cables!!!!!!!')
    end
%     bhvboundariestimes=[times(1,1),times(end,end)]
if length(sessStartsIdx) ==1 
    boundaridx=[find(t2 >=times(1,1), 1, 'first'), find(t2 >=times(end,end), 1, 'first')];
    trialstarts = find(abs(diff(Y_ds( boundaridx(1,1):boundaridx(1,2) )))>(max(diff(Y_ds( boundaridx(1,1):boundaridx(1,2) )))./2));
    trialstarts=[boundaridx(1),trialstarts+boundaridx(1)];
    trialstarts=sort(trialstarts)
    trialstops=[trialstarts(2:end),boundaridx(end)];
     hold off
    plot(diff(Y_ds))
    hold on 
     plot(trialstarts,zeros(length(trialstarts)),'o')
    if      askConfirmation('do you confirm right detection of trials?')
        disp('lets continue analysing')
    else
        return
    end
elseif length(sessStartsIdx)>1
    for idxsess = [1:length(sessStartsIdx)]
        boundaridxtemp=[find(t2 >=times(idxsess,1), 1, 'first'), find(t2 >=times(idxsess,2), 1, 'first')];
        boundaridx=[find(t2 >=times(1,1), 1, 'first'), find(t2 >=times(end,end), 1, 'first')];
        trialstartsess = find(abs(diff(Y_ds( boundaridxtemp(1,1):boundaridxtemp(1,2) )))>(max(diff(Y_ds( boundaridxtemp(1,1):boundaridxtemp(1,2) )))./2));
        sessIDS=find(abs(diff(Y_ds( boundaridx(1,1):boundaridx(1,2) )))  > max(diff(Y_ds( boundaridxtemp(1,1):boundaridxtemp(1,2) ))));
        figure(17)
        plot(Y_ds(( boundaridxtemp(1,1):boundaridxtemp(1,2))))
        if idxsess==1, 
            startid=boundaridxtemp(1) ;
            trialstarts=[startid,trialstartsess+startid] ;
            trialstops=[trialstarts(2:end),boundaridxtemp(end)];
        elseif any(strcmp(CondNam,'POMb'))
         idxMb=[find(strcmp(CondNam,'POMb'))].';
             for i = [idxMb]
                 boundaridxtemp=[find(t2 >=times(i,1), 1, 'first'), find(t2 >=times(i,2), 1, 'first')];
                 boundaridx=[find(t2 >=times(1,1), 1, 'first'), find(t2 >=times(end,end), 1, 'first')];
                 trialstartsess=find(abs(diff(X_ds( boundaridxtemp(1,1):boundaridxtemp(1,2) )))>(max(abs(diff(X_ds( boundaridxtemp(1,1):boundaridxtemp(1,2) ))))./2))+boundaridxtemp(1,1);
                 trialstarts =[trialstarts ,  trialstartsess];
                 trialstops=[trialstops,trialstartsess(2:end),boundaridxtemp(end)];
                 sessIDS=find(abs(diff(Y_ds( boundaridx(1,1):boundaridx(1,2) )))  > max(diff(Y_ds( boundaridxtemp(1,1):boundaridxtemp(1,2) ))));
             end
        else
            startid=boundaridxtemp(1) ;
         trialstarts=[startid,trialstarts,trialstartsess+startid]    ;
         trialstarts=sort(trialstarts)
         trialstops=[trialstops,trialstartsess+startid,boundaridxtemp(end)]
        end
%          trialstarts(diff(trialstarts)<=5000)=[]
    end
    
end
     hold off
    plot(diff(Y_ds))
    hold on 
     plot(trialstarts,zeros(length(trialstarts)),'o')
     plot(trialstops,zeros(length(trialstarts)),'o')
    if      askConfirmation('do you confirm right detection of trials?')
        disp('lets continue analysing')
    else
        return
    end

%     OLD VEE  : trialstarts=find(abs(diff(X_ds_n(boundaridx(1):boundaridx(1,2))))>50); %OLD (transition2))>0.03
    

    %% extract speed from wheel distance

    Whl_u=smooth(Whl_r,300);
    Whl_u=resample(Whl_u,1,300);
    RealSpeed=(diff(Whl_u).*1000)./300; % estimate speed in 300ms windows but in cm/sec
    RealSpeed(end-10:end)=0;
    RealSpeed(find(RealSpeed<=0))=0;
    WhlSpeed=interp1(t_ds(1:300:end-300),RealSpeed,t_ds); %put back at 1kHz...


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
    starts=sort(trialstarts);
    stops=sort(trialstops);

    if ~exist(strcat([AnlPath '\process\'] ))
        mkdir(strcat([AnlPath '\process\'] ))
    end
     save( [AnlPath '\process\' 'Phenosys.mat'],'t_ds','X_ds','X_ds_n','Y_ds','Y_r_ds', 'Whl_r','starts','stops','sessStartsIdx','SessStartsTimes','WhlSpeed','XSpeed','transitionThresh','CondNam')
    save([  OutputFolder '\' Namesession '_Phenosys.mat'],'t_ds','X_ds','X_ds_n','Y_ds','Y_r_ds', 'Whl_u','starts','stops','sessStartsIdx','SessStartsTimes','WhlSpeed','XSpeed','transitionThresh','CondNam')
    Phenodone=1
end

function [Traj, newW]= runTraj_vinca (t, X, transition, W, starts, stops, sessStartsIdx, Speed, XSpeed,transitionThresh,AnlPath, OutpuFolder, sesNam,CondNam,Circu)
Whl=W;
condition=1;
sessStartsIdx;
if Circu == 43
    edgesX=[0:2:200]; %NEW101bins 
elseif Circu==2*pi*15
    edgesX=[0:1.45:145]; %NEW101bins 
end
     highYlevel=mean(transition) ; %good
jj=1;
sf=1000;
Traj=struct ('transition',[],...
    'Cond', [],...
    'time',[],...
    'wheel',[],...
    'VRtraj',[],...
    'condition', [],...
    'Speed',[],...
    'XSpeed',[],...
    'binSpX',[],...
    'binSpW',[],...
    'WB',[],...
    'start',[],...
    'stop',[],...
    'tstart',[],...
    'tstop',[],...
    'endVR',[],...
    'meanSpeed',[],...%
    'maxSpeed',[],...%
    'meanSpeedVR',[],...%
    'maxSpeedVR',[],...%
    'FirstRunTime',[],...
    'pcentStopped',[],...
    'NumStops',[],...
    'distW',[], ...
    'distVR',[],...
    'dur',[],...
    'durVR',[],...
    'pcentStoppedVR',[],...
    'NumStopsVR',[],...
    'StopDurations',[],...
    'PauseTimes',[],...%
    'StopDurationsVR',[],...
    'PauseTimesVR',[],...%
    'PauseLocsVR',[],...
    'index', []);

NotATrial=[];

 for lap=1:length(starts)
%     if lap<length(starts) % the last start should not be a start actually...
        Traj(lap).tstart=t(starts(lap)+2);
        Traj(lap).start=(starts(lap)+2);
        Traj(lap).stop=(stops(lap)-2);
        Traj(lap).tstop=t(stops(lap)-2);
        Traj(lap).dur=Traj(lap).tstop-Traj(lap).tstart;
        figure(11)
        hold on
        plot(t(Traj(lap).start:Traj(lap).stop),X(Traj(lap).start:Traj(lap).stop))
     
        indexlap=find(t> Traj(lap).tstart & t< Traj(lap).tstop);
        temptraj=X(indexlap);
        [mx,imax]=max(temptraj);
        imax=find(temptraj>= mx-(0.01*mx),1,'first'); %when animal reaches end of track (1% under max)
        indexVR=indexlap(1):indexlap(imax); %think about adding the time after teleport
        
        Traj(lap).distVR=X(indexVR(end))-X(indexVR(1));
        Traj(lap).durVR=t(indexVR(end))-t(indexVR(1));
        Traj(lap).index=indexlap;%I will use this for jmaps
        Traj(lap).distW=W(indexlap(end))-W(indexlap(1));
        
        W(indexlap(1):end)=W(indexlap:end)-W(indexlap(1));
        W(indexlap(1):end)=W(indexlap(1):end)+abs(min(W(indexlap(1):end)));
        
        Traj(lap).time=t(indexlap)-Traj(lap).tstart;
        Traj(lap).trajTimes=t(indexlap);
        Traj(lap).transition=transition(indexlap);
        Traj(lap).wheel=Whl(indexlap)-Whl(indexlap(1)); %start from 0 at each start
        Traj(lap).VRtraj=X(indexlap);
        Traj(lap).Speed=Speed(indexlap);
        Traj(lap).XSpeed=XSpeed(indexlap);
        
        %compute speed in every VR bin
        [N,edges,bin] = histcounts(Traj(lap).VRtraj,edgesX);
        speed4bin=zeros(size(edgesX));
        for bb=1:length(edges)
            iid=find(bin==( bb));
            if ~isempty(iid)
                speed4bin((bb))=mean(Traj(lap).XSpeed(iid));
            else
                speed4bin((bb))=nan; % when I don't detect the animal inside
            end
        end
        figure(3)   
        plot(speed4bin)
        Traj(lap).binSpX=speed4bin; %REMPLACE PAR CETTE LIGNE SEULE
%         if Circu == 43  %JAI ENLEVE CA
%             speed4bin(edges<=10)=nan;% attention les valeurs sont un peu arbitraires
%             speed4bin(edges>max(edges)-10)=nan;
%             Traj(lap).binSpX=speed4bin;
%         elseif Circu == 2*pi*15
%             speed4bin(edges<=7)=nan;% attention les valeurs sont un peu arbitraires
%             speed4bin(edges>max(edges)-7)=nan;
%             Traj(lap).binSpX=speed4bin;
%         end
        
        
        Traj(lap).FirstRunTime=Traj(lap).trajTimes(find(Traj(lap).Speed>=0.2,1,'first'))-Traj(lap).trajTimes(1);
        runtim=find(Traj(lap).Speed>=0.2,1,'first');
        
        runperiod=Traj(lap).wheel(runtim:end); % this is for the whole time
        runperiodspeed=Traj(lap).Speed(runtim:end);
        tempT=Traj(lap).time(runtim:end);
        stopped=find(runperiodspeed<=0.5);
        Traj(lap).pcentStopped=(length(stopped)./length(runperiod))*100;
        Dtp=diff(stopped);
        [d, id] = getchunks(Dtp);
        %getchunksreturns an array of n elements, where n is the number
        %   of consecutive chunks (2 or more repetitions) in Dtp, and each element is
        %   the number of repetitions in each chunk.
        % I am still at 1KHz so I want at least 100ms i.e. 100 consecutive
        % elements
        
        
        if ~isempty (d)
            longones=find(d>500);
            if ~isempty(longones)
                Traj(lap).NumStops=length(d(longones));
                Traj(lap).StopDurations=d(longones)./sf;
                Pauses=stopped(id(longones));
                Traj(lap).PauseTimes=tempT(Pauses);
                
            else
                Traj(lap).NumStops=0;
                Traj(lap).StopDurations=[];
                Traj(lap).PauseTimes=[];
            end
        else
            Traj(lap).NumStops=0;
            Traj(lap).StopDurations=[];
            Traj(lap).PauseTimes=[];
        end
        
        Traj(lap).meanSpeed=nanmean(runperiodspeed);
        Traj(lap).maxSpeed=nanmax(runperiodspeed);
        
        % get stops in the VR period
        runperiodVR=X(indexVR);
        runperiodspeedVR=Speed(indexVR);
        tempT=t(indexVR);
        Traj(lap).endVR=t(indexVR(end));
        stoppedVR=find(runperiodspeedVR<=0.5);
        Traj(lap).pcentStoppedVR=(length(stoppedVR)./length(runperiodVR))*100;
        DtpVR=diff(stoppedVR);
        [d, id] = getchunks(DtpVR);
        
        if ~isempty(d)
            longones=find(d>500);
            if ~isempty(longones)
                Traj(lap).NumStopsVR=length(d(longones));
                Traj(lap).StopDurationsVR=d(longones)./sf;
                Pauses=stoppedVR(id(longones));
                Traj(lap).PauseTimesVR=tempT(Pauses);
                Traj(lap).PauseLocsVR=runperiodVR(Pauses);
            else
                Traj(lap).NumStopsVR=0;
                Traj(lap).StopDurationsVR=[];
                Traj(lap).PauseTimesVR=[];
                Traj(lap).PauseLocsVR=[];
            end
        else
            Traj(lap).NumStopsVR=0;
            Traj(lap).StopDurationsVR=[];
            Traj(lap).PauseTimesVR=[];
            Traj(lap).PauseLocsVR=[];
        end
        Traj(lap).meanSpeedVR=nanmean(runperiodspeed);
        Traj(lap).maxSpeedVR=nanmax(runperiodspeed);

        
        figure(5)
        plot(Traj(lap).transition)
        a=mean(Traj(lap).transition);

        if(mean(Traj(lap).transition)> highYlevel && Traj(lap).distVR>50)
            Traj(lap).WB='B';
                 
                    
        elseif (mean(Traj(lap).transition)<transitionThresh && Traj(lap).distVR>50)
            Traj(lap).WB='W';

        else
            Traj(lap).WB='Reset';
            NotATrial=[NotATrial; lap];
        end
      
        if starts(lap) <  sessStartsIdx(condition); %si le début du lap est avant le début de la condition = condition 
            condition = condition-1 ;

        elseif starts(lap) > sessStartsIdx(condition) 
            if numel(sessStartsIdx) >= condition + 1
                 if   starts(lap) > sessStartsIdx(condition+1) 
                    condition=condition+1;
                 elseif abs(starts(lap)-sessStartsIdx(condition+1)) <12000
                        condition=condition+1;
                    
                 else
                     condition=condition;
                 end
            end
        end
        Traj(lap).Cond=condition;
        Traj(lap).condition=CondNam{Traj(lap).Cond};
% 
%         if strcmp(Traj(lap).condition,'POM')
%             if mod(jj,2)
%                 Traj(lap).WB='B';
%             else
%                 Traj(lap).WB='W';
%             end
%             jj=jj+1;
%         elseif strcmp(Traj(lap).condition,'POnM')
%             if lap==1;
%                 Traj(lap).WB='W';
%             elseif lap==2;
%                 Traj(lap).WB='B';
%             elseif mod(jj,2)
%                 Traj(lap).WB='B';
%             else
%                 Traj(lap).WB='W';
%             end
%             jj=jj+1;
%                elseif strcmp(Traj(lap).condition,'POnM9')
%             if lap==1;
%                 Traj(lap).WB='W';
%             elseif lap==2;
%                 Traj(lap).WB='B';
%             elseif mod(jj,2)
%                 Traj(lap).WB='B';
%             else
%                 Traj(lap).WB='W';
%             end
%             jj=jj+1;
%         elseif  strcmp(Traj(lap).condition,'POMb')
%             Traj(lap).WB='B';
%         end
% 
%         

%     end
end
Traj(NotATrial)=[];
 for  lap=1:length(Traj)
     if any(Traj(lap).WB=='B') 
         Traj(lap).icondway_tr=Traj(lap).Cond*2;
     else 
         
         Traj(lap).icondway_tr=Traj(lap).Cond*2-1;
     
     end
 end



%% ici A REMETTRE si il y a le X wheel!!!!

 newW=W;

%figure(8)
%plot(newW)
%
edgesW=[0:5:max(newW)];
for lap= 1:length(Traj)
    %compute speed in each bin of ball distance
    figure(25)
    hold on
    plot(Traj(lap).wheel)
    [N,edges,bin] = histcounts(Traj(lap).wheel,edgesW);
    speed4bin=zeros(size(edgesW));
    for bb=1:length(edges)
        iid=find(bin== (bb));
        if ~isempty(iid)
            speed4bin(bb)=mean(Traj(lap).Speed(iid));
        else
            speed4bin(bb)=nan; % when I don't detect the animal inside
        end
    end
    Traj(lap).binSpW=speed4bin;
    figure(4)
    hold on
    plot(speed4bin)
    
end




%     end


%% make a table from the structure and save it as a csv (excel readable file)
fields2remove={'start','stop','transition','time', 'wheel', 'VRtraj','Speed','XSpeed', 'index','trajTimes','binSpX','binSpW'};
OutStruct=rmfield(Traj,fields2remove);
outTab=struct2table(OutStruct);
cd (OutpuFolder)
writetable(outTab,[sesNam,'.xls'])% output to an excell file



cd (AnlPath)
save([AnlPath '\process\' 'TrajData.mat'],'Traj','newW');% ask user input to re-write?
save([OutpuFolder '\' sesNam '_TrajData.mat'],'Traj','newW')

end

function [Spikes,phyData,goodClusListIdx,datacell,headercell,outputt]=getphyInfo(CircusPath)
% gets all info from phy

[Spikes,phyData,goodClusListIdx]=Load_phyResults(CircusPath);
cd (CircusPath)
[datacell, headercell,outputt]=tsvread2('cluster_info.tsv');
end

function [Spikes,phyData,goodClusListIdx]=Load_phyResults(phyOutputFolder)
curDir=cd; %keep current directory in memory

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


%go back to original directory
cd(curDir);
clusgroup= phyData.cluster_ids(:,2);
%find good ones
goodClusListIdx=cellfun(@(clusgroup) strcmp(clusgroup,'good'), phyData.cluster_ids(:,2));
goodClus=phyData.cluster_ids(goodClusListIdx,1);
goodClusIdx=ismember(phyData.spike_clusters,str2double(goodClus));
goodTemplates=unique(phyData.spike_clusters(goodClusIdx)+1); % FIX ME for some reason, need to add +1 to template number

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


end

function varargout = tsvread2( varargin )
%[data, header, raw] = tsvread( file ) reads in text file with tab-seperated variables. default value for data is nan.
%alternative input/output option is suppluying header strings
%[col1, col2, col3, ..., header, raw] = tsvread( file, header1, header2, header3, ... )
%header is the first row (assumed to have header names) and raw is the imported text
%if a vector is supplied, this specifies the number rows to be imported.
%examples:
%[col1, col2,col3,header,raw] = tsvread( 'example.tsv', 'header1', 'header2', 'header3', 1:5 )
%will import data from example.tsv, and cols corresponding to header1,
%header2, header3, cols 1 t0 5.
%if no outputs are requested, then a portion of the rquested table is
%displayed -- good idea to see how the import and header requests are
%working!
%If there is a tsv file in local directory, then just running "tsvread" at
%the command line will read the newest tsv file and display first ten rows to screen.
%sak 9/2/11
if nargin == 0
    disp( 'importing most recent tsv file in local directory' );
    d = dir( '*.tsv' );
    [~,i] = sort( [d.datenum], 'descend' );
    fprintf( 'tsvread( ''%s'', 1:10 )\n', d(i(1)).name );
    eval( sprintf( 'tsvread( ''%s'', 1:10 )', d(i(1)).name ) );
    return;
end;
fid = fopen( varargin{1}, 'r' );
if nargout == 0
    fprintf( 'loading %s, and displaying sideways\n', varargin{1} );
end
varargin(1) = [];
% stuff = textscan( fid, '%s', 'delimiter', '\n');  %OLD vinca change le 08/08/24 car délimiter faux pour moi 
 stuff = textscan( fid, '%s', 'delimiter', '\t ');
stuff = stuff{1,:};
fclose(fid);
j=1;




stuffs=[];
for i = 1:9:length(stuff);
    stuffs=[stuffs;[stuff(i:i+8)].'];
end


outputt=[];
outputt=[stuffs(:,1)];
outputt(:,2)=[stuffs(:,2)];

numrows = 1:size(stuff,1);
ind = cellfun( 'isclass', varargin, 'double' );
if any( ind )
    numrows = varargin{ind};
    varargin(ind) = [];
    if numel(numrows) == 1
        numrows = 1:min(numrows, size( stuff, 1) );
    end
end
numrows = intersect( numrows, 1:size( stuff, 1 ) );
header = regexprep( regexp( stuff{1}, '[^\t]*\t', 'match' ), '\t', '' );
raw = repmat( {}, numel(numrows), numel(header) );
for i=numrows
    stuff{i}(end+1) = 9;
    tmp = regexprep( regexp( stuff{i}, '[^\t]*\t', 'match' ), '\t', '');
    raw(i,1:numel(tmp)) = tmp;
end
header = raw(1,:);
data = nan(size(raw));
for i=numrows
    for j=1:size( raw, 2 )
        if ~isempty( raw{i,j} )
            [a, count, errmsg] = sscanf( raw{i,j}, '%f' );
            if ~isempty( a )
                data(i,j) = a;
            end
        end
    end
end
%%
    varargout = {data, header, raw,outputt}; %added stuffs (vinca 08/08/2024)
    
end

function [allcel]=make_ephy(cluster_infos_tsv,datacell,phyData,Spikes,goodClusListIdx,Namesession,AnlPath,whatprob,isoe,OutputFolder)
cd(AnlPath)

if isoe
          A=dir;
          idx = find(startsWith({A.name}, '10'), 1);
         if isempty(idx), idx = 4; end
         
% OPENEPHYS FORMAT
        if strcmp( A(idx).name(4:5),'_1') 
            prefix =  A(idx).name(1:4);
            isorder=0;
            channellist=[];
                for i=1:64
                Nam=strcat(prefix,string(i),'.continuous');
                file_path=fullfile(AnlPath, Nam);
                fid = fopen(file_path); 
                [~, file_name] = fileparts(file_path);
                string_nb = 1024;
                hdr = fread(fid, string_nb, 'char*1');
                eval(char(hdr')); 
                chh=header.channel;
                fclose(fid);
                num_str = regexp(chh, '\d+', 'match', 'once');  % extrait '123'
                num = str2double(num_str);
                channellist=[channellist;num];
                end
        elseif  strcmp(A(idx).name(4:7),'_CH1') 
            prefix =  A(idx).name(1:6);
            isorder=0;
            
        elseif strcmp(A(idx).name(4:8),'_ADC1')
            prefix = A(idx).name(1:4);
            isorder=1;
         end
end
 shkcol=find(strcmp(cluster_infos_tsv,'sh'));%looking in the header what column corresponds to shank
idcol=find(strcmp(cluster_infos_tsv,'cluster_id'));% and clusterID
% datacell(1,:)=[]; %%%datacell is the table with all the units caracteristics
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
shanksize=probechannels(whatprob,1);
allcel.waveform=nan(length([-20:30]),length(shanksize),100,length(allcel.id_cel));

for i=1:length(unique(Spikes.Units))
    disp(strcat('cell n°',string(i),'/',string(length(allcel.id_cel))))
    channeldata=[];
    idx =bestchannelsunitsgood(i);
    channel=bestchannelsgood(i);% logical indexing 
    mask=(double(Spikes.Units+1)==idx); %NEW -1 pour des questions de base zéro
    selectedtimespk= allcel.itime_spk(mask,1);
    channelshank=probechannels(whatprob,channel);
    timespktemp=allcel.time_spk(mask);
    NN=randperm(length(selectedtimespk),100);
    MM=sort(selectedtimespk(NN));
 %openephys format      
 if isoe
    for k=1:length(channelshank);    
   
            if isorder
             Nam=[prefix,'CH',num2str(channelshank(k)),'.continuous']; 
             
             max_records_to_read=ceil(MM(end) * 25000 / 1024);
             [dischanneldata, ~, ~] = fct_read_continuous_filepart(fullfile(AnlPath, Nam),max_records_to_read); 
%               dischanneldata = [dischanneldata].';
             channeldata=[channeldata;dischanneldata];   
            else
             chh=channelshank(k);
             chh=find(chh==channellist);
             Nam=[prefix,num2str(chh),'.continuous'];
             max_records_to_read=ceil(MM(end) * 25000 / 1024);
             [dischanneldata, ~, ~] = fct_read_continuous_filepart(fullfile(AnlPath, Nam),max_records_to_read); 
%              dischanneldata = [dischanneldata].';
             channeldata=[channeldata;dischanneldata];   
            end
    end
 else
        if isfile('continuous.dat')
            [channeldata, ~]=Read_OEP_Binary ('continuous.dat',channelshank, 0,MM(end)/25000, 1); 
        elseif ~isfile('continuous.dat')  %equivalent à si continuous.dat est ailleurs (pas sorti de rythmFGPA) ; ne marchera surement pas car val n'exisste pas si on a déja fait traj etc. j'ai une soluce en tête si besoin :) 
            path_dat = strcat(AnlPath, filesep, 'continuous\Rhythm_FPGA-', num2str(val.continuous.source_processor_id), '.', num2str(val.continuous.source_processor_sub_idx), filesep); %va chercher la localisation du .dat dans les sous-dossiers
            [dischanneldata, ~]=Read_OEP_Binary_altered(strcat(path_dat,'continuous.dat'),channelshank, 0,MM(end)/25000+10, sr); %la version "altérée consiste en la modification d'une ligne (D.timestamps) où j'ai ajouté le path_dat (comme au dessus) avant continuous.dat dans les inputs de ReadNPY pour qu'il trouve le fichier quin'est pas au même endroit que oebin ! c'est tout!
             channeldata=[channeldata;dischanneldata.'];   
             
   
    end
 end
    selectedamptimespk=selectedtimespk(1:100);
    for j=1:100 %j = nb de spike sélectionnés pour le sampling
        timewindow=[-20:30]; %fenetre pour le spike 
        timewindow=timewindow+double(selectedamptimespk(j));
           if timewindow(1) <= 0
               timewindow=timewindow+abs(timewindow(1))+1;
           end
%                   plot(channeldata(find(channel==channelshank),timewindow).')
            allcel.waveform(:,:,j,i)=channeldata(:,timewindow).';
        
    end
    allcel.bestchannelid(:,i)=(find(channel==channelshank));
end
allcel.bestchannel=bestchannelsgood;
allcel.meanwaveform=squeeze(mean(allcel.waveform,3));
for i=1:length(unique(Spikes.Units))
     allcel.bestwaveform(:,i)=allcel.meanwaveform(:,allcel.bestchannelid(i),i);
     plot(allcel.bestwaveform(:,i))
      allcel.bestswaveforms(:,:,i)=allcel.waveform(:,allcel.bestchannelid(i),:,i);
end
save([AnlPath  '\process\ePhy.mat'],'allcel')
save([OutputFolder '\' Namesession '_ePhy.mat'],'allcel')
end

function [data, time, header, record_nb] = fct_read_continuous_filepart(file_path, max_records_to_read)

fid = fopen(file_path);
[~, file_name] = fileparts(file_path);

% Tailles en octets
file_size = fct_getfilesize(fid);
header_size = 1024*1;
timestamps_size = 1*8;
sample_nb_size = 1*2;
record_ind_size = 1*2;
samples_size = 1024*2;
record_marker_size = 10*1;
record_size = timestamps_size + sample_nb_size + record_ind_size + samples_size + record_marker_size;

string_nb = 1024;

hdr = fread(fid, string_nb, 'char*1');
eval(char(hdr')); 
clear hdr

% Si pas de limite fournie, on lit tout
if nargin < 2
    max_records_to_read = inf;
end

i = 1;
while ftell(fid) + record_size <= file_size && i <= max_records_to_read
    timestamp(i) = fread(fid, 1, 'int64', 0, 'l');
    samples_nb(i) = fread(fid, 1, 'uint16',0,'l');

    if i > 1
        samples_nb_sum(i) =  samples_nb(i) + samples_nb_sum(i-1);
    else
        samples_nb_sum(i) = samples_nb(i);
    end

    record_ind(i) = fread(fid, 1, 'uint16');
    data((samples_nb_sum(i) - samples_nb(i) + 1):samples_nb_sum(i)) = fread(fid, samples_nb(i), 'int16', 0, 'b');
    time((samples_nb_sum(i) - samples_nb(i) + 1):samples_nb_sum(i)) = timestamp(i):timestamp(i)+samples_nb(i)-1;

    fread(fid, 10, 'char*1'); 
    i = i + 1;
end

fclose(fid);

record_nb = length(timestamp);

% Conversion à partir du Header
time = time / header.sampleRate; % en secondes
data = data * header.bitVolts; % en microVolts
end

function [channellist]=probechannels(probetype,channel);
%Select all the channels that are oon 1 shank of a seen spike
if probetype == 5
    if channel <= 12
        channellist=[1,2,3,4,5,6,7,8,9,10,11,12];
        
    elseif channel  <= 24 & channel >  12;
           channellist= [13,14,15,16,17,18,19,20,21,22,23,24];
     elseif channel <= 40 & channel >  24;
        channellist=[25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40];
    elseif channel <= 52 & channel >  40;
        channellist=[41,42,43,44,45,46,47,48,49,50,51,52];
    elseif  channel >  52;
        channellist=[53,54,55,56,57,58,59,60,61,62,63,64];
    end
elseif probetype == 64
    if channel <= 8
        channellist = [1:8];
    elseif channel  <= 16 & channel >  8;
        channellist = [9:16];
            elseif channel  <= 24 & channel > 16;
        channellist = [17:24];
            elseif channel  <= 32 & channel >  24;
        channellist = [25:32];
        elseif channel  <= 40 & channel >  32;
        channellist = [33:40];
            elseif channel  <= 48 & channel >  40;
        channellist = [41:48];
            elseif channel  <= 56 & channel >  48;
        channellist = [49:56];
            elseif channel  <= 64 & channel >  56;
        channellist = [57:64];
    end
elseif probetype == 77 
    if  ismember(channel,[16,21,18,23,20,24,27,30,19,29,17,25,26,32,28,22])
        channellist = [16,21,18,23,20,24,27,30,19,29,17,25,26,32,28,22];
    elseif ismember(channel,[11,1,13,3,15,5,31,7,14,9,12,10,8,6,4,2])
        channellist=[11,1,13,3,15,5,31,7,14,9,12,10,8,6,4,2];
    elseif ismember(channel,[54,64,52,62,50,60,34,58,51,56,53,55,57,59,61,63])
        channellist=[54,64,52,62,50,60,34,58,51,56,53,55,57,59,61,63];
    elseif ismember(channel,[49,44,47,42,45,41,38,35,46,36,48,40,39,33,37,43])
        channellist = [49,44,47,42,45,41,38,35,46,36,48,40,39,33,37,43];
    end
   
end
end

%% === Fonctions complémentaires === 
function confirmed = askConfirmation(prompt)
    while true
        response = input([prompt, ' (Y/N): '], 's');
        if strcmpi(response, 'Y')
            confirmed = true;
            return;
        elseif strcmpi(response, 'N')
            confirmed = false;
            return;
        else
            disp('Réponse invalide. Tapez Y ou N.');
        end
    end
end