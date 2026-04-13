function [Traj, newW]= runTraj (t, X, transition, W, starts, sessS, Speed, XSpeed,transitionThresh,AnlPath, OutpuFolder, sesNam,CondNam)
Whl=W;
%edgesX=[0:2.5:max(X)];
%edgesX=[0:2.5:145]; OLD59bin %1,8125 pour avoir 80 bins  
edgesX=[0:1.45:145]; %NEW101bins 
jj=1;
sf=1000;
Traj=struct ('transition',[],...
    'Cond', [],...
    'time',[],...
    'wheel',[],...
    'VRtraj',[],...
    'Speed',[],...
    'XSpeed',[],...
    'binSpX',[],...
    'binSpW',[],...
    'LR',[],...
    'start',[],...
    'stop',[],...
    'tstart',[],...
    'tstop',[],...
    'endVR',[],...
    'mnSpeed',[],...%
    'maxSpeed',[],...%
    'mnSpeedVR',[],...%
    'maxSpeedVR',[],...%
    'FirstRunTime',[],...
    'pcentStopped',[],...
    'NumStops',[],...
    'dist',[], ...
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
    'index', [],...
    'condition', []);

NotATrial=[];

if transitionThresh==2600 %%%%original
%if transitionThresh==2300   
    highYlevel=2900;
    highYlevel=2700; %for session MM1_2024-03-20_15-02-54 %METTRE 1750 pour sess ou manque retours par ex
else
    highYlevel=2.25 ; %good
%     highYlevel=-0.6 ;%vérifier
end
for ss=1:length(starts)
    if ss<length(starts) % the last start should not be a start actually...
        Traj(ss).tstart=t(starts(ss)+2);
        Traj(ss).start=(starts(ss)+2);
        Traj(ss).stop=(starts(ss+1)-2);
        Traj(ss).tstop=t(starts(ss+1)-2);
        Traj(ss).dur=Traj(ss).tstop-Traj(ss).tstart;
%     elseif ss==length(starts)
%         last=X(Traj(ss-1).start:end)
%         Traj(ss).stop=Traj(ss)find(abs(diff(last))>120);
%         Traj(ss).tstart=t(starts(ss)+2);
%         Traj(ss).start=(starts(ss)+2);
%         Traj(ss).stop=(starts(ss+1)-2);
%         Traj(ss).tstop=t(starts(ss+1)-2);
%         Traj(ss).dur=Traj(ss).tstop-Traj(ss).tstart;
%         
        index=find(t> Traj(ss).tstart & t< Traj(ss).tstop);
        tmptraj=X(index);
        [mx,ii]=max(tmptraj);
        ii=find(tmptraj>= mx-(0.01*mx),1,'first'); %xhen animal reaches end of track (1% under max)
        indexVR=index(1):index(ii); %think about adding the time after teleport
        
        Traj(ss).index=index;%I will use this for jmaps
        Traj(ss).dist=W(index(end))-W(index(1));
        %         Traj(ss).distfromwall= %this is ridiculous I know in VR it's a
        %         fixed distance
        W(index(1):end)=W(index:end)-W(index(1));
        W(index(1):end)=W(index(1):end)+abs(min(W(index(1):end)));
        Traj(ss).time=t(index)-Traj(ss).tstart;
        Traj(ss).trajT=t(index);
        Traj(ss).transition=transition(index);
        Traj(ss).wheel=Whl(index)-Whl(index(1)); %start from 0 at each start
        Traj(ss).VRtraj=X(index);
        Traj(ss).Speed=Speed(index);
        Traj(ss).XSpeed=XSpeed(index);
        
        %compute speed in every VR bin
        [N,edges,bin] = histcounts(Traj(ss).VRtraj,edgesX);
        speed4bin=zeros(size(edgesX));
        for bb=1:length(edges)
            iid=find(bin==( bb));
            if ~isempty(iid)
                speed4bin((bb))=mean(Traj(ss).XSpeed(iid));
            else
                speed4bin((bb))=nan; % when I don't detect the animal inside
            end
        end
        
        
        speed4bin(edges<=6)=nan;% attention les valeurs sont un peu arbitraires
        speed4bin(edges>max(edges)-6)=nan;
        Traj(ss).binSpX=speed4bin;
        
        
        
        Traj(ss).FirstRunTime=Traj(ss).trajT(find(Traj(ss).Speed>=0.2,1,'first'))-Traj(ss).trajT(1);
        runtim=find(Traj(ss).Speed>=0.2,1,'first');
        runperiod=Traj(ss).wheel(runtim:end); % this is for the whole time
        runperiodspeed=Traj(ss).Speed(runtim:end);
        tempT=Traj(ss).time(runtim:end);
        stopped=find(runperiodspeed<=0.5);
        Traj(ss).pcentStopped=(length(stopped)./length(runperiod))*100;
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
                Traj(ss).NumStops=length(d(longones));
                Traj(ss).StopDurations=d(longones)./sf;
                Pauses=stopped(id(longones));
                Traj(ss).PauseTimes=tempT(Pauses);
                
            else
                Traj(ss).NumStops=0;
                Traj(ss).StopDurations=[];
                Traj(ss).PauseTimes=[];
            end
        else
            Traj(ss).NumStops=0;
            Traj(ss).StopDurations=[];
            Traj(ss).PauseTimes=[];
        end
        
        Traj(ss).mnSpeed=nanmean(runperiodspeed);
        Traj(ss).maxSpeed=nanmax(runperiodspeed);
        
        % get stops in the VR period
        runperiodVR=X(indexVR);
        runperiodspeedVR=Speed(indexVR);
        tempT=t(indexVR);
        Traj(ss).endVR=t(indexVR(end));
        stoppedVR=find(runperiodspeedVR<=0.5);
        Traj(ss).pcentStoppedVR=(length(stoppedVR)./length(runperiodVR))*100;
        DtpVR=diff(stoppedVR);
        [d, id] = getchunks(DtpVR);
        
        if ~isempty(d)
            longones=find(d>500);
            if ~isempty(longones)
                Traj(ss).NumStopsVR=length(d(longones));
                Traj(ss).StopDurationsVR=d(longones)./sf;
                Pauses=stoppedVR(id(longones));
                Traj(ss).PauseTimesVR=tempT(Pauses);
                Traj(ss).PauseLocsVR=runperiodVR(Pauses);
            else
                Traj(ss).NumStopsVR=0;
                Traj(ss).StopDurationsVR=[];
                Traj(ss).PauseTimesVR=[];
                Traj(ss).PauseLocsVR=[];
            end
        else
            Traj(ss).NumStopsVR=0;
            Traj(ss).StopDurationsVR=[];
            Traj(ss).PauseTimesVR=[];
            Traj(ss).PauseLocsVR=[];
        end
        Traj(ss).mnSpeedVR=nanmean(runperiodspeed);
        Traj(ss).maxSpeedVR=nanmax(runperiodspeed);
        Traj(ss).distVR=X(indexVR(end))-X(indexVR(1));
        Traj(ss).durVR=t(indexVR(end))-t(indexVR(1));
        
        figure(5)
        plot(Traj(ss).transition)
        a=mean(Traj(ss).transition);

        if(mean(Traj(ss).transition)> highYlevel && Traj(ss).distVR>50)
            Traj(ss).LR='B';
                 
                    
        elseif (mean(Traj(ss).transition)<transitionThresh && Traj(ss).distVR>50)
            Traj(ss).LR='W';

        else
            Traj(ss).LR='Reset';
            NotATrial=[NotATrial; ss];
        end
        indSess=find(sessS> Traj(ss).start,1,'first');
        if indSess==1
                        Traj(ss).Cond=indSess;
                        Traj(ss).condition=CondNam{Traj(ss).Cond};
        elseif~isempty(indSess)

%             Traj(ss).Cond=(Traj(ss-1).Cond)+1
                Traj(ss).Cond=indSess;
                                        Traj(ss).condition=CondNam{Traj(ss).Cond};
        else
               Traj(ss).Cond=length(sessS)+1; %old bug 
%                 Traj(ss).Cond=(Traj(ss-1).Cond)+1 %still buggy
                        Traj(ss).condition=CondNam{Traj(ss).Cond}
        end
        if strcmp(Traj(ss).condition,'POM')
            if mod(jj,2)
                Traj(ss).LR='B';
            else
                Traj(ss).LR='W';
            end
            jj=jj+1;
        elseif strcmp(Traj(ss).condition,'POnM')
            if ss==1;
                Traj(ss).LR='W';
            elseif ss==2;
                Traj(ss).LR='B';
            elseif mod(jj,2)
                Traj(ss).LR='B';
            else
                Traj(ss).LR='W';
            end
            jj=jj+1;
               elseif strcmp(Traj(ss).condition,'POnM9')
            if ss==1;
                Traj(ss).LR='W';
            elseif ss==2;
                Traj(ss).LR='B';
            elseif mod(jj,2)
                Traj(ss).LR='B';
            else
                Traj(ss).LR='W';
            end
            jj=jj+1;
        elseif  strcmp(Traj(ss).condition,'POMb')
            Traj(ss).LR='B';
        end

        

    end
end
Traj(NotATrial)=[];
for  ss=1:length(Traj)
        if Traj(ss).Cond==1
        Traj(ss).icondway_tr=Traj(ss).Cond+(any(strcmp(Traj(ss).LR,'B'))*2)-(1*any(strcmp(Traj(ss).LR,'B')));
        elseif any(Traj(ss).LR=='B')
             Traj(ss).icondway_tr=Traj(ss).Cond+(any(strcmp(Traj(ss).LR,'B'))*2);
        else 
            Traj(ss).icondway_tr=Traj(ss).Cond+(any(strcmp(Traj(ss).LR,'W'))*2)-1;
        end
end





% plot(t,X)
% hold on
% for tt=1:length(Traj)
% plot(Traj(tt).time+t(Traj(tt).start), Traj(tt).transition,'k')
% end

%% ici A REMETTRE si il y a le X wheel!!!!

 newW=W;

%figure(8)
%plot(newW)
%
edgesW=[0:5:max(newW)];
for ss= 1:length(Traj)
    %compute speed in each bin of ball distance
    figure(5)
    plot(Traj(ss).wheel)
    [N,edges,bin] = histcounts(Traj(ss).wheel,edgesW);
    speed4bin=zeros(size(edgesW));
    for bb=1:length(edges)
        iid=find(bin== (bb));
        if ~isempty(iid)
            speed4bin(bb)=mean(Traj(ss).Speed(iid));
        else
            speed4bin(bb)=nan; % when I don't detect the animal inside
        end
    end
    Traj(ss).binSpW=speed4bin;
    
    
end




%     end


%% make a table from the structure and save it as a csv (excel readable file)
fields2remove={'start','stop','transition','time', 'wheel', 'VRtraj','Speed','XSpeed', 'index','trajT','binSpX','binSpW'};
OutStruct=rmfield(Traj,fields2remove);
outTab=struct2table(OutStruct);
cd (OutpuFolder)
writetable(outTab,[sesNam,'.xls'])% output to an excell file



cd (AnlPath)
save([AnlPath '\process\' 'TrajData.mat'],'Traj','newW');% ask user input to re-write?
save([OutpuFolder '\' sesNam '_TrajData.mat'],'Traj','newW')


