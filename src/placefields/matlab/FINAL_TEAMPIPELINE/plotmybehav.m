function []=plotmybehav (Traj,edgesX,edgesW, OutpuFolder,sesNam,CondNam)

nbsess= length(unique([Traj.Cond]));
plotind=1;

f0=figure('Name',['Behavior'],'NumberTitle','off');
for cc=unique([Traj.Cond])
    fun = @(x) Traj(x).Cond == cc; % useful for complicated fields
    tf2 = arrayfun(fun, 1:numel(Traj));
    index2 = find(tf2);
    
    
    for tt=1:length(index2)
        if strcmp(Traj(index2(tt)).LR,'RL')
            subplot(5,4,plotind)
            plot(Traj(index2(tt)).time, Traj(index2(tt)).VRtraj,'r')
            hold on
            title(CondNam(cc))
            xlabel('Time (s)')
            ylabel('VR Location (pixels)')
            
            subplot(5,4,plotind+1)
            plot(Traj(index2(tt)).time, Traj(index2(tt)).wheel,'r')
            hold on
            title('Wheel')
            xlabel('Time (s)')
            ylabel('Distance (cm)')
            
            subplot(5,4,plotind+2)
            plot(edgesX, Traj(index2(tt)).binSpX,'r')% utiliser la matrice utilisée pur la moyenne
            hold on
            xlabel('VR Location (pix)')
            ylabel('Speed (cm/s)')
            %             xlim([0 200])
            %             ylim([0 150])
            subplot(5,4,plotind+3)
            plot(edgesW, Traj(index2(tt)).binSpW,'r')
            hold on
            xlabel('Wheel Dist (pix)')
            ylabel('Speed (cm/s)')
            %             xlim([0 300])
            
        elseif strcmp(Traj(index2(tt)).LR,'LR')
            subplot(5,4,plotind)
            plot(Traj(index2(tt)).time, Traj(index2(tt)).VRtraj,'b')
            hold on
            title(CondNam(cc))
            xlabel('Time (s)')
            ylabel('VR Location (pixels)')
            
            subplot(5,4,plotind+1)
            plot(Traj(index2(tt)).time, Traj(index2(tt)).wheel,'b')
            hold on
            title('Wheel')
            xlabel('Time (s)')
            ylabel('Distance (cm)')
            
            subplot(5,4,plotind+2)
            plot(edgesX,Traj(index2(tt)).binSpX,'b')
            hold on
            xlabel('VR Location (pix)')
            ylabel('Speed (cm/s)')
            %              xlim([0 200])
            %              ylim([0 150])
            subplot(5,4,plotind+3)
            plot(edgesW, Traj(index2(tt)).binSpW,'b')
            hold on
            xlabel('Wheel Dist (pix)')
            ylabel('Speed (cm/s)')
            %             xlim([0 300])
        end
        
    end
    
    plotind=plotind+4;
    
end


cd (OutpuFolder)
OutputNam=[sesNam,'_Behavior','_Raw','.pdf'];
subplot(5,4,1)
legend({'RL', 'LR'},'Location','SouthEast')
f0.PaperUnits='Normalized';
f0.PaperPosition=[0 0 1 1];
%
saveas(f0,OutputNam,'pdf')
% close (f0);
f1=figure('Name',['BehaviorMeanVR'],'NumberTitle','off');

plotind=1;
for cond=1:nbsess
    subplot(5,2,plotind)
    
    index=find(([Traj.Cond]==cond) & strcmp({Traj.LR}, 'LR'));
    if~isempty(index)
        tempMatXLR=vertcat(Traj([index]).binSpX);
        
        %         errorbar(edgesX,mean(tempMatXLR),std(tempMatXLR),'r')
        MN=nanmean(tempMatXLR);
        SD=nanstd(tempMatXLR);
        Yd=[MN(4:end-4)+SD(4:end-4), fliplr(MN(4:end-4)-SD(4:end-4))];
        Xd=[edgesX(4:end-4), fliplr(edgesX(4:end-4))];
        fill(Xd, Yd,[0 0.5 0.8],'linestyle','none','FaceAlpha',0.3);
        
        hold on
        
        plot(edgesX(4:end-4),MN(4:end-4),'--b')% the nans prevent me to use fill
        ylabel('Speed (cm.s^-^1)')
%         LRobjects=[24, 48, 130];OLD
        LRobjects=[27, 46, 127];
        if strcmp(CondNam{cond}, ('NO-Obj'))
            plot(LRobjects,ones(size(LRobjects)),'Marker','^','MarkerEdgeColor',[0.8 0.8 0.8], 'MarkerSize',8,'linestyle','none')
        else
            plot(LRobjects,ones(size(LRobjects)),'^b','MarkerSize',8,'linestyle','none')
        end
        title([CondNam{cond}, ' -->'])
    end
    
    subplot(5,2,plotind+1)
    
    index2=find(([Traj.Cond]==cond) & strcmp({Traj.LR}, 'RL'));
    if~isempty(index2)
        tempMatXRL=vertcat(Traj([index2]).binSpX);
        MN=nanmean(tempMatXRL);
        SD=nanstd(tempMatXRL);
        Yd=[MN(4:end-4)+SD(4:end-4), fliplr(MN(4:end-4)-SD(4:end-4))];
        Xd=[edgesX(4:end-4), fliplr(edgesX(4:end-4))];
        fill(Xd, Yd,[0.9 0.6 0.6],'linestyle','none','FaceAlpha',0.3);
        
        hold on
        
        plot(edgesX(4:end-4),MN(4:end-4),'--r')% the nans prevent me to use fill
        
%         RLobjects=[15, 102, 120]; OLD
           RLobjects=[18, 99, 118]; 
        if strcmp(CondNam{cond}, ('NO-Obj'))
            plot(RLobjects,ones(size(RLobjects)),'Marker','^','MarkerEdgeColor',[0.8 0.8 0.8], 'MarkerSize',8,'linestyle','none')
        else
            plot(RLobjects,ones(size(RLobjects)),'Marker','^','MarkerEdgeColor','r', 'MarkerSize',8,'linestyle','none')
        end
        title([CondNam{cond}, ' <--'])
        set(gca,'xDir','reverse')
    end
    
    plotind=plotind+2;
    
end
text(0, -1*max(MN(~isnan(MN))), sesNam)
subplot(5,2,9)
xlabel('location in VR (cm)')
subplot(5,2,10)
xlabel('location in VR (cm)')


OutputNam1=[sesNam,'Behavior','_meanVR','.pdf'];

f2=figure('Name',['BehaviorMeanDistce'],'NumberTitle','off');
plotind=1;
for cond=1:nbsess
    subplot(5,2,plotind)
    
    index=find(([Traj.Cond]==cond) & strcmp({Traj.LR}, 'LR'));
    if~isempty(index)
        tempMatXLR=vertcat(Traj([index]).binSpW);
        MN=mean(tempMatXLR);
        SD=std(tempMatXLR);
        Yd=[MN(~isnan(MN))+SD(~isnan(MN)), fliplr(MN(~isnan(MN))-SD(~isnan(MN)))];
        Xd=[edgesW(~isnan(MN)), fliplr(edgesW(~isnan(MN)))];
        fill(Xd, Yd,[0 0.5 0.8],'linestyle','none','FaceAlpha',0.3);
        hold on
        plot(edgesW((~isnan(MN))),MN((~isnan(MN))),'--b')% the nans prevent me to use fill
        
        
        
        ylabel('Speed (cm.s^-^1)')
        title([CondNam{cond}, ' -->'])
    end
    
    subplot(5,2,plotind+1)
    
    index2=find(([Traj.Cond]==cond) & strcmp({Traj.LR}, 'RL'));
    if~isempty(index2)
        tempMatXRL=vertcat(Traj([index2]).binSpW);
        
        MN=mean(tempMatXRL);
        SD=std(tempMatXRL);
        Yd=[MN(~isnan(MN))+SD(~isnan(MN)), fliplr(MN(~isnan(MN))-SD(~isnan(MN)))];
        Xd=[edgesW(~isnan(MN)), fliplr(edgesW(~isnan(MN)))];
        fill(Xd, Yd,[0.9 0.6 0.6],'linestyle','none','FaceAlpha',0.3);
        
        hold on
        
        plot(edgesW(~isnan(MN)),MN(~isnan(MN)),'--r')% the nans prevent me to use fill
        
        title([CondNam{cond}, ' <--'])
    end
    plotind=plotind+2;
end
text(0, -1*max(MN(~isnan(MN))), sesNam)
subplot(5,2,9)
xlabel('Distance from teleport (cm)')
subplot(5,2,10)
xlabel('Distance from teleport (cm)')

OutputNam2=[sesNam,'Behavior','_meanDist','.pdf'];

f1.PaperUnits='Normalized';
f1.PaperPosition=[0 0 1 1];

saveas(f1,OutputNam1,'pdf')
f2.PaperUnits='Normalized';
f2.PaperPosition=[0 0 1 1];

%
saveas(f2,OutputNam2,'pdf')