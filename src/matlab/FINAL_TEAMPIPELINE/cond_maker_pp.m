%% get separation between sessions S1-S2-S3

hf=figure;
ax1=subplot(2,1,1);
plot(Sess,'k')
hold on
set(hf,'Units','normalized');
set(hf, 'Position',[0.05 0.05 0.8 0.5])

NSess=str2double(cell2mat(inputdlg('How Many Sessions?','Sessions',1,{'1'})));
prompts=cell(1,NSess);

prompts(1:end)=compose('Name for Sess %d', 1:NSess);
Defaults=cell(1,NSess);
Defaults(1:end)={'Obj'};
CondNam=inputdlg(prompts,'Names',1,Defaults);

close(hf)

hf=figure;
set(hf,'Units','normalized');
set(hf, 'Position',[0.05 0.05 0.8 0.5])

times=zeros(NSess,2);
sessS=zeros(NSess,1);

t2=t-t(1);

for c=1:NSess%tryin to plot a peace at a time, for precision
    idx=find(t2>times(c,1)& t2<times(c,1)+(900));
    plot(t2(idx),Sess(idx))
    xlim ([t2(idx(1)),t2(idx(end))])
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