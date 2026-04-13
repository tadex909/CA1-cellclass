%select_your data
fpath = 'N:\Vinca\MATLAB\Data\All_data';
main_folder='N:\Vinca\MATLAB\Data\'
h=readtable(strcat(main_folder,'select_sess.xlsx'))
keep_dat=h.session_name(logical(h.doselect))

%color infos
loadcolmathSFN
load colormaprom.mat

%load data
load(strcat(fpath, '\supool.mat'));

%groups
a17=find(supool.agegrp==17);
a18=find(supool.agegrp==18);
g1=sort([a17;a18]);
a19=find(supool.agegrp==19);
a20=find(supool.agegrp==20);
g2=sort([a19;a20]);
a21=find(supool.agegrp==21);
a22=find(supool.agegrp==22);
g3=sort([a21;a22]);
a23=find(supool.agegrp==23);
a24=find(supool.agegrp==24);
g4=sort([a23;a24]);
g5=find(supool.agegrp==25);

boxplot=0;
rankplot=1;
%Faire PO vs PNO plus tard quand supoolé plus de data. 
%Faire boucle avec 1 de g1à4 et j de g2à g5 pour toutes les comparaisons
%2à2
%% BOXPLOT

for i = 1:4;
    for j=2:5;
        if i>=j ;
            disp('already done')
        else
            
            grpi=eval(strcat('g',string(i)))
            grpj=eval(strcat('g',string(j)))
            
            if boxplot
fprm = fdt_fprmlist('boxplot');
fprm.x12 = {'PO17-18','PO19-20'};
fprm.x13 = {'PO17-18','PO21-22'};
fprm.x14 = {'PO17-18','PO23-24'};
fprm.x15 = {'PO17-18','PO25'};
fprm.x23 = {'PO19-20','PO21-22',};
fprm.x24 = {'PO19-20','PO23-24'};
fprm.x25 = {'PO19-20','PO25'};
fprm.x34 = {'PO21-22','PO23-24'};
fprm.x35 = {'PO21-22','PO25'};
fprm.x45 = {'PO23-24','PO25'};


%subplot(2, 5, 1)
figure()
fprm.ylab = 'Animal''s Velocity (cm /s)';
fprm.color = color_set;
fprm.ispaired = false;
%choose what to plot : 
fprm.x=eval(strcat('fprm.x',string(i),string(j)))


v1 = supool.vel_cond(grpi);
v2 = supool.vel_cond(grpj)
%v3 = supool.vel_cond(:, 3);
%v4 = supool.vel_cond(:, 4);
%v5 = supool.vel_cond(:, 5);
%v = [v1 v2 v3 v4 v5]
%boxplot(v)
fct_plot_boxplot (v1, v2 ,fprm)
%title(datalc.session_name, 'fontsize', 15)
%title(this_session, 'fontsize', 15)

saveas(gcf,strcat(fpath,'\velocity',string(fprm.x(1)),'_vs_',string(fprm.x(2)),'.png'));

%subplot(2, 5, 2)
figure()
fprm.ylab = 'Reward Frequency (rew /min)';
fprm.ispaired = true;
finderr=[grpi*2;grpi*2+1]
finder=[grpj*2;grpj*2+1]
 v2=supool.freq_rewardsess(finder)
v1=supool.freq_rewardsess(finderr)
% v1 = supool.freq_reward_cond(:, 1);
% v2 = supool.freq_reward_cond(:, 3);
%v3 = supool.freq_reward_cond(:, 3);
%v4 = supool.freq_reward_cond(:, 4);
%v5 = supool.freq_reward_cond(:, 5);
%fct_plot_boxplot(v1, v2,v3 ,v4, v5 ,fprm)
fct_plot_boxplot(v1, v2, fprm)


saveas(gcf,strcat(fpath,'\rew_fr',string(fprm.x(1)),'_vs_',string(fprm.x(2)),'.png'));


%subplot(2, 5, 3)
figure()
fprm.ylab = 'Active Cells (%)';
v1 = supool.pct_pyr_active_c(grpi);
v2 = supool.pct_pyr_active_c(grpj);
fct_plot_boxplot(v1, v2, fprm)
 %title(['N = ' num2str(allcel.nb_pyr_active_c(grpi)) ' vs ' num2str(allcel.nb_pyr_active_c(grpj))])
saveas(gcf,strcat(fpath,'active_cell',string(fprm.x(1)),'_vs_',string(fprm.x(2)),'.png'));

%subplot(2, 5, 4)
figure()
fprm.ylab = 'Place Cells (%)';
v1 = supool.pct_pyr_sm_c(grpi);
v2 = supool.pct_pyr_sm_c(grpj);
fct_plot_boxplot(v1, v2, fprm)
% title(['N = ' num2str(allcel.nb_pyr_sm_c(1)) ' vs ' num2str(allcel.nb_pyr_sm_c(2))])
saveas(gcf,strcat(fpath,'\%PC',string(fprm.x(1)),'_vs_',string(fprm.x(2)),'.png'));

% fprm.ylab = 'Bidirectional Cells (%)';
% v1 = supool.pct_pyr_bid_c(:, 1);
% v2 = supool.pct_pyr_bid_c(:, 2);
% fct_plot_boxplot(v1, v2, fprm)
% 
% fprm.ylab = 'Unidirectional Cells (%)';
% v1 = supool.pct_pyr_uni_c(:, 1);
% v2 = supool.pct_pyr_uni_c(:, 2);
% fct_plot_boxplot(v1, v2, fprm)
cellidi=[];
cellidj=[];
for k =1:length(grpi);
    idi=find(supool.keptsession_idcel==grpi(k));
    cellidi=sort([cellidi;idi]);
end
for k =1:length(grpj);
    idj=find(supool.keptsession_idcel==grpj(k));
    cellidj=sort([cellidj;idi]);
end
%subplot(2, 5, 5)
figure()
fprm.ylab = 'Firing Rate (Hz)';
v1 = supool.fr_uc(find(cellidi));
v2 = supool.fr_uc(find(cellidj));
fct_plot_boxplot(v1, v2, fprm)
%subplot(2, 5, 6)
saveas(gcf,strcat(fpath,'\FR',string(fprm.x(1)),'_vs_',string(fprm.x(2)),'.png'));


figure()
fprm.ylab = 'Out/In Field Firing Rate';
outovin_ucf = supool.outovin_ucf(:, :, 1);
v1 = outovin_ucf(find(cellidi));
v2 = outovin_ucf(find(cellidj));
fct_plot_boxplot(v1, v2, fprm)
saveas(gcf,strcat(fpath,'\out-in',string(fprm.x(1)),'_vs_',string(fprm.x(2)),'.png'));


%subplot(2, 5, 7)
figure()
fprm.ylab = 'PF Width (cm)';
size_ucf = 2*supool.size_ucf(:, :, 1);
v1 = size_ucf(find(cellidi));
v2 = size_ucf(find(cellidj));
fct_plot_boxplot(v1, v2, fprm)
saveas(gcf,strcat(fpath,'\width',string(fprm.x(1)),'_vs_',string(fprm.x(2)),'.png'));


%subplot(2, 5, 8)
figure()
fprm.ylab = 'PF Deviation (cm)';
lap_devwithmean_ucf = 2*supool.lap_devwithmean_ucf(:, :, 1);
v1 = lap_devwithmean_ucf(find(cellidi));
v2 = lap_devwithmean_ucf(find(cellidj));
fct_plot_boxplot(v1, v2, fprm)
saveas(gcf,strcat(fpath,'\deviation',string(fprm.x(1)),'_vs_',string(fprm.x(2)),'.png'));


%subplot(2, 5, 9)
figure()
fprm.ylab = 'Spatial Information (bit/spk)';
v1 = supool.si_uc(find(cellidi));
v2 = supool.si_uc(find(cellidj));
fct_plot_boxplot(v1, v2, fprm)
saveas(gcf,strcat(fpath,'\SI',string(fprm.x(1)),'_vs_',string(fprm.x(2)),'.png'));


%subplot(2, 5, 10)
figure()
fprm.ylab = 'Stability';
v1 = supool.stb_allpairs_uc(find(cellidi));
v2 = supool.stb_allpairs_uc(find(cellidj));


fct_plot_boxplot(v1, v2, fprm)
saveas(gcf,strcat(fpath,'\stability',string(fprm.x(1)),'_vs_',string(fprm.x(2)),'.png'));


set(gcf, 'color', 'w')

fct_fullscreen(gcf)

close all

%fct_save_figure(gcf, [datalc.save_path filesep 'Pooled_BoxPlot' i 'vs' j], 'jpg', 'eps')
            end
        end
    end

          
%% Slect data for rankfield and plot by age
if rankplot
load colormaprom.mat
allrmap.feat.ispyr_uc = kron(supool.pyr_uc, ones(1, 2));
ispyrpf = supool.ispf_ucf(:, :, 1) & allrmap.feat.ispyr_uc;

% supool_g1
for i = 1:5;
            grpi=eval(strcat('g',string(i)))
cellidi=[];

for k =1:length(grpi);
    idi=find(supool.keptsession_idcel==grpi(k));
    cellidi=sort([cellidi;idi]);
end
ispyrpftemp=ispyrpf(cellidi,1:2) 
ispyrtemp = kron(supool.pyr_uc, ones(1, 2))
% rankfields 
load colormaprom.mat

allrmap.feat.ispyr_uc = kron(supool.pyr_uc, ones(1, 2));
ispyrpf = supool.ispf_ucf(:, :, 1) & allrmap.feat.ispyr_uc;
fr_s_cxutemp=supool.fr_s_cxu(1:2,:,cellidi);
% nb_cond =length( supool.nbcond)*2
nb_cond =2
nb_cel = length(supool.session_idcel(cellidi));


A = cell(nb_cond, 1);
B = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
             if ispyrtemp(u, c)  & c==1%si c une place cell 
        A{c} = [A{c} ; fr_s_cxutemp(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
%          A{c+1} = [A{c+1} ; supool.fr_s_cxu(c+1, : , u)]; 
         elseif ispyrtemp(u, c)  & c==2%             A{c-1} = [A{c-1} ; fr_s_cxutemp(c-1, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
           A{c-1} = [A{c-1} ; supool.fr_s_cxu(c, : , u)]; 
           
%            A{c} = [A{c} ; fr_s_cxutemp(c, : , u)];
%         if ispyrpftemp(u, c)  & c==2 % si il s'agit d'une condition différente de la cond 1 (condw 1et2)
%              A{c} = [A{c} ; fr_s_cxutemp(c, : , u)]; 
% %         B{c} = [B{c} ; supool.fr_s_cxu(c, : , u)]; 
%         B{c + 2}  = [B{c + 2} ; supool.fr_s_cxu(c + 2, : , u)]; 
        end
        end
end

    if find(A{1}) | find(B{1})
%
fprm = fdt_fprmlist_ball('rankfield');
fprm.xt = fprm.xt - 10;
fprm.xlim = fprm.xlim - 10;
fprm.ylab = 'Place Field Number';
edgobj = [15 36 ; 85 90] - 10;

figure
subplot(221)
fprm.colormap = colormapbnl;
fct_plot_rankfield([A{1}], fprm)
%title(datalc.session_name, 'fontsize', 15)
title('PO Aller & retours 17_18')
subplot(223)
v = nanmean(fct_matnorm([A{1}]), 1);
plot(v)
hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
ylim=[0,1]
hold off


subplot(222)
fprm.colormap = colormapcnl;
fct_plot_rankfield([A{2}], fprm)
title('PO Retour')


subplot(224)
hold off
v = nanmean(fct_matnorm([A{2}]), 1);
plot(v)

hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');



saveas(gcf,strcat(fpath,(strcat('g',string(i))),'allAR_allpyr''.png'));
    end
end
end
end
%% rankfields essai vinca

 allrmap.feat.ispyr_uc = kron(supool.pyr_uc, ones(1, 2));
ispyrpf = supool.ispf_ucf(:, 1:2, 1) & allrmap.feat.ispyr_uc(:,1:2);

%  nb_cond = (length(supool.nbcond(:,1)));  OLD
nb_cond = (length(find(~isnan(supool.icond)))); %try
nb_cel = length(supool.session_idcel)


A = cell(nb_cond, 2);

 matcondw=supool.icond(~isnan(supool.icond));
supool.condw=[matcondw,matcondw];
 A(:)=num2cell(supool.condw');
%A(:,2)=num2cell(1)
B = cell(nb_cond, 2);

for u = 1:nb_cel
    for c = 1:2
        if ispyrpf(u, c) % si la cellule dans cette condition est un place cell ...
        A{c} = [ supool.fr_s_cxu(c, : , u)];  %l'index de cette condition dans A se voit append le firing rate de cette cellule pour cette session
       
        if c < 3  % tant que c est inférieur à 3 (donc = 0 1 ou 2) ...
        B{c} = [B{c} ; supool.fr_s_cxu(c, : , u)]; %l'index de cette condition dans B se voit append le firing rate de cette cellule pour cette session
        B{c + 2}  = [B{c + 2} ; supool.fr_s_cxu(c + 2, : , u)]; %l'index de valeur cette condition plus 2 se voit append le firing rate de cette +2 condition à B
        end
        end
    end
end

fprm = fdt_fprmlist('rankfield');
fprm.xt = fprm.xt - 10;
fprm.xlim = fprm.xlim - 10;
fprm.ylab = 'Place Field Number';
edgobj = [15 35 ; 80 90] - 10;

%créer A(condPO)w and b puis plot cela en 221

subplot(221)
fprm.colormap = colormapbnl;
fct_plot_rankfield([A{1} ; A{2}], fprm) % utiliser un n allant de 1 à 10 avec A(1) A(2)^puis A(3) A(4) et mettre une légende = à condw?
fct_plot_rankfield([B{1} ; B{2}], fprm) 
%title(datalc.session_name, 'fontsize', 15)

subplot(223)
v = nanmean(fct_matnorm([A{1} ; A{2}]), 1);
plot(v)
hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');


