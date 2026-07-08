%select_your data
fpath = 'X:\MATLAB\Data\All_data_final\process';
main_folder='X:\MATLAB\Data\All_data_final\'
h=readtable(strcat(main_folder,'Finalyses_Vinca.xlsx'))
keep_dat=h.session_name(logical(h.doselect_NOW))

%color infos
loadcolmathSFN
load colormaprom.mat

%load data
load(strcat(fpath, '\supool_final_all_newclass_halfwidth.mat'));

%groups
% 
% a15=find(supool.agegrp==15);
% a16=find(supool.agegrp==16);
% g1=sort([a15;a16]);
% a17=find(supool.agegrp==17);
% a18=find(supool.agegrp==18);
% g2=sort([a17;a18]);
% a19=find(supool.agegrp==19);
% a20=find(supool.agegrp==20);
% g3=sort([a19;a20]);
% a21=find(supool.agegrp==21);
% a22=find(supool.agegrp==22);
% g4=sort([a21;a22]);
% a23=find(supool.agegrp==23);
% a24=find(supool.agegrp==24);
% g5=sort([a23;a24]);
% a25=find(supool.agegrp==25);
% g6=sort([a25]);
%% 4 GROUPS

a15=find(supool.agegrp==15);
a16=find(supool.agegrp==16);
a17=find(supool.agegrp==17);

g1=sort([a15;a16;a17]);

a18=find(supool.agegrp==18);
a19=find(supool.agegrp==19);
a20=find(supool.agegrp==20);
g2=sort([a18;a19;a20]);

a22=find(supool.agegrp==22);
a21=find(supool.agegrp==21);
g3=sort([a21;a22]);

a23=find(supool.agegrp==23);
a24=find(supool.agegrp==24);
a25=find(supool.agegrp==25);
g4=sort([a23;a24;a25]);
%Faire PO vs PNO plus tard quand supoolé plus de data. 
%Faire boucle avec 1 de g1à4 et j de g2à g5 pour toutes les comparaisons
%2à2
%% BOXPLOT

%% Slect data for rankfield and plot by age

load colormaprom.mat
allrmap.feat.ispyr_uc = kron(supool.pyr_uc, ones(1, 2));
ispyrpf = supool.ispf_ucf(:, :, 1) & allrmap.feat.ispyr_uc;

%% PYRAMIDALES
for i = 1:4;
            gi=eval(strcat('g',string(i)))
cellidi=[];

for k =1:length(gi);
    idi=find(supool.keptsession_idcel==gi(k));
    cellidi=sort([cellidi;idi]);
end

% rankfields 
load colormaprom.mat

allrmap.feat.ispyr_uc = kron(supool.pyr_uc, ones(1, 2));
ispyrpf = supool.ispf_ucf(:, :, 1) & allrmap.feat.ispyr_uc;
fr_s_cxutemp=supool.fr_s_cxu(1:2,:,cellidi);
ispyrpftemp=ispyrpf(cellidi,1:2) 
ispyrtemp = kron(supool.pyr_uc, ones(1, 2))
% nb_cond =length( supool.nbcond)*2
nb_cond =2
nb_cel = length(supool.session_idcel(cellidi));


A = cell(nb_cond, 1);
B = cell(nb_cond, 1);

A = cell(nb_cond, 1);
B = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
             if ispyrtemp(u, c)  & c==1%si c une place cell 
        A{c} = [A{c} ; fr_s_cxutemp(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
%         A{c+1} = [A{c+1} ; supool.fr_s_cxu(c+1, : , u)]; 
         elseif ispyrtemp(u, c)  & c==2%             A{c-1} = [A{c-1} ; fr_s_cxutemp(c-1, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
           A{c-1} = [A{c-1} ; supool.fr_s_cxu(c, : , u)]; 
             end
    end
end

%
fprm = fdt_fprmlist('rankfield');
fprm.xt = fprm.xt - 10;
fprm.xlim = fprm.xlim - 10;
fprm.ylab = 'Cell Number';
edgobj = [15 35 ; 80 90] -10 ;

figure
subplot(211)
fprm.colormap = colormapbnl;
fct_plot_rankfield([A{1}], fprm)
%title(datalc.session_name, 'fontsize', 15)
title(strcat('PO Aller et Retour Pyr only Groupe',num2str(i)))
subplot(212)
v = nanmean(fct_matnorm([A{1}]), 1);
plot(v)
hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
ylim=[0,1]
hold off

% 
% subplot(222)
% fprm.colormap = colormapcnl;
% fct_plot_rankfield([A{2}], fprm)
% title('PO Retour')
% 
% 
% subplot(224)
% hold off
% v = nanmean(fct_matnorm([A{2}]), 1);
% plot(v)
% 
% hold on
% line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
% line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
% line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
% line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');



 saveas(gcf,strcat(fpath,(strcat('\g',string(i))),'allAR_allpyr_sm_and_nonsm_4groups''.png'));
close all
end
%% ALL PF IN et PC
for i = 1:4;
    
            gi=eval(strcat('g',string(i)))
cellidi=[];

for k =1:length(gi);
    idi=find(supool.keptsession_idcel==gi(k));
    cellidi=sort([cellidi;idi]);
end

ispf = supool.ispf_ucf(:, :, 1) ;
fr_s_cxutemp=supool.fr_s_cxu(1:2,:,cellidi);
ispftemp=ispf(cellidi,1:2) 
% nb_cond =length( supool.nbcond)*2
nb_cond =2
nb_cel = length(supool.session_idcel(cellidi));


A = cell(nb_cond, 1);
B = cell(nb_cond, 1);

A = cell(nb_cond, 1);
B = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
             if ispftemp(u, c)  & c==1%si c une place cell 
        A{c} = [A{c} ; fr_s_cxutemp(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
%         A{c+1} = [A{c+1} ; supool.fr_s_cxu(c+1, : , u)]; 
         elseif ispftemp(u, c)  & c==2%             A{c-1} = [A{c-1} ; fr_s_cxutemp(c-1, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
           A{c-1} = [A{c-1} ; supool.fr_s_cxu(c, : , u)]; 
             end
    end
end

%
fprm = fdt_fprmlist('rankfield');
fprm.xt = fprm.xt - 10;
fprm.xlim = fprm.xlim - 10;
fprm.ylab = 'Cell Number';
edgobj = [15 35 ; 80 90] -10 ;

figure
subplot(211)
fprm.colormap = colormapbnl;
fct_plot_rankfield([A{1}], fprm)
%title(datalc.session_name, 'fontsize', 15)
title(strcat('PO Aller et Retour all pf (IN and PC) Groupe',num2str(i)))
subplot(212)
v = nanmean(fct_matnorm([A{1}]), 1);
plot(v)
hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
ylim=[0,1]
hold off

 saveas(gcf,strcat(fpath,(strcat('\g',string(i))),'allAR_PCandIN_4groups''.png'));
  close all
end
%% INTERNEURONSS PF
a=kron(supool.pyr_uc, ones(1, 2))
isIN=a(:,1:2)+1
isIN(isIN==2)=0
isIN=[isIN,zeros(size(isIN,1),8)]
isINpf=isIN&supool.ispf_ucf(:, :, 1) 
for i = 1:4;
    
            gi=eval(strcat('g',string(i)))
cellidi=[];

for k =1:length(gi);
    idi=find(supool.keptsession_idcel==gi(k));
    cellidi=sort([cellidi;idi]);
end

fr_s_cxutemp=supool.fr_s_cxu(1:2,:,cellidi);
isINpftemp=isINpf(cellidi,1:2) 
% nb_cond =length( supool.nbcond)*2
nb_cond =2
nb_cel = length(supool.session_idcel(cellidi));


A = cell(nb_cond, 1);
B = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
             if isINpftemp(u, c)  & c==1%si c une place cell 
        A{c} = [A{c} ; fr_s_cxutemp(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
%         A{c+1} = [A{c+1} ; supool.fr_s_cxu(c+1, : , u)]; 
         elseif isINpftemp(u, c)  & c==2
%             A{c-1} = [A{c-1} ; fr_s_cxutemp(c-1, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
           A{c-1} = [A{c-1} ; supool.fr_s_cxu(c, : , u)]; 
             end
    end
end

%
fprm = fdt_fprmlist('rankfield');
fprm.xt = fprm.xt - 10;
fprm.xlim = fprm.xlim - 10;
fprm.ylab = 'Cell Number';
edgobj = [15 35 ; 80 90] -10 ;

figure
subplot(211)
fprm.colormap = colormapbnl;
fct_plot_rankfield([A{1}], fprm)
%title(datalc.session_name, 'fontsize', 15)
title(strcat('PO Aller et Retour all IN PF Only Groupe',num2str(i)))
subplot(212)
v = nanmean(fct_matnorm([A{1}]), 1);
plot(v)
hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
ylim=[0,1]
hold off

 saveas(gcf,strcat(fpath,(strcat('\g',string(i))),'allAR_IN_Only_4groups''.png'));
  close all
end
%% pyra pc aller vs retours

load colormaprom.mat
allrmap.feat.ispyr_uc = kron(supool.pyr_uc, ones(1, 2));
ispyrpf = supool.ispf_ucf(:, :, 1) & allrmap.feat.ispyr_uc;


for i = 1:4;
            gi=eval(strcat('g',string(i)))
cellidi=[];

for k =1:length(gi);
    idi=find(supool.keptsession_idcel==gi(k));
    cellidi=sort([cellidi;idi]);
end

% rankfields 
load colormaprom.mat

allrmap.feat.ispyr_uc = kron(supool.pyr_uc, ones(1, 2));
ispyrpf = supool.ispf_ucf(:, :, 1) & allrmap.feat.ispyr_uc;
fr_s_cxutemp=supool.fr_s_cxu(1:2,:,cellidi);
ispyrpftemp=ispyrpf(cellidi,1:2) 
ispyrtemp = kron(supool.pyr_uc, ones(1, 2))
% nb_cond =length( supool.nbcond)*2
nb_cond =2
nb_cel = length(supool.session_idcel(cellidi));



A = cell(nb_cond, 1);
B = cell(nb_cond, 1);

A = cell(nb_cond, 1);
B = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
             if ispyrpftemp(u, c)  & c==1%si c une place cell 
        A{c} = [A{c} ; fr_s_cxutemp(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
%         A{c+1} = [A{c+1} ; supool.fr_s_cxu(c+1, : , u)]; 
         elseif ispyrpftemp(u, c)  & c==2%             A{c-1} = [A{c-1} ; fr_s_cxutemp(c-1, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
           A{c-1} = [A{c-1} ; supool.fr_s_cxu(c, : , u)]; 
             end
    end
end

%
fprm = fdt_fprmlist('rankfield');
fprm.xt = fprm.xt - 10;
fprm.xlim = fprm.xlim - 10;
fprm.ylab = 'Place Field Number';
edgobj = [15 35 ; 80 90] -10 ;

figure
subplot(211)
fprm.colormap = colormapbnl;
fct_plot_rankfield([A{1}], fprm)
%title(datalc.session_name, 'fontsize', 15)
title(strcat('PO Aller et Retour Place cells Groupe',num2str(i)))
subplot(212)
v = nanmean(fct_matnorm([A{1}]), 1);
plot(v)
hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
ylim=[0,1]
hold off

% 
% subplot(222)
% fprm.colormap = colormapcnl;
% fct_plot_rankfield([A{2}], fprm)
% title('PO Retour')
% 
% 
% subplot(224)
% hold off
% v = nanmean(fct_matnorm([A{2}]), 1);
% plot(v)
% 
% hold on
% line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
% line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
% line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
% line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');



 saveas(gcf,strcat(fpath,(strcat('\g',string(i))),'allAR_allPF_CSI_4groups''.png'));
  close all
%   figure
% for n=1:size(A{1},1)
% subplot(size(A{1},1),1,n)
% plot(A{1}(n,:))
% end
%%


allrmap.feat.ispyr_uc = kron(supool.pyr_uc, ones(1, 2));
ispyrpf = supool.ispf_ucf(:, :, 1) & allrmap.feat.ispyr_uc;
% ispyrpf = supool.ispf_ucf(:, :, 1) & (allrmap.feat.ispyr_uc==0);
fr_s_cxutemp=supool.fr_s_cxu(1:2,:,cellidi);
ispyrpftemp=ispyrpf(cellidi,1:2) 
ispyrtemp = kron(supool.pyr_uc, ones(1, 2))
% nb_cond =length( supool.nbcond)*2
nb_cond =2
nb_cel = length(supool.session_idcel(cellidi));


A = cell(nb_cond, 1);
B = cell(nb_cond, 1);

A = cell(nb_cond, 1);
B = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
             if ispyrpftemp(u, c)  & c==1%si c une place cell 
        A{c} = [A{c} ; fr_s_cxutemp(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
%         A{c+1} = [A{c+1} ; supool.fr_s_cxu(c+1, : , u)]; 
         elseif ispyrpftemp(u, c)  & c==2%             A{c-1} = [A{c-1} ; fr_s_cxutemp(c-1, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
           A{c} = [A{c} ; supool.fr_s_cxu(c, : , u)]; 
             end
    end
end

%
fprm = fdt_fprmlist('rankfield');
fprm.xt = fprm.xt - 10;
fprm.xlim = fprm.xlim - 10;
fprm.ylab = 'Place Field Number';
edgobj = [15 35 ; 80 90] -10 ;

figure
subplot(221)
fprm.colormap = colormapbnl;
fct_plot_rankfield([A{1}], fprm)
%title(datalc.session_name, 'fontsize', 15)
title(strcat('PO Aller et Retour Place cells Groupe',num2str(i)))
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
fprm.colormap = colormapbnl;
fct_plot_rankfield([A{2}], fprm) % utiliser un n allant de 1 à 10 avec A(1) A(2)^puis A(3)A(4) et mettre une légende = à condw? fct_plot_rankfield([B{1} ; B{2}],fprm) %title(datalc.session_name, 'fontsize', 15)

subplot(224) 
v = nanmean(fct_matnorm([A{2}]), 1); plot(v) 
hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'r','linestyle', '--'); 
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)],'color', 'r', 'linestyle', '--'); 
line([edgobj(2, 1) edgobj(2, 1)],[min(v) max(v)], 'color', 'r', 'linestyle', '--'); 
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');

hold off
saveas(gcf,strcat(fpath,(strcat('\g',string(i))),'A&R_allpc_4groups''.png'));
close all

end
% % 

% 
% 
% rankfields essai vinca

 allrmap.feat.ispyr_uc = kron(supool.pyr_uc, ones(1, 2));
ispyrpf = supool.ispf_ucf(:, 1:2, 1) & allrmap.feat.ispyr_uc(:,1:2);

 nb_cond = (length(supool.nbcond(:,1)));  %OLD
nb_cond = (length(find(~isnan(supool.icond)))); %try
nb_cel = length(supool.session_idcel)


A = cell(nb_cond, 2);

 matcondw=supool.icond(~isnan(supool.icond));
supool.condw=[matcondw,matcondw];
 A(:)=num2cell(supool.condw');
A(:,2)=num2cell(1)
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
fprm.ylab = ' Place Field Number';
edgobj = [15 35 ; 80 90] - 10;

% créer A(condPO)w and b puis plot cela en 221

subplot(221)
fprm.colormap = colormapbnl;
fct_plot_rankfield([A{1} ; A{2}], fprm) % utiliser un n allant de 1 à 10 avec A(1) A(2)^puis A(3) A(4) et mettre une légende = à condw?
fct_plot_rankfield([B{1} ; B{2}], fprm) 
% title(datalc.session_name, 'fontsize', 15)

subplot(223)
v = nanmean(fct_matnorm([A{1} ; A{2}]), 1);
plot(v)
hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');


