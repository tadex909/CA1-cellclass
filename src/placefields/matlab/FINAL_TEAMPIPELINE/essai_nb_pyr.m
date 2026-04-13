%select_your data
fpath = 'N:\Vinca\MATLAB\Data\All_data';
main_folder='N:\Vinca\MATLAB\Data\'
h=readtable(strcat(main_folder,'select_sess.xlsx'))
keep_dat=h.session_name(logical(h.doselect))

%color infos
loadcolmathSFN
load colormaprom.mat

%load data
load(strcat(fpath, '\supool_v2.mat'));

%groups

a15=find(supool.agegrp==15);
a16=find(supool.agegrp==16);
g1=sort([a15;a16]);
a17=find(supool.agegrp==17);
a18=find(supool.agegrp==18);
g2=sort([a17;a18]);
a19=find(supool.agegrp==19);
a20=find(supool.agegrp==20);
g3=sort([a19;a20]);
a21=find(supool.agegrp==21);
a22=find(supool.agegrp==22);
g4=sort([a21;a22]);
a23=find(supool.agegrp==23);
a24=find(supool.agegrp==24);
g5=sort([a23;a24]);
a25=find(supool.agegrp==25);
g6=sort([a25]);
boxplots=0;
rankplot=1;



load colormaprom.mat
allrmap.feat.ispyr_uc = kron(supool.pyr_uc, ones(1, 2));
ispyrpf = supool.ispf_ucf(:, :, 1) & allrmap.feat.ispyr_uc;
%% ALLER ET RETOURS SEPARES
% supool_g1
for i =[1,2,3,4,5,6];
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
% ispf_cxutemp = supool.ispf_cxu(1:2,:, cellidi)
bin_nb_pyr=[]
bin_nb_pyr_b=[]
% nb_cond =length( supool.nbcond)*2
nb_cond =2
nb_cel = length(supool.session_idcel(cellidi));
n=1

A = cell(nb_cond, 1);
B = cell(nb_cond, 1);

A = cell(nb_cond, 1);
B = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
             if ispyrtemp(u, c)  & c==1%si c une place cell 
                A{c} = [A{c} ; fr_s_cxutemp(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
                      bin_nb_pyr(n,1:80)=zeros(1,80);
                    
                      
                     bin_nb_pyr(n,find(fr_s_cxutemp(c,:,u)==max(fr_s_cxutemp(c,:,u))))=1;



               n=n+1
%         A{c+1} = [A{c+1} ; supool.fr_s_cxu(c+1, : , u)]; 
         elseif ispyrtemp(u, c)  & c==2%             A{c-1} = [A{c-1} ; fr_s_cxutemp(c-1, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
           A{c} = [A{c} ; supool.fr_s_cxu(c, : , u)]; 
                      bin_nb_pf_b(n,1:80)=zeros(1,80);

                     bin_nb_pf_b(n,find(fr_s_cxutemp(c,:,u)==max(fr_s_cxutemp(c,:,u))))=1;

           n=n+1
             end
    end
end
fprm = fdt_fprmlist('rankfield');
fprm.xt = fprm.xt - 10;
fprm.xlim = fprm.xlim - 10;
fprm.ylab = 'Place Field Number';
edgobj = [15 35 ; 80 90] -10 ;

figure
subplot(321)
fprm.colormap = colormapbnl;
fct_plot_rankfield([A{1}], fprm)
%title(datalc.session_name, 'fontsize', 15)
title(strcat('PO Way Pyramidal cells Groupe',num2str(i)))
hold on
line([edgobj(1, 1) edgobj(1, 1)], [0 size(A{1},1)+0.5], 'color', 'm', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [0 size(A{1},1)+0.5], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [0 size(A{1},1)+0.5], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [0 size(A{1},1)+0.5], 'color', 'm', 'linestyle', '--');
ylim=[0,1]
hold off

subplot(322)
fprm.colormap = colormapbnl;
fct_plot_rankfield([A{2}], fprm)
%title(datalc.session_name, 'fontsize', 15)
title(strcat('PO Back Pyramidal cells Groupe',num2str(i)))
hold on
line([edgobj(1, 1) edgobj(1, 1)], [0 size(A{2},1)+0.5], 'color', 'm', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [0 size(A{2},1)+0.5], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [0 size(A{2},1)+0.5], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [0 size(A{2},1)+0.5], 'color', 'm', 'linestyle', '--');
ylim=[0,1]
hold off

subplot(323)
v = nanmean(fct_matnorm([A{1}]), 1);
plot(v)
hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
ylim=[0,1]
hold off

A_1=[]
[A_n, mat_divison] = fct_matnorm([A{1}], 2);
[ind_perm, ~] = find(abs([A{1}] - mat_divison) < 5*eps);
A_np = A_n(ind_perm, :);
for n=1:length(ind_perm)
A_1(n,1:80)=zeros(1,80)
A_1(n,find(A_np(n,:)==1))=1
end
A_2=[]
[A_n, mat_divison] = fct_matnorm([A{2}], 2);
[ind_perm, ~] = find(abs([A{2}] - mat_divison) < 5*eps);
A_np = A_n(ind_perm, :);
for n=1:length(ind_perm)
A_2(n,1:80)=zeros(1,80)
A_2(n,find(A_np(n,:)==1))=1
end

bin_allpf=sum(A_1);
bin_allpf_b=sum(A_2);
sumbin=[];
sumbin_b=[];
% if sum(ispyrpftemp)>=2
for j=1:10:71
   sumbin(j:j+9)=sum(bin_allpf(j:j+9))
   sumbin_b(j:j+9)=sum(bin_allpf_b(j:j+9))
end
% end

subplot(325)
stairs(sumbin)
% plot(bin_allpf)
hold on
line([edgobj(1, 1) edgobj(1, 1)], [0 max(sumbin)+0.1], 'color', 'm', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [0 max(sumbin)+0.1], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [0 max(sumbin)+0.1], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [0 max(sumbin)+0.1], 'color', 'm', 'linestyle', '--');
ylim=[0,1]
hold off

subplot(324)
v = nanmean(fct_matnorm([A{2}]), 1);
plot(v)
hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
ylim=[0,1]
hold off

subplot(326)
stairs(sumbin_b)
hold on
line([edgobj(1, 1) edgobj(1, 1)], [0 max(sumbin_b)+0.1], 'color', 'm', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [0 max(sumbin_b)+0.1], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [0 max(sumbin_b)+0.1], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [0 max(sumbin_b)+0.1], 'color', 'm', 'linestyle', '--');
ylim=[0,1]
% hold off
% % plot(bin_allpf_b)
% subplot(326)
% fprm = fdt_fprmlist('histline');
% fprm.xlim = [10.5 89.5];
% fprm.ylim = [-0.5 6];
% fprm.xlab = '';
% fprm.ylab = 'Place Fields (%)';
% fprm.xt = [10.5 30 50 70 89.5];
% fprm.xtl = 20:40:180;
% fct_plot_hist_line( sumbin_b,5, fprm);
% subplot(326)
% fprm = fdt_fprmlist('histline');
% fprm.xlim = [10.5 80.5];
% fprm.ylim = [-0.5 4];
% fprm.xlab = '';
% fprm.ylab = 'Place Fields (%)';
% fprm.xt = [10.5 30 50 70 89.5];
% fprm.xtl = 20:40:180;
% fct_plot_hist_line( sumbin_b,5, fprm);
% hold on
% line([edgobj(1, 1) edgobj(1, 1)], [0 5], 'color', 'm', 'linestyle', '--');
% line([edgobj(1, 2) edgobj(1, 2)], [0 5], 'color', 'm', 'linestyle', '--');
% line([edgobj(2, 1) edgobj(2, 1)], [0 5], 'color', 'm', 'linestyle', '--');
% line([edgobj(2, 2) edgobj(2, 2)], [0 5], 'color', 'm', 'linestyle', '--');
% ylim=[0,1]
% hold off
end

%% Allers et retours ensembles

for i =[1,2,3,4,5,6];
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
% ispyrpftemp=ispyrpf(cellidi,1:2) 
ispyrtemp = kron(supool.pyr_uc, ones(1, 2))
% ispf_cxutemp = supool.ispf_cxu(1:2,:, cellidi)
bin_nb_pyr=[]
bin_nb_pyr_b=[]
% nb_cond =length( supool.nbcond)*2
nb_cond =2
nb_cel = length(supool.session_idcel(cellidi));
n=1

A = cell(nb_cond, 1);
B = cell(nb_cond, 1);

A = cell(nb_cond, 1);
B = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
             if ispyrtemp(u, c)  & c==1%si c une place cell 
                A{c} = [A{c} ; fr_s_cxutemp(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
                      bin_nb_pyr(n,1:80)=zeros(1,80);
                     bin_nb_pyr(n,find(fr_s_cxutemp(c,:,u)==max(fr_s_cxutemp(c,:,u))))=1;
               n=n+1
%         A{c+1} = [A{c+1} ; supool.fr_s_cxu(c+1, : , u)]; 
         elseif ispyrtemp(u, c)  & c==2%             A{c-1} = [A{c-1} ; fr_s_cxutemp(c-1, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
           A{c-1} = [A{c-1} ; supool.fr_s_cxu(c, : , u)]; 
                      bin_nb_pyr(n,1:80)=zeros(1,80);
                     bin_nb_pyr(n,find(fr_s_cxutemp(c,:,u)==max(fr_s_cxutemp(c,:,u))))=1;
       
           n=n+1
             end
    end
end
fprm = fdt_fprmlist('rankfield');
fprm.xt = fprm.xt - 10;
fprm.xlim = fprm.xlim - 10;
fprm.ylab = 'Cell Number';
edgobj = [15 35 ; 80 90] -10 ;

figure
subplot(311)
fprm.colormap = colormapbnl;
fct_plot_rankfield([A{1}], fprm)
%title(datalc.session_name, 'fontsize', 15)
title(strcat('PO Pyramidal cells Groupe',num2str(i)))
hold on
line([edgobj(1, 1) edgobj(1, 1)], [0.5 size(A{1},1)+0.5], 'color', 'm', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [0 size(A{1},1)+0.5], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [0 size(A{1},1)+0.5], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [0 size(A{1},1)+0.5], 'color', 'm', 'linestyle', '--');
ylim=[0,1]
hold off
% 
% subplot(312)
% fprm.colormap = colormapbnl;
% fct_plot_rankfield([A{2}], fprm)
% %title(datalc.session_name, 'fontsize', 15)
% title(strcat('PO Retours Place cells Groupe',num2str(i)))
% hold on
% line([edgobj(1, 1) edgobj(1, 1)], [0 size(A{2},1)+0.5], 'color', 'm', 'linestyle', '--');
% line([edgobj(1, 2) edgobj(1, 2)], [0 size(A{2},1)+0.5], 'color', 'm', 'linestyle', '--');
% line([edgobj(2, 1) edgobj(2, 1)], [0 size(A{2},1)+0.5], 'color', 'm', 'linestyle', '--');
% line([edgobj(2, 2) edgobj(2, 2)], [0 size(A{2},1)+0.5], 'color', 'm', 'linestyle', '--');
% ylim=[0,1]
% hold off

subplot(312)
v = nanmean(fct_matnorm([A{1}]), 1);
plot(v)
hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
ylim=[0,1]
hold off

A_1=[]
[A_n, mat_divison] = fct_matnorm([A{1}], 2);
[ind_perm, ~] = find(abs([A{1}] - mat_divison) < 5*eps);
A_np = A_n(ind_perm, :);
for n=1:length(ind_perm)
A_1(n,1:80)=zeros(1,80)
A_1(n,find(A_np(n,:)==1))=1
end
% A_2=[]
% [A_n, mat_divison] = fct_matnorm([A{2}], 2);
% [ind_perm, ~] = find(abs([A{2}] - mat_divison) < 5*eps);
% A_np = A_n(ind_perm, :);
% for n=1:length(ind_perm)
% A_2(n,1:80)=zeros(1,80)
% A_2(n,find(A_np(n,:)==1))=1
% end

bin_allpf=sum(A_1);
% bin_allpf_b=sum(A_2);
sumbin=[];
% sumbin_b=[];
% if sum(ispyrpftemp)>=2
for j=1:10:71
   sumbin(j:j+9)=sum(bin_allpf(j:j+9))
%    sumbin_b(j:j+9)=sum(bin_allpf_b(j:j+9))
end
% end

subplot(313)
stairs(sumbin)
% plot(bin_allpf)
hold on
line([edgobj(1, 1) edgobj(1, 1)], [0 max(sumbin)], 'color', 'm', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [0 max(sumbin)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [0 max(sumbin)], 'color', 'm', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [0 max(sumbin)], 'color', 'm', 'linestyle', '--');
ylim=[0,1]
hold off
% 
% subplot(312)
% v = nanmean(fct_matnorm([A{2}]), 1);
% plot(v)
% hold on
% line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
% line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
% line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
% line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'm', 'linestyle', '--');
% ylim=[0,1]
% hold off

% subplot(326)
% plot(sumbin_b)
% hold on
% line([edgobj(1, 1) edgobj(1, 1)], [min(sumbin_b) max(sumbin_b)], 'color', 'm', 'linestyle', '--');
% line([edgobj(1, 2) edgobj(1, 2)], [min(sumbin_b) max(sumbin_b)], 'color', 'm', 'linestyle', '--');
% line([edgobj(2, 1) edgobj(2, 1)], [min(sumbin_b) max(sumbin_b)], 'color', 'm', 'linestyle', '--');
% line([edgobj(2, 2) edgobj(2, 2)], [min(sumbin_b) max(sumbin_b)], 'color', 'm', 'linestyle', '--');
% ylim=[0,1]
end