
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



v1 =sum( supool.nb_pyr_sm_c(g1));
v2 = sum(supool.nb_pyr_sm_c(g2))
v3=sum(supool.nb_pyr_sm_c(g3))
v4 =sum( supool.nb_pyr_sm_c(g4));
v5 = sum(supool.nb_pyr_sm_c(g5))
v6=sum(supool.nb_pyr_sm_c(g6))

X = [v1; v2; v3;v4;v5;v6];
grp = [ones(size(v1)); 2.*ones(size(v2)); 3.*ones(size(v3));4.*ones(size(v4));5.*ones(size(v5));6.*ones(size(v6))];
boxplot(X, grp,'boxstyle','filled')
title('nb pc PO')