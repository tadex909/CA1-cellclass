%select_your data
fpath = 'N:\Vinca\MATLAB\Data\All_data';
main_folder='N:\Vinca\MATLAB\Data\'
h=readtable(strcat(main_folder,'select_sess.xlsx'))
keep_dat=h.session_name(logical(h.doselect))
cd=fpath
%load data
load(strcat(fpath, '\supool_v3.mat'));

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
boxplot=1;
rankplot=1;
%Faire PO vs PNO plus tard quand supoolé plus de data. 
%Faire boucle avec 1 de g1à4 et j de g2à g5 pour toutes les comparaisons
%2à2
%% BOXPLOT

%Velocity
vel1 = supool.vel_cond(g1);
vel2 = supool.vel_cond(g2)
vel3 = supool.vel_cond(g3);
vel4 = supool.vel_cond(g4);
vel5 = supool.vel_cond(g5);
vel6 = supool.vel_cond(g6);

%Frequebce reward per sessions %REFLECHIR SUR CE DATASET
fr_rw1=supool.freq_rewardsess([g1*2; g1*2+1])
fr_rw2=supool.freq_rewardsess([g2*2 ;g2*2+1])
fr_rw3=supool.freq_rewardsess([g3*2; g3*2+1])
fr_rw4=supool.freq_rewardsess([g4*2; g4*2+1])
fr_rw5=supool.freq_rewardsess([g5*2; g5*2+1])
fr_rw6=supool.freq_rewardsess([g6*2; g6*2+1])

% pourcentage active pyramidal cells
ac1 = supool.pct_pyr_active_c(g1);
ac2 = supool.pct_pyr_active_c(g2);
ac3 = supool.pct_pyr_active_c(g3);
ac4 = supool.pct_pyr_active_c(g4);
ac5 = supool.pct_pyr_active_c(g5);
ac6 = supool.pct_pyr_active_c(g6);

%pourcentage place cells
pct_pf1 = supool.pct_pyr_sm_c(g1);
pct_pf2 = supool.pct_pyr_sm_c(g2);
pct_pf3 = supool.pct_pyr_sm_c(g3);
pct_pf4 = supool.pct_pyr_sm_c(g4);
pct_pf5 = supool.pct_pyr_sm_c(g5);
pct_pf6 = supool.pct_pyr_sm_c(g6);
% fprm.ylab = 'Bidirectional Cells (%)';
% v1 = supool.pct_pyr_bid_c(:, 1);
% fprm.ylab = 'Unidirectional Cells (%)';
% v1 = supool.pct_pyr_uni_c(:, 1);

     cell1=[]
    for k =1:length(g1);
    idi=find(supool.keptsession_idcel==g1(k));
    cell1=sort([cell1;idi])
    end
fr1 = supool.fr_uc(find(cell1));
'Out/In Field Firing Rate';
outovin_ucf = supool.outovin_ucf(:, :, 1);
o1 = outovin_ucf(find(cell1));
'PF Width (cm)';
size_ucf = 2*supool.size_ucf(:, :, 1);
s1 = size_ucf(find(cell1));
'PF Deviation (cm)';
lap_devwithmean_ucf = 2*supool.lap_devwithmean_ucf(:, :, 1);
l1 = lap_devwithmean_ucf(find(cell1));
'Spatial Information (bit/spk)';
si1 = supool.si_uc(find(cell1));
'Stability';
st1 = supool.stb_allpairs_uc(find(cell1));

'place fields distribution';
allrmap.feat.ispyr_uc = kron(supool.pyr_uc, ones(1, 2));
ispyrpf = supool.ispf_ucf(:, :, 1) & allrmap.feat.ispyr_uc;
fr_s_cxutemp=supool.fr_s_cxu(1:2,:,cell1);
ispyrpftemp=ispyrpf(cell1,1:2) 
nb_cond =2
nb_cel = length(supool.session_idcel(cell1));
A1 = cell(nb_cond, 1);
B1 = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
             if ispyrpftemp(u, c)  & c==1%si c une place cell 
        A1{c} = [A1{c} ; fr_s_cxutemp(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
         elseif ispyrpftemp(u, c)  & c==2%             A{c-1} = [A{c-1} ; fr_s_cxutemp(c-1, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
           A1{c} = [A1{c} ; supool.fr_s_cxu(c, : , u)]; 
             end
    end
end

     cell2=[]
    for k =1:length(g2);
    idi=find(supool.keptsession_idcel==g2(k));
    cell2=sort([cell2;idi])
    end
'FR';
fr2 = supool.fr_uc(find(cell2));
'Out/In Field Firing Rate';
o2 = outovin_ucf(find(cell2));
'PF Width (cm)';
s2 = size_ucf(find(cell2));
'PF Deviation (cm)';
l2 = lap_devwithmean_ucf(find(cell2));
'Spatial Information (bit/spk)';
si2 = supool.si_uc(find(cell2));
'Stability';
st2 = supool.stb_allpairs_uc(find(cell2));

'place fields distribution';
fr_s_cxutemp=supool.fr_s_cxu(1:2,:,cell2);
ispyrpftemp=ispyrpf(cell2,1:2) 
nb_cel = length(supool.session_idcel(cell2));
A2 = cell(nb_cond, 1);
B2 = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
             if ispyrpftemp(u, c)  & c==1%si c une place cell 
        A2{c} = [A2{c} ; fr_s_cxutemp(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
         elseif ispyrpftemp(u, c)  & c==2%             A{c-1} = [A{c-1} ; fr_s_cxutemp(c-1, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
           A2{c} = [A2{c} ; supool.fr_s_cxu(c, : , u)]; 
             end
    end
end

     cell3=[]
    for k =1:length(g3);
    idi=find(supool.keptsession_idcel==g3(k));
    cell3=sort([cell3;idi])
    end
'FR';
fr3 = supool.fr_uc(find(cell3));
'Out/In Field Firing Rate';
o3 = outovin_ucf(find(cell3));
'PF Width (cm)';
s3 = size_ucf(find(cell3));
'PF Deviation (cm)';
l3 = lap_devwithmean_ucf(find(cell3));
'Spatial Information (bit/spk)';
si3 = supool.si_uc(find(cell3));
'Stability';
st3 = supool.stb_allpairs_uc(find(cell3));

'place fields distribution';
fr_s_cxutemp=supool.fr_s_cxu(1:2,:,cell3);
ispyrpftemp=ispyrpf(cell3,1:2) 
nb_cel = length(supool.session_idcel(cell3));
A3 = cell(nb_cond, 1);
B3 = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
             if ispyrpftemp(u, c)  & c==1%si c une place cell 
        A3{c} = [A3{c} ; fr_s_cxutemp(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
         elseif ispyrpftemp(u, c)  & c==2%             A{c-1} = [A{c-1} ; fr_s_cxutemp(c-1, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
           A3{c} = [A3{c} ; supool.fr_s_cxu(c, : , u)]; 
             end
    end
end

     cell4=[]
    for k =1:length(g4);
    idi=find(supool.keptsession_idcel==g4(k));
    cell4=sort([cell4;idi])
    end
'FR';
fr4 = supool.fr_uc(find(cell4));
'Out/In Field Firing Rate';
o4 = outovin_ucf(find(cell4));
'PF Width (cm)';
s4 = size_ucf(find(cell4));
'PF Deviation (cm)';
l4 = lap_devwithmean_ucf(find(cell4));
'Spatial Information (bit/spk)';
si4 = supool.si_uc(find(cell4));
'Stability';
st4 = supool.stb_allpairs_uc(find(cell4));

'place fields distribution';
fr_s_cxutemp=supool.fr_s_cxu(1:2,:,cell4);
ispyrpftemp=ispyrpf(cell4,1:2) 
nb_cel = length(supool.session_idcel(cell4));
A4 = cell(nb_cond, 1);
B4 = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
             if ispyrpftemp(u, c)  & c==1%si c une place cell 
        A4{c} = [A4{c} ; fr_s_cxutemp(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
         elseif ispyrpftemp(u, c)  & c==2%             A{c-1} = [A{c-1} ; fr_s_cxutemp(c-1, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
           A4{c} = [A4{c} ; supool.fr_s_cxu(c, : , u)]; 
             end
    end
end

     cell5=[]
    for k =1:length(g5);
    idi=find(supool.keptsession_idcel==g5(k));
    cell5=sort([cell5;idi])
    end
'FR';
fr5 = supool.fr_uc(find(cell5));
'Out/In Field Firing Rate';
o5 = outovin_ucf(find(cell5));
'PF Width (cm)';
s5 = size_ucf(find(cell5));
'PF Deviation (cm)';
l5 = lap_devwithmean_ucf(find(cell5));
'Spatial Information (bit/spk)';
si5 = supool.si_uc(find(cell5));
'Stability';
st5 = supool.stb_allpairs_uc(find(cell5));

'place fields distribution';
fr_s_cxutemp=supool.fr_s_cxu(1:2,:,cell5);
ispyrpftemp=ispyrpf(cell5,1:2) 
nb_cel = length(supool.session_idcel(cell5));
A5 = cell(nb_cond, 1);
B5 = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
             if ispyrpftemp(u, c)  & c==1%si c une place cell 
        A5{c} = [A5{c} ; fr_s_cxutemp(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
         elseif ispyrpftemp(u, c)  & c==2%             A{c-1} = [A{c-1} ; fr_s_cxutemp(c-1, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
           A5{c} = [A5{c} ; supool.fr_s_cxu(c, : , u)]; 
             end
    end
end

     cell6=[]
    for k =1:length(g6);
    idi=find(supool.keptsession_idcel==g6(k));
    cell6=sort([cell6;idi])
    end
'FR';
fr6 = supool.fr_uc(find(cell6));
'Out/In Field Firing Rate';
o6 = outovin_ucf(find(cell6));
'PF Width (cm)';
s6 = size_ucf(find(cell6));
'PF Deviation (cm)';
l6 = lap_devwithmean_ucf(find(cell6));
'Spatial Information (bit/spk)';
si6 = supool.si_uc(find(cell6));
'Stability';
st6 = supool.stb_allpairs_uc(find(cell6));

'place fields distribution';
fr_s_cxutemp=supool.fr_s_cxu(1:2,:,cell6);
ispyrpftemp=ispyrpf(cell6,1:2) 
nb_cel = length(supool.session_idcel(cell6));
A6 = cell(nb_cond, 1);
B6 = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
             if ispyrpftemp(u, c)  & c==1%si c une place cell 
        A6{c} = [A6{c} ; fr_s_cxutemp(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
         elseif ispyrpftemp(u, c)  & c==2%             A{c-1} = [A{c-1} ; fr_s_cxutemp(c-1, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
           A6{c} = [A6{c} ; supool.fr_s_cxu(c, : , u)]; 
             end
    end
end

clear a15 a16 a17 a18 a19 a20 a21 a22 a23 a24 a25 allrmap ans boxplot c fpath fr_s_cxu_temp h idi ispyrpf ispyrpftemp keep_dat outovin_ucf lap_devwithmean_ucf main_folder nb_cel nb_cond rankplot size_ucf supool u fr_s-cxutemp

save('for_fig_fedwa')