main_folder='N:\Vinca\MATLAB\Data\';
h=readtable(strcat(main_folder,'select_sess.xlsx'));
keep_dat=h.session_name(logical(h.doselect));
age=h.Age(logical(h.doselect));
clear h;
outputfolder=strcat(main_folder,'All_data\');
cd=outputfolder;

for g=1:length(keep_dat)

load(strcat(string(keep_dat(g)),'_TrajData.mat'))
load(strcat(string(keep_dat(g)),'_Phenosys.mat'))
load(strcat(string(keep_dat(g)),'_Ratemap.mat'))
%% allcel 
clear allcel2
%allcel
allcel2.time_spk=double(allcel.time_spk2)/25000;
allcel2.id_spk=double(allcel.id_spk);
allcel2.itime_spk=double(allcel.itime_spk);
allcel2.id_cel=allcel.id_cel;
allcel2.samples_chnl_pts_spk_cel=allcel.waveform;
 allcel2.samples_chnl_pts_cel=allcel.meanwaveform;
 allcel2.nb_cel=allcel.nb_cel;
 allcel2.idtime_spk=[];
allcel2.state_cel_condw=allcel.state_uc(1:size(allrmap.feat.fr_uc,1),1:size(allrmap.feat.fr_uc,2));
 allcel2.fr_nr_cel=allcel.fr_u;
 
 %allrmap.feat
allcel2.fr_nr_cel_condw=allrmap.feat.fr_uc;
allcel2.stb_nr_cel_condw=allrmap.feat.stb_oddeven_uc;
allcel2.stb2_nr_cel_condw=allrmap.feat.stb_meancorr_uc;
allcel2.stb3_nr_cel_condw=allrmap.feat.stb_allpairs_uc;
allcel2.sparsity_nr_cel_condw=allrmap.feat.sparsity_uc;
 allcel2.si_mean_nr_cel_condw=allrmap.feat.si_uc;
 allcel2.mean_si_nr_cel_condw=allrmap.feat.si_meanlap_uc;
 
 %allpf.feat?/allrmap.feat?
 allcel2.pf_pctpfstb_nr_cel_condw=[];
 allcel2.pf_nb_nr_cel_condw=[];
 allcel2.pf_ispf_nr_cel_condw=[];
 allcel2.pf_nos_nr_cel_condw=[];
 allcel2.pf_noscorr_nr_cel_condw=[];
 allcel2.pf_size_nr_cel_condw=[];
 allcel2.pf_height_nr_cel_condw=[];
 allcel2.pf_icent_nr_cel_condw=[];
 allcel2.pf_spread_nr_cel_condw=[];
 allcel2.pf_maxdev_nr_cel_condw=[];
 allcel2.pf_spreadm_nr_cel_condw=[];
 allcel2.pf_meansize_nr_cel_condw=[];
 allcel2.pf_stdsize_nr_cel_condw=[];
 allcel2.pf_diffsize_nr_cel_condw=[];
 allcel2.pf_pctpf_nr_cel_condw=[];
 allcel2.pf_nr_cel_condw=[];
 
 for i=(1:2:size(allcel2.fr_nr_cel_condw,2))
     j=i+1;
      allcel2.state_cel_cond(:,j/2)=round(mean(allcel.state_uc(:,i:j),2));
      allcel2.state_cel_cond(: , j/2)=round(mean(allcel.state_uc(:,i:j),2));
       allcel2.fr_nr_cel_cond(: , j/2)=mean(allrmap.feat.fr_uc(:,i:j),2);
       allcel2.stb_nr_cel_cond(: , j/2)=mean(allrmap.feat.stb_oddeven_uc(:,i:j),2);
     allcel2.stb_nr_cel_cond(: , j/2)=mean(allrmap.feat.stb_oddeven_uc(:,i:j),2);
     allcel2.stb2_nr_cel_cond(: , j/2)=mean(allrmap.feat.stb_meancorr_uc(:,i:j),2);
     allcel2.stb3_nr_cel_cond(: , j/2)=mean(allrmap.feat.stb_allpairs_uc(:,i:j),2);
     allcel2.sparsity_nr_cel_cond(: , j/2)=mean(allrmap.feat.sparsity_uc(:,i:j),2);
     allcel2.si_mean_nr_cel_cond(: , j/2)=mean(allrmap.feat.si_uc(:,i:j),2);
     allcel2.mean_si_nr_cel_cond(: , j/2)=mean(allrmap.feat.si_meanlap_uc(:,i:j),2);
 end
 

 allcel2.pf_pctpfstb_nr_cel_cond=[];
 allcel2.pf_nb_nr_cel_cond=[];
 allcel2.pf_ispf_nr_cel_cond=[];
 allcel2.pf_nos_nr_cel_cond=[];
 allcel2.pf_noscorr_nr_cel_cond=[];
 allcel2.pf_size_nr_cel_cond=[];
 allcel2.pf_height_nr_cel_cond=[];
 allcel2.pf_icent_nr_cel_cond=[];
 allcel2.pf_spread_nr_cel_cond=[];
 allcel2.pf_maxdev_nr_cel_cond=[];
 allcel2.pf_spreadm_nr_cel_cond=[];
 allcel2.pf_meansize_nr_cel_cond=[];
 allcel2.pf_stdsize_nr_cel_cond=[];
 allcel2.pf_diffsize_nr_cel_cond=[];
 allcel2.pf_pctpf_nr_cel_cond=[];
 allcel2.pf_nr_cel_cond=[];
 
 %allrmap.feat
 allcel2.scd_nr_cel_cond=allrmap.feat.scd_nr_uc(1:size(allcel2.state_cel_cond,1),1:size(allcel2.state_cel_cond,2));
 allcel2.scp_nr_cel_cond=allrmap.feat.scp_nr_uc(1:size(allcel2.state_cel_cond,1),1:size(allcel2.state_cel_cond,2));
 allcel2.ipp_nr_cel_cond=allrmap.feat.inddis_nr_uc(1:size(allcel2.state_cel_cond,1),1:size(allcel2.state_cel_cond,2));
 allcel2.ipd_nr_cel_cond=allrmap.feat.indpos_nr_uc(1:size(allcel2.state_cel_cond,1),1:size(allcel2.state_cel_cond,2));
 
 %allcel
 allcel2.acg10_bin_cel=allcel.acg10_bin_u;
 allcel2.acg1_bin_cel=allcel.acg1_bin_u;
 allcel2.ref_per_cel=allcel.ref_per_u;
 allcel2.burst_cel=allcel.burst_u;
 allcel2.duration_cel=allcel.duration_u;
 allcel2.asymmetry_cel=allcel.asymmetry_u;
 allcel2.peaktrough_cel=allcel.peaktrough_u;
 allcel2.type_cel=allcel.type_u;
 
 allcel2.pyr_cel_cond=allcel.pyr_uc(1:size(allcel2.state_cel_cond,1),1:size(allcel2.state_cel_cond,2));
 allcel2.pyr_active_cel_cond=allcel.pyr_active_uc(1:size(allcel2.state_cel_cond,1),1:size(allcel2.state_cel_cond,2));
 allcel2.pyr_inactive_cel_cond=allcel.pyr_inactive_uc(1:size(allcel2.state_cel_cond,1),1:size(allcel2.state_cel_cond,2));
 allcel2.pyr_sm_cel_cond=allcel.pyr_sm_uc(1:size(allcel2.state_cel_cond,1),1:size(allcel2.state_cel_cond,2));
 allcel2.pyr_bid_cel_cond=allcel.pyr_bid_uc(1:size(allcel2.state_cel_cond,1),1:size(allcel2.state_cel_cond,2));
 allcel2.pyr_uni_cel_cond=allcel.pyr_uni_uc(1:size(allcel2.state_cel_cond,1),1:size(allcel2.state_cel_cond,2));
 allcel2.pyr_nonsm_cel_cond=allcel.pyr_nonsm_uc(1:size(allcel2.state_cel_cond,1),1:size(allcel2.state_cel_cond,2));
 
 %allcel
 allcel2.nb_pyr_cond= allcel.nb_pyr_c(~isnan(allcel.pct_pyr_active_c));
 allcel2.nb_pyr_active_cond=allcel.nb_pyr_active_c(~isnan(allcel.pct_pyr_active_c));
 allcel2.nb_pyr_inactive_cond=allcel.nb_pyr_inactive_c(~isnan(allcel.pct_pyr_active_c));
 allcel2.pct_pyr_active_cond=allcel.pct_pyr_active_c(~isnan(allcel.pct_pyr_active_c));
 allcel2.pct_pyr_inactive_cond=allcel.pct_pyr_inactive_c(~isnan(allcel.pct_pyr_inactive_c));
 allcel2.nb_pyr_sm_cond=allcel.nb_pyr_sm_c(~isnan(allcel.pct_pyr_active_c)); 
 allcel2.nb_pyr_bid_cond=allcel.nb_pyr_bid_c(~isnan(allcel.pct_pyr_active_c));
 allcel2.nb_pyr_uni_cond=allcel.nb_pyr_uni_c(~isnan(allcel.pct_pyr_active_c));
 allcel2.nb_pyr_nonsm_cond=allcel.nb_pyr_nonsm_c(~isnan(allcel.pct_pyr_active_c));
 allcel2.pct_pyr_sm_cond=allcel.pct_pyr_sm_c(~isnan(allcel.pct_pyr_active_c));
 allcel2.pct_pyr_bid_cond=allcel.pct_pyr_bid_c(~isnan(allcel.pct_pyr_active_c));
 allcel2.pct_pyr_uni_cond=allcel.pct_pyr_uni_c(~isnan(allcel.pct_pyr_active_c));
 allcel2.pct_pyr_nonsm_cond=allcel.pct_pyr_nonsm_c(~isnan(allcel.pct_pyr_active_c));

 %% cel 
cel(j).frmap_nr.rate_s_tr_bin(j)=1
cel(j).frmap_nr.rate_s_tr_bin(j)=[]
 
 j=1;
for i=[allcel.id_cel]'
        
        spk_cel=find(allcel.id_spk==i);
        cel(j).idtime_spk=double(allcel.itime_spk(spk_cel));
        cel(j).frmap_nr.rate_s_tr_bin(j)=1
        cel(j).frmap_nr.rate_s_tr_bin(j)=[]
        cel(j).frmap_nr.rate_s_tr_bin=allrmap.fr_txu(:,:,j);
        cel(j).frmap_nr.rate_s_condw_bin=(allrmap.fr_cxu(:,:,j))
        j=j+1;
end
% cel.frmap_nr=[];
% cel.pssf_condw_bound=[]
% cel.pf_nr_tr_bin=[]
% cel.pf_nr_condw_bin=[]
% cel.one_pf_nr_condw_bin=[]
% cel.one_pf_nr_tr_bin=[]
 %% bhv
 
bhv.time_d=t_ds;
bhv.p_x_d = X_ds;
bhv.p_y_d = Y_ds;
bhv.p_x_ds = X_ds_n.';
bhv.v_x_d = XSpeed.';
bhv.v_x_ds =XSpeed.';
bhv.sdata =size(t_ds,2);
bhv.way_tr =strcmp({Traj.LR},'W').' ;
bhv.idtrack_tr(:,1) =[Traj(1:end).start].';
bhv.idtrack_tr(:,2) =[Traj(1:end).stop].';
bhv.idreward_tr=[Traj(1:end).stop].';
bhv.nb_tracks =size(Traj,2);
bhv.nb_rewards= size(Traj,2);
bhv.v_x_ds2 =[];
bhv.ts_imm_tr=arrayfun(@(x)(sum(x.StopDurations)),Traj);
bhv.vel_tr =arrayfun(@(x)(sum(x.XSpeed)),Traj);
bhv.ts_objr_tr= [];
bhv.freq_reward_cond= [];
bhv.ts_imm_cond =[];
bhv.vel_cond =[];
bhv.ts_objr_cond =[];
bhv.freq_reward_condf= [];
bhv.ts_imm_condf =[];
bhv.vel_condf =[];
bhv.ts_objr_condf =[];
k=1;
condw=zeros(1,length(unique([Traj.Cond]))*2)
for n=[1:length(unique([Traj.Cond]))];
    n
    for m=[1,2]
        m
    WB=['W','B'];
    WB=WB(m);
       for i=[1:size(Traj,2)]
           i
       a=Traj(i).Cond==n&Traj(i).LR==WB
       condw(k)=condw(k)+a;
       
       end
       k=k+1;
    end
end
bhv.nb_tracks_condw =condw;
bhv.nb_tracks_condfw =[];
bhv.ts_rev_tr_bin =[];
bhv.vel_rev_tr_bin =[];
bhv.ts_rev_cond_bin =[];
bhv.vel_rev_cond_bin =[];
bhv.ts_rev_condf_bin =[];
bhv.vel_rev_condf_bin =[];

 %% eprm.icond 
 eprm.icond=unique([Traj.Cond]) %PAS VRAIMENT JUSTE : devrait correspondre à la condition codifiée genre PO=1 et PNO =2 et POnM=3 POnMb=4 PONOTREE=5 etc. A faire.
 %% end things
 allcel = allcel2;
clear allcel2;

save(strcat(  outputfolder , string(keep_dat(g)) ,'_CellClass.mat'),'bhv','allcel','cel','eprm');
end
 