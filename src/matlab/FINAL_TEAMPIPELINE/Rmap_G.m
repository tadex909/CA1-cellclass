%%%%You need to load _TrajData.mat / _Phenosys.mat / _ePhy.mat %%%%%%%%%


loadcolmathSFN
load colormaprom.mat

savepath = OutputFolder;

filesep = '/';
datalc.session_name = Namesession;


%%%%You need to use the Traj structure of PP

maxmaze = max(X_ds_n); %X_ds?
X_ds_n = X_ds_n*100/maxmaze;
maxmaze = max(X_ds_n);

%%%index of condition and trial  %OLD ONLY WORKS FOR 2Conds
% A=size(Traj);
% cond=ones(A(2),1);
% Cond2 = find([Traj(:).Cond]>1);
% cond(Cond2)=2;
% way=strcmp({Traj.LR}, 'RL');
% way=+way;
% way(Cond2) = way(Cond2)+ones(length(Cond2),1)';
% condway = cond+way';

%% index of condition and trial
A=size(Traj);
condway=ones(A(2),1);
% Cond2 = find([Traj(:).Cond]==1);
% cond(Cond2)=2;
way=(strcmp({Traj.WB},'W' )).';
back=(strcmp({Traj.WB}, 'B')).';
for i=1:max([Traj(:).Cond])
    maskingall=([Traj(:).Cond]==i).';
    maskingodd=find(back==maskingall & back);
    %cond(maskingall)=way(maskingall)+i;
    
    condway(maskingodd)=back(maskingodd)+i;
    if i > 1 
        maskingway=(way==maskingall & way);
        condway(maskingway)=2*i-1;
        condway(maskingodd)=2*i;
    end
% way(Cond2) = way(Cond2)+ones(length(Cond2),1)';
% condway = cond+way';
end

%% index of trial (for geoffrey it's 100Hz, for us it's 1000Hz)%%%%

indextr = [];
for tt=1:A(2)
    tp = [Traj(tt).start Traj(tt).stop];
    indextr = [indextr;tp];
end



bhv.prm.freq = 25000;

prm.rmap.idcond_t = condway;
prm.rmap.tbin = indextr;
prm.rmap.idx_rem = []; %%%let like that for now
prm.rmap.xbin = 0:floor(maxmaze/100):maxmaze;
prm.rmap.xbin_rem = 10; %%%normally it's 10
%prm.rmap.nb_cond = 10; OLDJULIE
nbsess=length(unique([Traj.Cond]));
prm.rmap.nb_cond = nbsess*2;
prm.rmap.freq = 1000;
prm.rmap.ismooth = 10; 

prm.rmap.frmax = 1; %1.5 %1  %%%%% 
prm.rmap.frmean = 0.3; %0.5 %0.3
prm.rmap.pctspklap = 0.50;  %%% ancien 0.50

prm.pf.mtd = 'random'; %random %poisson %circular_shift

prm.pf.nb_rep = 1000; %500 %1000

prm.pf.min_len = 3;
prm.pf.max_len = 45; %Max length %OLDrom = 45

prm.pf.max_btw_pf = 3; %OLD 5 %max dist between 2 pour merged
prm.pf.max_ext_pf = 5; %oldrom=5 %etendre pf  si pval=0.30

prm.pf.pval_crit = 0.01; % OLD 0.01 pval max pour sig
prm.pf.pval_edge_pf = 0.30; %OLD 0.30 augmentation taille pf pour vals adjacentes inf � 

%corerelation spatiale entre pf
prm.pf.pctstab = 0.4; %%%% 40%des essais %OLDROM=0.4
prm.pf.stab = 0.6; %%% 60% min

prm.pf.frmax = prm.rmap.frmax;  %%
prm.pf.frmean = prm.rmap.frmean;
prm.pf.pctspklap = prm.rmap.pctspklap;

%cvlassification pyr pour pyrpf 
prm.cellclass.burstind = 0;
prm.cellclass.frmaxinterneuron = 20; %20Hz max pour pyr
prm.cellclass.duration = 0.4;

maxc = 5;


%%%%IMPORTANT%%%%%%
%%%%NB of condition%%%%%%
eprm.nb_cond = nbsess;
%.nb_cond = nbsess;

NbCell = length(allcel.id_cel);

for g = 1:NbCell
    g;
    idx_spk = allcel.itime_spk(allcel.id_spk == allcel.id_cel(g));
    idx_spk = fct_ifreq_swap(idx_spk,bhv.prm.freq,1000);
    [rmap(g),pf(g)] = fct_mapandfields(X_ds_n,idx_spk,prm);
    
end


%%%%to save the plots
%
%path_save=('X:\julie\DATA\MM1\data processed\MM1_2024-03-19_15-05-51')
fct_plot_rmap(rmap,pf,'save',AnlPath)
%}
%%%%%%%%%


[allrmap, allpf] = fct_mapandfields_reshape(rmap, pf);  %D�tection PF

allcel.nb_cel = NbCell;

%% Codage Distance/Position

allrmap.feat.scd_nr_uc = NaN(allcel.nb_cel, maxc);
allrmap.feat.scp_nr_uc = NaN(allcel.nb_cel, maxc);

allrmap.feat.inddis_nr_uc = NaN(allcel.nb_cel, maxc);
allrmap.feat.indpos_nr_uc = NaN(allcel.nb_cel, maxc);



for g = 1:allcel.nb_cel
    for c = 1:eprm.nb_cond
        ind_c = (2*(c-1) + 1):(2*(c-1) + 2);
        [allrmap.feat.scd_nr_uc(g, c), allrmap.feat.scp_nr_uc(g, c)] = fct_spatial_corr(rmap(g).fr_s_cx(ind_c(1), :), rmap(g).fr_s_cx(ind_c(2), :));
        [allrmap.feat.inddis_nr_uc(g, c), allrmap.feat.indpos_nr_uc(g, c)] = fct_distposindexpeaktopeak(squeeze(allpf.feat.ifrmax_ucf(g, ind_c(1), :)), squeeze(allpf.feat.ifrmax_ucf(g, ind_c(2), :))); 
    end
end

%% Spatial correlation avant/apr�s muscimol : on corr�le allers sans muscimol avec allers sous muscimol, de m�me pour les retours 
% allrmap.feat.scp_aftmusci_uc = sc allers/allers_musci ; sc retours/retours_musci; sc allers_musci/allers_musci2 ; sc retours_musci/retours_musci2
% allrmap.feat.scd_aftmusci_uc : on s'en fout (codage de distance pas pertinent dans cette config)

allrmap.feat.scd_aftmusci_uc = NaN(allcel.nb_cel, maxc);
allrmap.feat.scp_aftmusci_uc = NaN(allcel.nb_cel, maxc);

for g = 1:allcel.nb_cel
    for c = 1:(2*(eprm.nb_cond - 2) + 2) 
        ind_c = [((c-1) + 1),((c-1) + 3)];
        [allrmap.feat.scd_aftmusci_uc(g, c), allrmap.feat.scp_aftmusci_uc(g, c)] = fct_spatial_corr(rmap(g).fr_s_cx(ind_c(1), :), rmap(g).fr_s_cx(ind_c(2), :));
    end
end
%% AutoCorrelogram
%allcel.itime_spk = allcel.itime_spk./30000;


 [allcel.acg10_bin_u, allcel.acg10_window] = ACG(allcel.itime_spk, allcel.id_spk, allcel.id_cel, bhv.prm.freq, 10, 700);
 [allcel.acg1_bin_u, allcel.acg1_window] = ACG(allcel.itime_spk, allcel.id_spk, allcel.id_cel, bhv.prm.freq, 1, 50);

%% Refractory Period, Burst Index & Mean Firing Rate

for g = 1:allcel.nb_cel
    allcel.ref_per_u(g) = ACGrefractoryT_new(allcel.acg1_bin_u(:, g));
    allcel.burst_u(g) = fct_burst_index(allcel.ref_per_u(g), allcel.acg1_bin_u(:, g));
    allcel.fr_u(g) = nanmean(rmap(g).fr_s_tx(:));
end


% Spike Features

% %%%find waveform 
% A = size(allcel.meanwaveform);
% allwave = [];
% for ii=1:NbCell
%     allmin=[];
%     for tt=1:A(2) %timewindow/sample
%         tp = allcel.meanwaveform(:,tt,ii); 
%         tpmin = min(tp);
%         allmin=[allmin,tpmin];
%     end
%     [good,ind]=min(allmin);
%     goodwave = allcel.meanwaveform(:,ind,ii);
%     allwave = [allwave;goodwave]; %les spikes min de chacun des 23 timestamps des n cellules 
% end


%[allcel.duration_u, allcel.asymmetry_u, allcel.peaktrough_u] = WaveformFeature_SC(allwave, bhv.prm.freq);  


[allcel.duration_u, allcel.asymmetry_u, allcel.peaktrough_u] = WaveformFeatureVinca(allcel.waveform, allcel.meanwaveform, allcel.bestwaveform,allcel.bestswaveforms,bhv.prm.freq);  

%% Type de Cellule (Pyramidale / Interneurone)

allcel.type_u = ones(length(allcel.id_cel), 1);

ind_tmp1 = find(nansum(squeeze(nanmean(allrmap.fr_s_cxu, 2)) > prm.cellclass.frmaxinterneuron, 1) > 0);
ind_tmp2 = find((allcel.burst_u < prm.cellclass.burstind & allcel.duration_u*1000 < prm.cellclass.duration));
ind_tmp = [ind_tmp1 ind_tmp2];
ind_tmp = unique(ind_tmp);

if ~isempty(ind_tmp)
    allcel.type_u(ind_tmp) = 0;   
end

%% Active/Silent

allcel.state_uc = NaN(allcel.nb_cel, maxc*2);
bhv.icondw_tr = condway;

for g = 1:allcel.nb_cel
    for c = 1:(eprm.nb_cond*2)
        is_c = bhv.icondw_tr == c;
        ratem =  rmap(g).fr_s_tx(is_c, :);  
        allcel.state_uc(g, c) = double(fct_isactive(ratem, prm.rmap));
    end
end


% D?termination des crit?res

crit_pyr_u = (allcel.type_u == 1);
crit_pyr_uc = zeros(allcel.nb_cel, maxc);
crit_pyr_uc(:, 1:eprm.nb_cond) = repmat(crit_pyr_u, 1, eprm.nb_cond);

crit_itn_u = (allcel.type_u == 0);                  %%%%
crit_itn_uc = zeros(allcel.nb_cel, maxc);
crit_itn_uc(:, 1:eprm.nb_cond) = repmat(crit_itn_u, 1, eprm.nb_cond);

crit_pf_stable_uc = allpf.feat.ispf_ucf(:, :, 1);

crit_active_uc = false(allcel.nb_cel, maxc);
crit_inactive_uc = false(allcel.nb_cel, maxc);
crit_bid_uc = false(allcel.nb_cel, maxc);
crit_uni_uc = false(allcel.nb_cel, maxc);
crit_nonsm_uc = false(allcel.nb_cel, maxc);

crit_itn_active_uc = false(allcel.nb_cel, maxc);    %%%%
crit_itn_inactive_uc = false(allcel.nb_cel, maxc);
crit_itn_bid_uc = false(allcel.nb_cel, maxc);
crit_itn_uni_uc = false(allcel.nb_cel, maxc);
crit_itn_nonsm_uc = false(allcel.nb_cel, maxc);
allrmap.feat.fr_c=[]; 

for c = 1:eprm.nb_cond
    ind_c = (2*(c-1) + 1):(2*(c-1) + 2);
    crit_active_uc(:, c) = (sum(allcel.state_uc(:, ind_c), 2) >= 1); % Active dans au moins un sens
    crit_inactive_uc(:, c) = (sum(allcel.state_uc(:, ind_c), 2) == 0);  % Non active dans les deux sens
    
    crit_bid_uc(:, c) = crit_pf_stable_uc(:, ind_c(1)) & crit_pf_stable_uc(:, ind_c(2));
    crit_uni_uc(:, c) = xor(crit_pf_stable_uc(:, ind_c(1)), crit_pf_stable_uc(:, ind_c(2)));
    crit_nonsm_uc(:, c) = ~crit_pf_stable_uc(:, ind_c(1)) & ~crit_pf_stable_uc(:, ind_c(2));
    fr_c =  nanmean(allrmap.feat.fr_uc ( : , c : c+1 ),2);
    allrmap.feat.fr_c= [allrmap.feat.fr_c fr_c];
end

% Recalcul nombre cellules -> diff�rencier cellules inactives des cellules que l'on perd (pyr_off/itn_off)

%JAI AUTOMATISE POUR QUE CA VARIE AVEC NB DE CONDITIONS VOIR AU DESSUS
% p = nanmean(allrmap.feat.fr_uc ( : , 1 : 2 ),2);
% q = nanmean(allrmap.feat.fr_uc ( : , 3 : 4 ),2);
% % r = nanmean(allrmap.feat.fr_uc ( : , 5 : 6 ),2);
% % s = nanmean(allrmap.feat.fr_uc ( : , 7 : 8 ),2);
% % t = nanmean(allrmap.feat.fr_uc ( : , 9 : 10 ),2);
% allrmap.feat.fr_c = [p q r s t];

lost_cell_c = zeros(size(allrmap.feat.fr_c));
lost_cell_c(find(allrmap.feat.fr_c == 0)) = 1;

crit_present_cell_uc = (lost_cell_c == 0);
crit_lost_cell_uc = (lost_cell_c == 1);
crit_present_cell_uc = [crit_present_cell_uc,zeros(NbCell,5-size(crit_present_cell_uc ,2))];
crit_lost_cell_uc=[crit_lost_cell_uc,zeros(NbCell,5-size(crit_lost_cell_uc ,2))];
%si n�cessaire : crit_present_cell_uc=[crit_present_cell_uc, zeros(30,3)]

% D?termination des cellules 
allcel.pyr_uc = crit_pyr_uc & crit_present_cell_uc;
allcel.pyr_active_uc = crit_pyr_uc & crit_present_cell_uc & crit_active_uc;
allcel.pyr_inactive_uc = crit_pyr_uc & crit_present_cell_uc & crit_inactive_uc;
allcel.pyr_off_uc = crit_pyr_uc & crit_lost_cell_uc;

allcel.pyr_sm_uc = crit_pyr_uc & crit_present_cell_uc & crit_active_uc & (crit_bid_uc | crit_uni_uc); % active is not necessary since bid and uni cells are based on a place field which is never determined for silent cell
allcel.pyr_bid_uc = crit_pyr_uc & crit_present_cell_uc & crit_active_uc & crit_bid_uc;
allcel.pyr_uni_uc = crit_pyr_uc & crit_present_cell_uc & crit_active_uc & crit_uni_uc;
allcel.pyr_nonsm_uc = crit_pyr_uc & crit_present_cell_uc & crit_active_uc & crit_nonsm_uc;

allcel.itn_uc = crit_itn_uc & crit_present_cell_uc;
allcel.itn_active_uc = crit_itn_uc & crit_present_cell_uc & crit_active_uc;
allcel.itn_inactive_uc = crit_itn_uc & crit_present_cell_uc & crit_inactive_uc;
allcel.itn_off_uc = crit_itn_uc & crit_lost_cell_uc;

allcel.itn_sm_uc = crit_itn_uc & crit_present_cell_uc & crit_active_uc & (crit_bid_uc | crit_uni_uc); % active is not necessary since bid and uni cells are based on a place field which is never determined for silent cell
allcel.itn_bid_uc = crit_itn_uc & crit_present_cell_uc & crit_active_uc & crit_bid_uc;
allcel.itn_uni_uc = crit_itn_uc & crit_present_cell_uc & crit_active_uc & crit_uni_uc;
allcel.itn_nonsm_uc = crit_itn_uc & crit_present_cell_uc & crit_active_uc & crit_nonsm_uc;

% Nombre de cellules et pourcentages 
allcel.nb_cel_uc = allcel.nb_cel*ones(allcel.nb_cel , 5);
allcel.nb_cel_on = mean(allcel.nb_cel_uc) - sum(allcel.pyr_off_uc + allcel.itn_off_uc); %allcel.nb_cel = nb cell enregistr�es sur tout l'enregistrement

allcel.nb_pyr_c = sum(allcel.pyr_uc); 
allcel.nb_pyr_active_c = sum(allcel.pyr_active_uc);
allcel.nb_pyr_inactive_c = sum(allcel.pyr_inactive_uc);
allcel.nb_pyr_off_c = sum(allcel.pyr_off_uc);
allcel.pct_pyr_c = (allcel.nb_pyr_c ./ allcel.nb_cel_on)*100; 
allcel.pct_pyr_active_c = (allcel.nb_pyr_active_c ./ allcel.nb_pyr_c)*100;
allcel.pct_pyr_inactive_c = (allcel.nb_pyr_inactive_c ./ allcel.nb_pyr_c)*100;

allcel.nb_pyr_sm_c = sum(allcel.pyr_sm_uc);
allcel.nb_pyr_bid_c = sum(allcel.pyr_bid_uc);
allcel.nb_pyr_uni_c = sum(allcel.pyr_uni_uc);
allcel.nb_pyr_nonsm_c = sum(allcel.pyr_nonsm_uc);
allcel.pct_pyr_sm_c = (allcel.nb_pyr_sm_c ./ allcel.nb_pyr_active_c)*100;
allcel.pct_pyr_bid_c = (allcel.nb_pyr_bid_c ./ allcel.nb_pyr_sm_c)*100;
allcel.pct_pyr_uni_c = (allcel.nb_pyr_uni_c ./ allcel.nb_pyr_sm_c)*100;
allcel.pct_pyr_nonsm_c = (allcel.nb_pyr_nonsm_c ./ allcel.nb_pyr_active_c)*100;

allcel.nb_itn_c = sum(allcel.itn_uc); 
allcel.nb_itn_active_c = sum(allcel.itn_active_uc);
allcel.nb_itn_inactive_c = sum(allcel.itn_inactive_uc);
allcel.nb_itn_off_c = sum(allcel.itn_off_uc);
allcel.pct_itn_c = (allcel.nb_itn_c ./ allcel.nb_cel_on)*100; 
allcel.pct_itn_active_c = (allcel.nb_itn_active_c ./ allcel.nb_itn_c)*100;
allcel.pct_itn_inactive_c = (allcel.nb_itn_inactive_c ./ allcel.nb_itn_c)*100;

allcel.nb_itn_sm_c = sum(allcel.itn_sm_uc);
allcel.nb_itn_bid_c = sum(allcel.itn_bid_uc);
allcel.nb_itn_uni_c = sum(allcel.itn_uni_uc);
allcel.nb_itn_nonsm_c = sum(allcel.itn_nonsm_uc);
allcel.pct_itn_sm_c = (allcel.nb_itn_sm_c ./ allcel.nb_itn_active_c)*100;
allcel.pct_itn_bid_c = (allcel.nb_itn_bid_c ./ allcel.nb_itn_sm_c)*100;
allcel.pct_itn_uni_c = (allcel.nb_itn_uni_c ./ allcel.nb_itn_sm_c)*100;
allcel.pct_itn_nonsm_c = (allcel.nb_itn_nonsm_c ./ allcel.nb_itn_active_c)*100;

allcel.prm = prm.cellclass;



save([savepath filesep datalc.session_name '_Ratemap'], 'allrmap', 'allpf', 'allcel', 'rmap', 'pf');


%yo=1




