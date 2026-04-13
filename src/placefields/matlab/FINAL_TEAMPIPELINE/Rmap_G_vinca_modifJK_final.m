%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%Cell classification / Ratemap / PC detection / AND MORE%%%%%%%%%%%%%%
%%%All Geoffrey's variables%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%The key difference in this code is the X you give to
%%%%%%'fct_mapandfields' function!!!!!!!!

%%%%%%Old version (2024) in X : all tracks aligned to the starting point
%%%%%%This one : X = all tracks aligned to the position (like old Geoff)

%% SEE THIS PART : change of the signal to get position instead of distance
asdsad
%%%%I also added the speed threshold for the rate maps, you didn't have it
%%%For me it works but I was at 100Hz but normally here it should work,
%%%just do a test first. 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Give a folder + session name(s), then files are discovered automatically
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

folder = '\\Epsztein-nas02\TEAM\Vinca\MATLAB\Data\Data_thèse_Vinca';   % dossier contenant les .mat de session
session_name = '';  % ex: 'MM1_2024-03-19_15-05-51' ; vide = toutes les sessions
run_only_rmap_pf = true;  % true: stop after ratemap + place-field outputs
overwrite_existing = false; % false: skip sessions with existing output file
skip_sessions = {'VS129_2025-11-25_16-51-07'}; % e.g. {'VS100_2024-11-03_14-40-11'}

if isempty(session_name)
    traj_files = [dir(fullfile(folder, '*_TrajData.mat')); dir(fullfile(folder, '*_trajData.mat'))];
    if isempty(traj_files)
        error('No *_TrajData.mat or *_trajData.mat files found in folder: %s', folder);
    end
    session_list = cell(numel(traj_files), 1);
    for ii = 1:numel(traj_files)
        session_list{ii} = regexprep(traj_files(ii).name, '(_TrajData|_trajData)\.mat$', '');
    end
    session_list = unique(session_list, 'stable');
else
    if ischar(session_name) || isstring(session_name)
        session_list = cellstr(session_name);
    elseif iscell(session_name)
        session_list = session_name;
    else
        error('session_name must be a char, string, or cell array of names.');
    end
end

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(fileparts(fileparts(script_dir))));
for k = 1:length(session_list)
    Namesession = strtrim(char(session_list{k}));

    if any(strcmpi(Namesession, skip_sessions))
        disp(['Skipping session (manual skip): ' Namesession])
        continue
    end

    sess_parts = strsplit(Namesession, '_');
    if ~isempty(sess_parts), mouse_name = sess_parts{1}; else, mouse_name = 'unknown_mouse'; end
    savepath = fullfile(repo_root, 'data', 'matlab', mouse_name, Namesession);
    if ~exist(savepath, 'dir'), mkdir(savepath); end

    if run_only_rmap_pf
        out_file = fullfile(savepath, [Namesession '_ratemap_circular.mat']);
    else
        out_file = fullfile(savepath, [Namesession '_Ratemap_final_CASE1.mat']);
    end
    if ~overwrite_existing && exist(out_file, 'file') == 2
        disp(['Skipping session (output exists): ' Namesession])
        continue
    end

    try
        traj_path = find_session_file(folder, Namesession, {'_TrajData.mat','_trajData.mat'});
        ephy_path = find_session_file(folder, Namesession, {'_ePhy.mat','_ePhys.mat','_ephys.mat','_EPhys.mat'});
        pheno_path = find_session_file(folder, Namesession, {'_Phenosys.mat','_phenosys.mat'});
    catch ME
        warning('Skipping session %s (missing input triplet): %s', Namesession, ME.message);
        continue
    end

    disp(['Session : ' Namesession])
    disp(['  TrajData : ' traj_path])
    disp(['  ePhys    : ' ephy_path])
    disp(['  Phenosys : ' pheno_path])

    try
        load(traj_path);
        load(ephy_path);
        load(pheno_path);
    catch ME
        warning('Skipping session %s (load failed): %s', Namesession, ME.message);
        continue
    end
%       load(fnames)
    % ... ton traitement ici
loadcolmathSFN
% load colormaprom.mat  % handled (optionally) by loadcolmathSFN

filesep = '/';
datalc.session_name = Namesession;

%%%%%Downsample freq%%%%%%
prm.freq_d=1000;


%%%%You need to use the Traj structure of PP
maxmaze = max(X_ds_n); 
X_ds_n = X_ds_n*100/maxmaze;
maxmaze = max(X_ds_n);

bhv.d_x_1000 = X_ds_n;


%{
%%%%%%%NOT USED NOW%%%%%%%%%%%
%%%%%%%I kept it anyway%%%%%%%
%SEED 1000Hz
Speed = diff(X_ds_n)./0.001;
Speed = [Speed(1);Speed];
Speed_s = fct_smoothgauss(Speed, 500);
%}


%% index of condition and trial
A=size(Traj);
condway=ones(A(2),1);
way=(strcmp({Traj.WB},'W' )).';
back=(strcmp({Traj.WB}, 'B')).';
bhv.way_tr=double(way);
for i=1:max([Traj(:).Cond])
    maskingall=([Traj(:).Cond]==i).';
    maskingodd=find(back==maskingall & back);
    condway(maskingodd)=back(maskingodd)+i;
    if i > 1 
        maskingway=(way==maskingall & way);
        condway(maskingway)=2*i-1;
        condway(maskingodd)=2*i;
    end
end


%% index of trial (for geoffrey it's 100Hz, for us it's 1000Hz)%%%%
indextr = [];
for tt=1:A(2)
    tp = [Traj(tt).start Traj(tt).stop];
    indextr = [indextr;tp];
end

bhv.itrack_tr_1000 = indextr;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%IMPORTANT TO DO YOUR SPEED THRESHOLD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Downsampling / Velocity / Smooth
prm.freq_d = 1000;
prm.idsmooth = (prm.freq_d / 2);
%% Velocity
bhv.d_x_1000s = bhv.d_x_1000;
bhv.v_x_1000_s2 = zeros(1, length(bhv.d_x_1000)- 1);
% for k = 1:length(indextr_d(:,1)) OLD JULIE MISTAKE?
 for k = 1:length(indextr(:,1)) 
    idx = bhv.itrack_tr_1000(k, 1):bhv.itrack_tr_1000(k, 2);
    bhv.d_x_1000s(idx) = fct_smoothgauss(bhv.d_x_1000(idx), prm.idsmooth); %%%smooth position
%     C = diff(bhv.d_x_100s(idx))./(1/prm.freq_d); %%%%Speed %OLD JULIE MISTAKE?
    C = diff(bhv.d_x_1000s(idx))./(1/prm.freq_d); %%%%Speed
    bhv.v_x_1000_s2(idx(1):idx(end-1)) = fct_smoothgauss(C, prm.idsmooth); %%%%smooth speed
end
bhv.v_x_1000_s2 = [bhv.v_x_1000_s2(1) bhv.v_x_1000_s2];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Change of the signal to get position instead of distance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Changing Distance to Position')
bhv.p_x_d = bhv.d_x_1000;
bhv.p_x_ds = bhv.d_x_1000s;

if bhv.way_tr(1) == 0
    for k = find(bhv.way_tr == 1)'
        idx = bhv.itrack_tr_1000(k, 1):bhv.itrack_tr_1000(k, 2);
            bhv.p_x_d(idx) = bhv.d_x_1000(bhv.itrack_tr_1000(k - 1, 2)) + bhv.d_x_1000(idx(1)) - bhv.d_x_1000(idx);
            bhv.p_x_ds(idx) = bhv.d_x_1000s(bhv.itrack_tr_1000(k - 1, 2)) + bhv.d_x_1000s(idx(1)) - bhv.d_x_1000s(idx);
    end

else

    for k = find(bhv.way_tr == 0)'
           idx = bhv.itrack_tr_1000(k, 1):bhv.itrack_tr_1000(k, 2);
            bhv.p_x_d(idx) = bhv.d_x_1000(bhv.itrack_tr_1000(k - 1, 2)) + bhv.d_x_1000(idx(1)) - bhv.d_x_1000(idx);
            bhv.p_x_ds(idx) = bhv.d_x_1000s(bhv.itrack_tr_1000(k - 1, 2)) + bhv.d_x_1000s(idx(1)) - bhv.d_x_1000s(idx);
    end

end

%plot(bhv.p_x_d)
%hold on

bhv.nb_tracks = length(bhv.way_tr);
bhv.nb_tracks = length(bhv.itrack_tr_1000(:, 1));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%YOU CAN ADD THAT IT'S WORKING BUT NO IDEA WHAT IS DOING%%%%
%%%%%%FROM GEOFFREY ORIGINAL CODE%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% just for aesthetic purpose (no idea what he means by that)
for k = 1:(bhv.nb_tracks - 1)
    %crap = (bhv.itrack_tr_100(k + 1, 1) - bhv.itrack_tr_100(k, 2));
    if (bhv.itrack_tr_1000(k + 1, 1) - bhv.itrack_tr_1000(k, 2)) < 5
        idx = bhv.itrack_tr_1000(k, 2):bhv.itrack_tr_1000(k + 1, 1);
        bhv.p_x_d(idx) = bhv.p_x_d(idx(1));
        bhv.p_x_ds(idx) = bhv.p_x_ds(idx(1));
    end
end



%%%%%SAMPLING RATE%%%%%
bhv.prm.freq = 25000;

prm.rmap.idcond_t = condway;
prm.rmap.tbin = indextr;
%prm.rmap.idx_rem = []; %%%let like that for now

%%%NEW NEW NEW%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
prm.rmap.idx_rem = find(bhv.v_x_1000_s2 < 2); %%%NEW NEW NEW speed threshold for ratemaps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

prm.rmap.xbin = 0:floor(maxmaze/100):maxmaze;
prm.rmap.xbin_rem = 0 %%%normally it's 10



nbsess=length(unique([Traj.Cond]));
%%%%%%%%%%%OLD%%%%%%%%%%%%
%prm.rmap.nb_cond = nbsess*2;

%%%%%%%%%%%NEW NEW NEW NEW%%%%%%%%%%%%
prm.rmap.nb_cond = 10; %%%%NB : Mathilde �tait toujours sur 10 l� dessus m�me si on avait que 2 ou 3 conditions. Si tu as au max 6 cond tu te mets � 12 (c'est ce que PP a fait)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



prm.rmap.freq = 1000;
prm.rmap.ismooth = 10; 

prm.rmap.frmax = 1; %1.5 %1  %%%%% 
prm.rmap.frmean = 0.3; %0.5 %0.3
prm.rmap.pctspklap = 0.50;  %%% ancien 0.50

prm.pf.mtd = 'circular_shift'; %random %poisson %circular_shift

prm.pf.nb_rep = 500; %500 %1000

prm.pf.min_len = 3;
prm.pf.max_len = 45; %Max length %OLDrom = 45

prm.pf.max_btw_pf = 5; %max dist between 2 pour merged
prm.pf.max_ext_pf = 3; %oldrom=5 %etendre pf  si pval=0.30

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
prm.cellclass.frmaxinterneuron = 20; %OLD : 20Hz max pour pyr
prm.cellclass.duration = 0.6; %OLD 0,4

maxc = 5;  %%%%%THIS IS GOOD (PUT 6 IF MAX 6 CONDITIONS)


%%%%IMPORTANT%%%%%%
%%%%NB of condition%%%%%%
eprm.nb_cond = nbsess; %%%%%THIS IS GOOD it is the REAL NUMBER OF CONDITIONS


NbCell = length(allcel.id_cel);
if run_only_rmap_pf
    rmap = cell(NbCell, 1);
    pf = cell(NbCell, 1);
else
    rmap = struct([]);
    pf = struct([]);
end
for g = 1:NbCell
    g
    idx_spk = allcel.itime_spk(allcel.id_spk == allcel.id_cel(g));
    idx_spk = fct_ifreq_swap(idx_spk,bhv.prm.freq,1000);
    %[rmap(g),pf(g)] = fct_mapandfields(X_ds_n,idx_spk,prm);

    %%%%%%%%%%%%%%%%%%%%%NEW NEW NEW%%%%%%%%%%%%%%%%%%%%%%%%
    [rmap_g, pf_g] = fct_mapandfields(bhv.p_x_ds,idx_spk,prm);
    if run_only_rmap_pf
        rmap{g} = rmap_g;
        pf{g} = pf_g;
    else
        rmap(g) = rmap_g;
        pf(g) = pf_g;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
end


%%%%to save the plots
%
%path_save=('X:\julie\DATA\MM1\data processed\MM1_2024-03-19_15-05-51')
%fct_plot_rmap(rmap,pf,'save',AnlPath)

%%%%%%%%%

close all
allcel.nb_cel = NbCell;
if run_only_rmap_pf
    % In batch PF-only mode we can skip reshape (some sessions/cells may have
    % heterogeneous rmap field sizes that break cat(3) in fct_mapandfields_reshape).
    save(fullfile(savepath, [datalc.session_name '_ratemap_circular']), 'allcel', 'rmap', 'pf','bhv');
    disp('Saved ratemap/pf only; skipping reshape + downstream cell-feature/classification blocks.')
    close all
    continue
end

[allrmap, allpf] = fct_mapandfields_reshape(rmap, pf);  %D�tection PF
%% Codage Distance/Position

% allrmap.feat.scd_nr_uc = NaN(allcel.nb_cel, maxc);
% allrmap.feat.scp_nr_uc = NaN(allcel.nb_cel, maxc);
% 
% allrmap.feat.inddis_nr_uc = NaN(allcel.nb_cel, maxc);
% allrmap.feat.indpos_nr_uc = NaN(allcel.nb_cel, maxc);
% 
% 
% 
% for g = 1:allcel.nb_cel
%     for c = 1:eprm.nb_cond
%         ind_c = (2*(c-1) + 1):(2*(c-1) + 2);
%         [allrmap.feat.scd_nr_uc(g, c), allrmap.feat.scp_nr_uc(g, c)] = fct_spatial_corr(rmap(g).fr_s_cx(ind_c(1), :), rmap(g).fr_s_cx(ind_c(2), :));
%         [allrmap.feat.inddis_nr_uc(g, c), allrmap.feat.indpos_nr_uc(g, c)] = fct_distposindexpeaktopeak(squeeze(allpf.feat.ifrmax_ucf(g, ind_c(1), :)), squeeze(allpf.feat.ifrmax_ucf(g, ind_c(2), :))); 
%     end
% end

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
% for g= 1:allcel.nb_cel
%  subplot(1,2,1) 
%   bar(allcel.acg1_bin_u(:,g))
% 
% subplot(1,2,2)
%  bar(allcel.acg10_bin_u(:,g))
%  
%  title(strcat('cell ', string(g)))
% %saveas(gcf, fullfile(AnlPath, ['ACG_' char(sess) '_cell_' num2str(g) '.png']));
% end
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

scatter([1:length(allcel.duration_u)],allcel.duration_u)
nRows = ceil(sqrt(allcel.nb_cel));      % nombre de lignes
nCols = ceil(allcel.nb_cel / nRows);    % nombre de colonnes
%{
figure;
for i = 1:allcel.nb_cel
    subplot(nRows, nCols, i);
    % Ton code de plot ici, par exemple :
     plot(allcel.bestwaveform(:,i))
    title(['Cellule ' num2str(i)]);
    
end
saveas(gcf,strcat(AnlPath,'\waveforms_new'),'png')
%}
%% Type de Cellule (Pyramidale / Interneurone) OLD
% 
allcel.type_u = ones(length(allcel.id_cel), 1);

ind_tmp1 = find(nansum(squeeze(nanmean(allrmap.fr_s_cxu, 2)) > prm.cellclass.frmaxinterneuron, 1) > 0);
 ind_tmp2 = find((allcel.burst_u < prm.cellclass.burstind & allcel.duration_u*1000 < prm.cellclass.duration));
% ind_tmp2 = find((allcel.burst_u < prm.cellclass.burstind &allcel.halfwidth < 0.0004));
ind_tmp = [ind_tmp1 ind_tmp2];
ind_tmp = unique(ind_tmp);

if ~isempty(ind_tmp)
    allcel.type_u(ind_tmp) = 0;   
end
%NEW
allcel.type_u_old=allcel.type_u;
allcel.type_u={}
% bursthappenat=[]
dur_u=[]
% burstindex=[]
% burstindexIN=[]
% for i = 1:allcel.nb_cel
%     [peaki,peakid]=max(allcel.acg1_bin_u(:,i));
%     bursthappenat=[bursthappenat;51-peakid];
%     baseline=mean(allcel.acg1_bin_u(1:11,i))
%     peakacg=mean(allcel.acg1_bin_u(51:56,i))
%     peakacgIN=mean(allcel.acg1_bin_u(51:65,i))
%     if peakacg-baseline>=0
%         burstindex=[burstindex;(peakacg-baseline)/peakacg]
%         burstindexIN=[burstindexIN;(peakacgIN-baseline)/peakacgIN]
%     else
%         burstindex=[burstindex;(peakacg-baseline)/baseline]
%         burstindexIN=[burstindexIN;(peakacgIN-baseline)/baseline]
%     end
% end
% 
% %subtracting (baseline) from the peak.negative amplitudes were normalized to the baseline to obtain indexes ranging from ?1 to 1.
% isburstindex=burstindex>=0
% isburstindexIN=burstindexIN>=0.3
% isburstearly=bursthappenat>=1.5 & bursthappenat<=5;
% isfr5=allcel.fr_u<=10;
% isfr20=allcel.fr_u>=20;
% isfr0=allcel.fr_u<5;
% isfr520=allcel.fr_u>=5&allcel.fr_u<=20;
isdur3=allcel.duration_u*1000> prm.cellclass.duration
%%
tau_rise_all = [];
tau_decay_all = [];
A_all = [];
C_all = [];

spkid = unique(allcel.id_spk);
%%
CV=[]
for g = 1:allcel.nb_cel
    spike_times = allcel.time_spk(allcel.id_spk == spkid(g));
% 
%     % --- Appel de la fonction de fit (d�finie plus bas dans ton script ou dans un fichier s�par�) ---
%  [beta, x_fit, y_fit, y_model, t0_auto,tau_rise] = fit_acg_risetime3_shift(spike_times);
% 
% tau_rise_all(g)  = tau_rise;
% tau_decay1_all(g)= beta(3)*1000;
% tau_decay2_all(g)= beta(5)*1000;
% 
% t0_all(g)= t0_auto * 1000;  % en ms
%     % Affiche pour suivi
%     fprintf('Cellule %d : tau_rise = %.3f ms, tau_decay = %.3f ms\n', ...
%         g, tau_rise, beta(3)*1000);

ISI = diff(spike_times);
CV2 = 2 * abs(ISI(2:end) - ISI(1:end-1)) ./ (ISI(2:end) + ISI(1:end-1));
CV2_mean = mean(CV2);
CV=[CV;CV2_mean];

% 
% figure(850)
% scatter(tau_rise_all,[1:length(tau_rise_all)],'o')
% set(gca, 'XScale', 'log')
%  line([2 2],[0 length(tau_rise_all)])
% title('Tau Rise ACG')
% figure(851)
% plot(CV,'o')
% line([0 length(CV)], [1 1])

% tau_rise_mask=tau_rise_all<=5;
CV2_mask=CV>=1;
 [peak1wf_id,peak1wf]=LocalMaxima(allcel.bestwaveform(1:22,g));
    if isempty(peak1wf)
        peak1wf_id=1;
        peak1wf=allcel.bestwaveform(peak1wf_id,g);
    else
        peak1wf_id=peak1wf_id(end)-1;
        peak1wf=peak1wf(end);
    end
    
    [throughwf,throughwf_id]=min(allcel.bestwaveform(peak1wf_id:30,g));
    throughwf_id=throughwf_id+peak1wf_id-1;
    
    [peak2wf_id,peak2wf]=LocalMaxima(allcel.bestwaveform(throughwf_id:end,g));
    if isempty(peak2wf_id)
        peak2wf_id=size(allcel.bestwaveform,1);
        peak2wf=allcel.bestwaveform(peak2wf_id,g);
    else
        peak2wf_id=peak2wf_id(1)-1;
        peak2wf=peak2wf(1);
        peak2wf_id = peak2wf_id+throughwf_id;
    end
    peak2through_dur=peak2wf_id-throughwf_id;
        peak2through_all=[peak2through_all;peak2through_dur];
end

    shortwf=peak2through_all<10;
    longwf=~shortwf;
% halfwidthmask=[]
% for i = 1:allcel.nb_cel
%     [max1,max1t]=max(allcel.bestwaveform(15:25,i));
%     max1t=max1t(1)+15;
%     [trough,trought]=min(allcel.bestwaveform(:,i));
%     [max2t,max2]=LocalMaxima(allcel.bestwaveform(trought:end,i));
%    % halfwidtht=ismember(allcel.bestwaveform(1:trough,i),max1t./2)
%     [halfinc halfinct] = min(abs(allcel.bestwaveform(1:trought,i)-max1t./2))
%     [halfdec halfdect] = min(abs(allcel.bestwaveform(trought:end,i)-max1t./2))
%     %halfwidth=allcel.bestwaveform(halfwidtht) % pyramidal cells generate action potentials (APs) with a half-width of about 2 ms, which are longer than those of most gamma
%     halfwidth=halfdect-halfinct;
%     halfwidthmask=[halfwidthmask;halfwidth>=3];
%     if isempty(max2t)
%         max2t=51;
%         max2=allcel.bestwaveform(end,i)
%     end
%     max2t=max2t(1)+trought;
%     dur_u=[dur_u;max2t-max1t];
%     plot(allcel.bestwaveform(:,i))
% end
% isdur28=dur_u>=2&dur_u<=8
% 
% isfr0=allcel.fr_u<5;
% isfr520=allcel.fr_u>=5&allcel.fr_u<=20;
% isfrIN=isfr520|isfr20
% 
isfrlow=allcel.fr_u<=15
ispyr = zeros(allcel.nb_cel,1);
ispyrlogic=CV2_mask&isfrlow.'&longwf;
figure(852)
scatter(allcel.fr_u(1,ispyrlogic),[1:length(find(ispyrlogic==1))],'r','filled')
hold on
scatter(allcel.fr_u(1,~ispyrlogic),[1:length(find(ispyrlogic==0))],'k')
title('FR mean of int and pyr cells')


inds = find(ispyrlogic == 1);
ncols = ceil(sqrt(numel(inds)));
nrows = ceil(numel(inds) / ncols);
figure(853)
hold off
for j = 1:numel(inds)
    f = inds(j);
    subplot(nrows, ncols, j);
    hold on
    plot(allcel.acg1_bin_u(:, f));
    title(sprintf('Cell Pyr %d', f));
end
saveas(gcf,strcat(fnames(1:end-27),'_pyrACG.png'))

inds = find(ispyrlogic == 0);
ncols = ceil(sqrt(numel(inds)));
nrows = ceil(numel(inds) / ncols);
figure(854)
hold off
for j = 1:numel(inds)
    f = inds(j);
    subplot(nrows, ncols, j);
    hold on
    plot(allcel.acg1_bin_u(:, f));
    title(sprintf('Cell Int %d', f));
end
saveas(gcf,strcat(fnames(1:end-27),'_intACG.png'))
% ispyr(isfr0)=0;
% % ispyrlogic=isdur28.'&isdur3&isfr5&isburstearly.'&isburstindex.'
% conds = [isdur28(:), isdur3(:), isfr5(:), isburstearly(:), isburstindex(:),halfwidthmask(:)];
% % somme du nombre de conditions vraies par cellule
% nb_true = sum(conds, 2);
% % crit�re : au moins 3 conditions vraies
% ispyrlogic = nb_true >= 5;
% ispyr(isfrIN)=3;
% ispyr(isfrIN&isburstindexIN.')=2;
ispyr(ispyrlogic)=1;
% 
allcel.type_u=ispyr;
%% Active/Silent

allcel.state_uc = NaN(allcel.nb_cel, maxc*2);
bhv.icondw_tr = condway;

% for g = 1:allcel.nb_cel
%     for c = 1:(eprm.nb_cond*2)
%         is_c = bhv.icondw_tr == c;
%         ratem =  rmap(g).fr_s_tx(is_c, :);  
%         allcel.state_uc(g, c) = double(fct_isactive(ratem, prm.rmap));
%     end
% end
% 

% D?termination des crit?res

crit_pyr_u = (allcel.type_u == 1);
% crit_pyr_u=allcel.type_u
crit_pyr_uc = zeros(allcel.nb_cel, maxc);
crit_pyr_uc(:, 1:eprm.nb_cond) = repmat(crit_pyr_u, 1, eprm.nb_cond);

crit_itn_u = (allcel.type_u == 0|allcel.type_u == 2|allcel.type_u == 3);                  
crit_itn_uc = zeros(allcel.nb_cel, maxc);
crit_itn_uc(:, 1:eprm.nb_cond) = repmat(crit_itn_u, 1, eprm.nb_cond);

crit_pf_stable_uc = allpf.feat.ispf_ucf(:, :, 1);

crit_active_uc = false(allcel.nb_cel, maxc);
crit_inactive_uc = false(allcel.nb_cel, maxc);
crit_bid_uc = false(allcel.nb_cel, maxc);
crit_uni_uc = false(allcel.nb_cel, maxc);
crit_nonsm_uc = false(allcel.nb_cel, maxc);

crit_itn_active_uc = false(allcel.nb_cel, maxc);    
crit_itn_inactive_uc = false(allcel.nb_cel, maxc);
crit_itn_bid_uc = false(allcel.nb_cel, maxc);
crit_itn_uni_uc = false(allcel.nb_cel, maxc);
crit_itn_nonsm_uc = false(allcel.nb_cel, maxc);
allrmap.feat.fr_c=[]; 

% Recalcul nombre cellules -> diff�rencier cellules inactives des cellules que l'on perd (pyr_off/itn_off)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%OLD OLD OLD OLD%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for c = 1:eprm.nb_cond
    ind_c = (2*(c-1) + 1):(2*(c-1) + 2);
    crit_active_uc(:, c) = (sum(allcel.state_uc(:, ind_c), 2) >= 1); % Active dans au moins un sens
    crit_inactive_uc(:, c) = (sum(allcel.state_uc(:, ind_c), 2) == 0);  % Non active dans les deux sens
    
    crit_bid_uc(:, c) = crit_pf_stable_uc(:, ind_c(1)) & crit_pf_stable_uc(:, ind_c(2));
    crit_uni_uc(:, c) = xor(crit_pf_stable_uc(:, ind_c(1)), crit_pf_stable_uc(:, ind_c(2)));
    crit_nonsm_uc(:, c) = ~crit_pf_stable_uc(:, ind_c(1)) & ~crit_pf_stable_uc(:, ind_c(2));
%     fr_c =  nanmean(allrmap.feat.fr_uc ( : , c : c+1 ),2);
%     allrmap.feat.fr_c= [allrmap.feat.fr_c fr_c]
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%NEW NEW NEW NEW%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%WE PUT THAT BACK!!!!!!! IF 6 COND MAX YOU ADD an extra 'u' 

%JAI AUTOMATISE POUR QUE CA VARIE AVEC NB DE CONDITIONS VOIR AU DESSUS
p = nanmean(allrmap.feat.fr_uc ( : , 1 : 2 ),2);
q = nanmean(allrmap.feat.fr_uc ( : , 3 : 4 ),2);
r = nanmean(allrmap.feat.fr_uc ( : , 5 : 6 ),2);
s = nanmean(allrmap.feat.fr_uc ( : , 7 : 8 ),2);
t = nanmean(allrmap.feat.fr_uc ( : , 9 : 10 ),2);
%u = nanmean(allrmap.feat.fr_uc ( : , 11 : 12 ),2);

allrmap.feat.fr_c = [p q r s t];
%allrmap.feat.fr_c = [p q r s t u];

lost_cell_c = zeros(size(allrmap.feat.fr_c));
lost_cell_c(find(allrmap.feat.fr_c == 0)) = 1;

crit_present_cell_uc = (lost_cell_c == 0);
crit_lost_cell_uc = (lost_cell_c == 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%NEW NEW NEW NEW%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%THIS NEXT 2 LINES ARE NOT NECESSARY%%%%%%%%%%%%%%%%%

%crit_present_cell_uc = [crit_present_cell_uc,zeros(NbCell,5-size(crit_present_cell_uc ,2))];
%crit_lost_cell_uc=[crit_lost_cell_uc,zeros(NbCell,5-size(crit_lost_cell_uc ,2))];

%si n�cessaire : crit_present_cell_uc=[crit_present_cell_uc, zeros(30,3)]

% D?termination des cellules 
allcel.pyr_uc = crit_pyr_uc ==1& crit_present_cell_uc;
allcel.pyr_active_uc = crit_pyr_uc==1 & crit_present_cell_uc & crit_active_uc;
allcel.pyr_inactive_uc = crit_pyr_uc ==1& crit_present_cell_uc & crit_inactive_uc;
allcel.pyr_off_uc = crit_pyr_uc==1 & crit_lost_cell_uc;

allcel.pyr_sm_uc = crit_pyr_uc==1 & crit_present_cell_uc & crit_active_uc & (crit_bid_uc | crit_uni_uc); % active is not necessary since bid and uni cells are based on a place field which is never determined for silent cell
allcel.pyr_bid_uc = crit_pyr_uc==1 & crit_present_cell_uc & crit_active_uc & crit_bid_uc;
allcel.pyr_uni_uc = crit_pyr_uc==1 & crit_present_cell_uc & crit_active_uc & crit_uni_uc;
allcel.pyr_nonsm_uc = crit_pyr_uc==1 & crit_present_cell_uc & crit_active_uc & crit_nonsm_uc;

allcel.itn_uc = crit_itn_uc & crit_present_cell_uc;
allcel.itn_active_uc = crit_itn_uc & crit_present_cell_uc & crit_active_uc;
allcel.itn_inactive_uc = crit_itn_uc & crit_present_cell_uc & crit_inactive_uc;
allcel.itn_off_uc = crit_itn_uc & crit_lost_cell_uc;

allcel.itn_sm_uc = crit_itn_uc & crit_present_cell_uc & crit_active_uc & (crit_bid_uc | crit_uni_uc); % active is not necessary since bid and uni cells are based on a place field which is never determined for silent cell
allcel.itn_bid_uc = crit_itn_uc & crit_present_cell_uc & crit_active_uc & crit_bid_uc;
allcel.itn_uni_uc = crit_itn_uc & crit_present_cell_uc & crit_active_uc & crit_uni_uc;
allcel.itn_nonsm_uc = crit_itn_uc & crit_present_cell_uc & crit_active_uc & crit_nonsm_uc;

% Nombre de cellules et pourcentages 
allcel.nb_cel_uc = allcel.nb_cel*ones(allcel.nb_cel , maxc);
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



%save([savepath filesep datalc.session_name '_Ratemap'], 'allrmap', 'allpf', 'allcel', 'rmap', 'pf');

%%%%%%NEW NEW NEW%%%%%%%%%%%%%%
%%%%%%Add BHV%%%%%%%%%%%%%%%%%%%
save(fullfile(savepath, [datalc.session_name '_Ratemap_final_CASE1']), 'allrmap', 'allpf', 'allcel', 'rmap', 'pf','bhv');
clearvars -except folder session_list 
close all
end
%yo=1




function filepath = find_session_file(folder, session_name, suffixes)
    filepath = '';

    % 1) Exact candidates first
    for ii = 1:numel(suffixes)
        candidate = fullfile(folder, [session_name suffixes{ii}]);
        if exist(candidate, 'file') == 2
            filepath = candidate;
            return
        end
    end

    % 2) Fallback: search any *.mat starting with session name + matching suffix
    files = dir(fullfile(folder, [session_name '*.mat']));
    for jj = 1:numel(files)
        fname = files(jj).name;
        for ii = 1:numel(suffixes)
            sfx = suffixes{ii};
            if length(fname) >= length(sfx)
                tail = fname(end-length(sfx)+1:end);
                if strcmpi(tail, sfx)
                    filepath = fullfile(folder, fname);
                    return
                end
            end
        end
    end

    error('Missing file for session "%s" in "%s". Expected suffixes: %s', ...
        session_name, folder, strjoin(suffixes, ', '));
end

%% --- Fonction utilitaire pour calculer l�autocorr�logramme ---
function [acf, lags] = autocorr_spikes(spike_times, binSize, maxLag)
    edges = -maxLag:binSize:maxLag;
    acf = zeros(size(edges)-[0 1]);
    for i = 1:length(spike_times)
        diffs = spike_times - spike_times(i);
        diffs(abs(diffs) < binSize) = []; % supprime le z�ro central
        acf = acf + histcounts(diffs, edges);
    end
    lags = edges(1:end-1) + binSize/2;
    acf = acf / length(spike_times); % normalisation
end

%%
function [beta, x_fit, y_fit, y_model, t0_auto, tau_rise] = fit_acg_risetime3_shift(spike_times, showPlot)
% FIT_ACG_RISETIME3_SHIFT - Fit de l'autocorr�logramme d'un train de spikes
% avec estimation du t0 automatique et de tau_rise (temps pour atteindre 63% de l'amplitude max)
%
% Entr�es :
%   spike_times : vecteur des temps de spikes (en secondes)
%   showPlot    : (optionnel) true/false pour afficher le fit
%
% Sorties :
%   beta      : param�tres du fit (1x5)
%   x_fit     : axe temporel du fit
%   y_fit     : donn�es normalis�es (partie positive du lag)
%   y_model   : mod�le ajust�
%   t0_auto   : temps de d�but de mont�e d�tect� (s)
%   tau_rise  : temps de mont�e 63% (ms)

    if nargin < 2
        showPlot =true;
    end

    % --- 1. Param�tres ---
    binSize = 0.001; % 1 ms
    maxLag  = 0.1;   % 100 ms

    % --- 2. Calcul de l'autocorr�logramme ---
    [acf, lags] = autocorr_spikes(spike_times, binSize, maxLag);

    % Partie positive seulement
    idx = lags > 0;
    x = lags(idx);
    y = acf(idx);
    y = y / max(y); % normalisation � 1

    % --- 3. D�tection automatique du t0 ---
    y_smooth = smoothdata(y, 'gaussian', 5); % lissage l�ger pour �viter le bruit
    threshold = 0.1 * max(y_smooth);
    idx_start = find(y_smooth > threshold, 1, 'first');

    if isempty(idx_start)
        t0_auto = x(1);
    else
        t0_auto = x(idx_start);
    end

    % --- 4. Calcul du tau_rise (temps pour atteindre 63% du max) ---
    target_amplitude = 0.63 * max(y);
    idx_above_threshold = find(x > t0_auto & y >= target_amplitude, 1, 'first');

    if ~isempty(idx_above_threshold)
        if idx_above_threshold > 1
            x_before = x(idx_above_threshold - 1);
            y_before = y(idx_above_threshold - 1);
        else
            x_before = x(idx_above_threshold);
            y_before = y(idx_above_threshold);
        end

        x_after = x(idx_above_threshold);
        y_after = y(idx_above_threshold);

        % interpolation lin�aire entre les deux points autour du 63%
        tau_rise_sec = x_before + ...
            (target_amplitude - y_before) * (x_after - x_before) / (y_after - y_before);

        tau_rise = (tau_rise_sec - t0_auto) * 1000; % conversion en ms

        % s�curit� : valeurs n�gatives ou aberrantes mises � NaN
        if tau_rise < 0 || tau_rise > 100
            tau_rise = NaN;
        end
    else
        tau_rise = NaN;
    end

    % --- 5. Fit du mod�le ---
    modelFun = @(b, x) (x > t0_auto) .* ( ...
        b(1) * ((1 - exp(-(x - t0_auto)./max(b(2),1e-6))) .* ...
                exp(-(x - t0_auto)./max(b(3),1e-6))) + ...
        b(4) * exp(-(x - t0_auto)./max(b(5),1e-6)) ...
    );

    % bornes et conditions initiales
    b0 = [1, 0.004, 0.02, 0.5, 0.08];
    lb = [0, 1e-4, 1e-4, 0, 1e-4];
    ub = [Inf, 0.02, 0.2, Inf, 0.2];

    opts = optimoptions('lsqcurvefit', 'Display', 'off', 'MaxFunctionEvaluations', 1e5);
    beta = lsqcurvefit(modelFun, b0, x, y, lb, ub, opts);

    % --- 6. R�sultat ---
    x_fit = linspace(0, max(x), 300);
    y_model = modelFun(beta, x_fit);
    y_fit = y;

    % --- 7. Affichage optionnel ---
    if showPlot
        figure; 
        plot(x*1000, y, 'k', 'LineWidth', 1.2);hold on;
        plot(x_fit*1000, y_model, 'r--', 'LineWidth', 2);

        % lignes de rep�re
        line([t0_auto*1000, t0_auto*1000], ylim, 'Color', [0 0.5 0], 'LineStyle', '--', 'LineWidth', 1.5);
        if ~isnan(tau_rise)
            line([t0_auto*1000 + tau_rise, t0_auto*1000 + tau_rise], ylim, 'Color', 'b', 'LineStyle', '--', 'LineWidth', 1.5);
        end

        xlabel('Lag (ms)');
        ylabel('Autocorr�logramme normalis�');
        title('Fit de l''ACG avec t_0 et \tau_{rise}');
        legend({'ACG r�el', 'Fit', 't_0', '\tau_{rise}'}, 'Location', 'best');
        grid on; hold off;
    end
end
