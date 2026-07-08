
%% Se mettre dans fichier cible où on a continuous.dat et structure oebin 
%MAKE LFP FILE AT 1250HZ
[lfp, ~]=Read_OEP_Binary ('continuous.dat',[1:64], 0,-1, 1250);
save(['lfp.mat'],'lfp')

%% Se mettre dans fichier cible contenant lfp, ratemap etc
Session = 'VS99_2024-11-02_15-41-34';
Mouse =' VS99';
prob=64

LFP_path = [cd, '\lfp.mat'];
ePhy_path = [cd  '\' num2str(Session) '_ePhy.mat'];
Rmap_path = [cd '\' num2str(Session) '_Ratemap.mat'];

Save_path = 'Z:\EPSZTEIN\Vinca\MATLAB\Data\DeepSup_analyses';
mkdir(Save_path, Session);
Save_folder = strcat(Save_path, '\', Session);

lfp_samplingrate = 1250;
bestchannel = 57; %Choisir sur neuroscope le channel où le ripple power semble le plus fort

chorder_buz64 = [5 4 6 3 7 2 8 1];
chorder_buz5 = [7 6 8 5 12 1 11 2 10 3 9 4];
chorder_buzlin = [25 40 26 39 27 38 28 37 29 36 30 35 31 34 32 33];
chorder_assy=[]
if prob == 64
   best_shk= bestshk_64(bestchannel); 
   chorder=chorder_buz64
  %Détection temps ripples et calcul power
   [eeg, rippleT] = detectripple(LFP_path,bestchannel,lfp_samplingrate);
   rippleT = cleanripple(Session,Save_folder,lfp_samplingrate,best_shk,chorder,eeg,rippleT); %Entrer n° ripple à jeter sous forme de vecteur
elseif prob == 5
    best_shk= bestshk_5(bestchannel);
    %Détection temps ripples et calcul power
    [eeg, rippleT] = detectripple(LFP_path,bestchannel,lfp_samplingrate);
    if length(best_shk)== 16
        best_shk = best_shk-24
       %rippleT = cleanripple_buzlin(Session,Save_folder,lfp_samplingrate,best_shk,chorder_buzlin-24,eeg,rippleT); %Entrer n° ripple à jeter sous forme de vecteur
    elseif length(best_shk)== 12
        %rippleT = cleanripple_buzlin(Session,Save_folder,lfp_samplingrate,best_shk,chorder_buz5,eeg,rippleT); %Entrer n° ripple à jeter sous forme de vecteur
    end
end

%%working until line 42
[PowerRippleMat, eeg_reorder] = Ripplelayout2(Save_folder,Session,eeg,lfp_samplingrate,chorder,rippleT,prob); %Figure power ripples sur l'ensemble des shanks
%PowerRippleMat et eeg_reorder -> pour chaque shank l'ordre des canaux suit chorder (du plus profond (channel 1) au plus superficiel (channel 8))

figure_ripples(eeg,best_shk,chorder,lfp_samplingrate,rippleT,Save_folder,Session) %Enregistrer figures ripples 

%Position cellules %ICI ! on en est à cellXYZ ligne 80 !! mais il nous faut
%un VRAI rec buzlin pour debugger. 
CellXYZg = CellPosition_vee(prob,ePhy_path,chorder,lfp_samplingrate,Save_folder,Session); %Enregistrer figures position cellule 
DeepSupMat = DeepvsSup_cell(CellXYZg,PowerRippleMat); % -1 = Superficial ; 1 = Deep ; 0 = au milieu

%Seulement 4 ripples détectés sur 3ème session MA84
%Seulement 9 ripples détectés sur 1ère session MA90
%Trop peu de ripples et trop petits pour MA91 pour être détectés (on s'en fout ->1PC et pas de bid faaaatch)

%% Nombre cellules deep/superficial

load(Rmap_path);

deepsup.cell.ispyr = allcel.pyr_uc;
deepsup.cell.ispyract = allcel.pyr_active_uc;
deepsup.cell.ispyrsm = allcel.pyr_sm_uc;
deepsup.cell.ispyrbid = allcel.pyr_bid_uc;
deepsup.cell.ispyruni = allcel.pyr_uni_uc;

for i = 1:3
    deepsup.cell.nb_pyr_sup(i) = sum(DeepSupMat(find(allcel.pyr_uc(:,i) == 1))== -1);
    deepsup.cell.nb_pyr_deep(i) = sum(DeepSupMat(find(allcel.pyr_uc(:,i) == 1))== 1);
    deepsup.cell.nb_pyr_mid(i) = sum(DeepSupMat(find(allcel.pyr_uc(:,i) == 1))== 0);
    
    deepsup.cell.nb_actpyr_sup(i) = sum(DeepSupMat(find(allcel.pyr_active_uc(:,i) == 1))== -1);
    deepsup.cell.nb_actpyr_deep(i) = sum(DeepSupMat(find(allcel.pyr_active_uc(:,i) == 1))== 1);
    deepsup.cell.nb_actpyr_mid(i) = sum(DeepSupMat(find(allcel.pyr_active_uc(:,i) == 1))== 0);
    
    deepsup.cell.nb_pc_sup(i) = sum(DeepSupMat(find(allcel.pyr_sm_uc(:,i) == 1))== -1);
    deepsup.cell.nb_pc_deep(i) = sum(DeepSupMat(find(allcel.pyr_sm_uc(:,i) == 1))== 1);
    deepsup.cell.nb_pc_mid(i) = sum(DeepSupMat(find(allcel.pyr_sm_uc(:,i) == 1))== 0);

    deepsup.cell.nb_bid_sup(i) = sum(DeepSupMat(find(allcel.pyr_bid_uc(:,i) == 1))== -1);
    deepsup.cell.nb_bid_deep(i) = sum(DeepSupMat(find(allcel.pyr_bid_uc(:,i) == 1))== 1);
    deepsup.cell.nb_bid_mid(i) = sum(DeepSupMat(find(allcel.pyr_bid_uc(:,i) == 1))== 0);

    deepsup.cell.nb_uni_sup(i) = sum(DeepSupMat(find(allcel.pyr_uni_uc(:,i) == 1))== -1);
    deepsup.cell.nb_uni_deep(i) = sum(DeepSupMat(find(allcel.pyr_uni_uc(:,i) == 1))== 1);
    deepsup.cell.nb_uni_mid(i) = sum(DeepSupMat(find(allcel.pyr_uni_uc(:,i) == 1))== 0);
end

deepsup.ripple.rippleTime = rippleT;
deepsup.ripple.PowerRippleMat = PowerRippleMat;
deepsup.ripple.CellXYZg = CellXYZg;
deepsup.ripple.DeepSupMat = DeepSupMat;
deepsup.prm.Mouse = Mouse;
deepsup.prm.Session = Session;
deepsup.prm.samplingrate = lfp_samplingrate;
deepsup.prm.bestchannel = bestchannel;
deepsup.prm.bestshank = best_shk;
deepsup.prm.chorder = chorder;

save([Save_folder filesep Session '_DeepSup'], 'deepsup');
copyfile ([Save_folder filesep Session '_DeepSup.mat'], [save_path filesep Session])

%% ALL SESSIONS
% P{1} = 'X:\EPSZTEIN\Mathilde\Septum Inactivation\DeepvsSuperficial_analyses\MA47_2019-09-12_14-27-11\MA47_2019-09-12_14-27-11_Ripples.mat';
% P{2} = 'X:\EPSZTEIN\Mathilde\Septum Inactivation\DeepvsSuperficial_analyses\MA47_2019-09-13_15-21-43\MA47_2019-09-13_15-21-43_Ripples.mat';
% P{3} = 'X:\EPSZTEIN\Mathilde\Septum Inactivation\DeepvsSuperficial_analyses\MA47_2019-09-14_15-48-37\MA47_2019-09-14_15-48-37_Ripples.mat';
% P{4} = 'X:\EPSZTEIN\Mathilde\Septum Inactivation\DeepvsSuperficial_analyses\MA49_2019-09-30_15-36-29\MA49_2019-09-30_15-36-29_Ripples.mat';
% P{5} = 'X:\EPSZTEIN\Mathilde\Septum Inactivation\DeepvsSuperficial_analyses\MA53_2019-11-28_13-18-53\MA53_2019-11-28_13-18-53_Ripples.mat';
% P{6} = 'X:\EPSZTEIN\Mathilde\Septum Inactivation\DeepvsSuperficial_analyses\MA63_2020-06-05_15-15-05\MA63_2020-06-05_15-15-05_Ripples.mat';
% P{7} = 'X:\EPSZTEIN\Mathilde\Septum Inactivation\DeepvsSuperficial_analyses\MA85_2021-02-12_19-15-39\MA85_2021-02-12_19-15-39_Ripples.mat';
% P{8} = 'X:\EPSZTEIN\Mathilde\Septum Inactivation\DeepvsSuperficial_analyses\MA89_2021-03-18_12-25-56\MA89_2021-03-18_12-25-56_Ripples.mat';
% 
% supoolDS.Position_CellXYZ = [];
% supoolDS.DeepvsSup = [];
% 
% for s = 1 : length(P)
%     
%     load(P{s});
%        
%     supoolDS.Position_CellXYZ = cat(1,supoolDS.Position_CellXYZ, CellXYZg);
%     supoolDS.DeepvsSup = cat(1,supoolDS.DeepvsSup, DeepSupMat);
% end
% 
% save([Save_path filesep 'Supool_deepVSuperficial_PNO'], 'supoolDS')
% 
