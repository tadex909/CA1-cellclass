Session = 'VS104_2024-11-10_16-39-46'
Save_folder = cd
%% Se mettre dans fichier cible où on a continuous.dat et structure oebin
%MAKE LFP FILE AT 1250HZ
[lfp, ~]=Read_OEP_Binary ('continuous.dat',[1:64], 0,-1, 20);
for k=1:64
dataset(k:64:length(lfp)*64) = int16(lfp(k,:));
dataset2(k:64:length(lfp)*64) = (lfp(k,:));
end
save(['lfp.mat'],'lfp')
%Save lfp.dat 1250HZ et 64channels.


%% FIND BEST CHANNEL ON NEUROSCOPE
bestchannel = 33
lfp_samplingrate=1250
prob=[77] %5,64,77?
%% DETECTION RIPPLE
disp('START ripple detection')
eeg = double(lfp);
eeg_bestch = eeg(bestchannel, :);

[yo,fo,to,do]=spectrogram(eeg_bestch,50,25,1:300,lfp_samplingrate);

tresh = 3;

nn=find(mean(do(121:230,:),1)>tresh*std(mean(do(121:230,:),1))&mean(do(121:230,:),1)>mean(do(251:300,:),1));
dd=find([2 diff(nn)>1]);
ndx=nn(dd);

rippleT=to(ndx);

disp('END ripple detection')

%% Clean ripple 
chorder_buz64 = [5 4 6 3 7 2 8 1; 13 12 14 11 15 10 16 9; 21 20 22 19 23 18 24 17; 29,28,30,27,31,26,32,25; 37,36,38,35,39,34,40,33; 45,44,46,43,47,42,48,41;53,52,54,51,55,50,56,49;61,60,62,59,63,58,64,57];
chorder_buz5 = [7 6 8 5 12 1 11 2 10 3 9 4;19,18,20,17,24,13,23,14,22,15,21,16;47,46,48,45,52,41,51,42,50,43,49,44;59,58,60,57,64,53,63,54,62,55,61,56];
% chorder_buzlin = [25 40 26 39 27 38 28 37 29 36 30 35 31 34 32 33];
chorder_assy=[22,28,32,26,25,17,29,19,30,27,24,20,23,18,21,16;2,4,6,8,10,12,9,14,7,31,5,15,3,13,1,11;63,61,59,57,55,53,56,51,58,34,60,50,62,62,64,54;43,37,33,39,40,48,36,46,35,38,41,45,42,47,44,49];
if prob == 64
   [best_shk,shk]= bestshk_64(bestchannel); 
   chorder=chorder_buz64(shk,:)
   rippleT = cleanripple_new(Session,Save_folder,lfp_samplingrate,best_shk,chorder,eeg,rippleT); %Entrer n° ripple à jeter sous forme de vecteur
elseif prob == 5
    [best_shk,shk]= bestshk_5(bestchannel);
       chorder=chorder_buz5(shk,:)
    rippleT = cleanripple_new(Session,Save_folder,lfp_samplingrate,best_shk,chorder,eeg,rippleT); %Entrer n° ripple à jeter sous forme de vecteur
elseif prob == 77
    [best_shk,shk]=bestshk_77(bestchannel)
       chorder=chorder_assy(shk,:)
    rippleT = cleanripple_new(Session,Save_folder,lfp_samplingrate,best_shk,chorder,eeg,rippleT); %Entrer n° ripple à jeter sous forme de vecteur
end

