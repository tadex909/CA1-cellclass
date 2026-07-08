function [samples_chnl_pts_spk_cell, samples_chnl_pts_cell] = Waveform_sample(Wave_chnl_pts_spk, id_spk )
%%%select the 20 first waveform of all cells
%nbchperEl = number of channels per shank

%Wave_chnl_pts_spk = matrice with number of recording sites * voltage (32 ms) * spike
%samples

%the 32 ms are decided in the params.prm file used to run klusta


%id_spk = index of cell identity (ordered by the timestamps)

%waveform.samples (channel#, 32 points, 20 spikes, cell index)
%waveform.average (channel#, 32 points, cell index)



S = size(Wave_chnl_pts_spk);
nSamplePerSpk = S(2); 

samples=[];
average=[];

%number of cells
G = unique(id_spk);
%Wave_chnl_pts_spk=permute(Wave_chnl_pts_spk,[3 1 2]);

%waveSampl=[];
%waveMean=[];

for t = 1:length(G) % Parcours le nombre de Clusters
    
   spkndx = find(id_spk == G(t));
   tp = Wave_chnl_pts_spk(:, :, spkndx);
   
   if length(spkndx) > 20     %%% get 20 first spikes for each cell
      tp = tp(:, :, 1:20);
   else
      tp = tp(:, :, 1:end);
   end
   
   
   samples_chnl_pts_spk_cell(:,:,:,t)=tp;
   samples_chnl_pts_cell(:,:,t) = mean(tp,3); % Plot 8 courbes (autant que de nombre de canaux ) sur la moyenne de chaque waveform pour un shk donné. Moyenne effectuée sur au plus 20 spikes
   
%    plot((mean(tp,3))')
    
end




end

