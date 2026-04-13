clear path
path='X:\MATLAB\Data\All_data_final\process\OLD_ratemap';

cd(path);
B=dir();
cd('X:\MATLAB\Data\All_data_final');
h=readtable('Finalyses_Vinca.xlsx');
cd('X:\MATLAB\Data\All_data_final\process');
o=readtable('C:\Users\suire\Desktop\Waveform_CR.xlsx');
isgood=logical(o.is_good);
sessnames=char(o.Var1(~isgood));
%  for i=3:length(B)
for i=14:size(sessnames,1)
     sess=sessnames(i,2:end-14)

    idx = cellfun(@(x) ~isempty(strfind(x, sess)), h.session_name);
   AnlPath=char(string(h.path(idx)))
    load(strcat(sess,'_Phenosys.mat'));
    load(strcat(sess,'_TrajData.mat'));
    load(strcat(sess,'_ePhy.mat'));
    
 OutputFolder=cd;
% 
%     Namesession=sess;
%     nRows = ceil(sqrt(length(allcel.id_cel)));      % nombre de lignes
% nCols = ceil(length(allcel.id_cel) / nRows);    % nombre de colonnes
% 
% figure;
% for i = 1:length(allcel.id_cel)
%     subplot(nRows, nCols, i);
%     % Ton code de plot ici, par exemple :
%      plot(allcel.bestwaveform(:,i))
%     title(['Cellule ' num2str(i)]);
%     
% end
% saveas(gcf,strcat(AnlPath,'\waveforms_new'),'png')
    Rmap_G_vinca_modifJK_final
    clearvars -except B h sessnames

end