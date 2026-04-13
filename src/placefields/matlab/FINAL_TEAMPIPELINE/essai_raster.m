fpath = 'N:\Vinca\MATLAB\Data\All_data';
main_folder='N:\Vinca\MATLAB\Data\'
h=readtable(strcat(main_folder,'select_sess.xlsx'))
keep_dat=h.session_name(logical(h.doselect))


load('N:\Vinca\MATLAB\Data\All_data\VS68_2023-02-22_16-56-07_Bhv0.mat')
load('N:\Vinca\MATLAB\Data\All_data\VS68_2023-02-22_16-56-07_Ratemap.mat')
for i = 16
behav=strcat(fpath,'\',string(keep_dat(i)),'_Bhv0.mat')    
ratemap=strcat(fpath,'\',string(keep_dat(i)),'_Ratemap.mat')
load(behav)
load(ratemap)    
for c=1:max(bhv.icondw_tr)
    for g =1:size(allrmap.nbspk_txu,3)

       
% 
% g =1; 
% c = 2;

vel_th = 2;

icond = bhv.icondw_tr == c;

idspk = allcel.itime_spk(allcel.id_spk == allcel.id_cel(g));
idspk = fct_ifreq_swap(idspk, 25000, 100);
idspk(bhv.v_x_ds2(idspk) < vel_th) = [];
idspk=double(idspk)
idtr = bhv.idtrack_tr(icond, :);
ispktr = fct_discretize(idspk, idtr, 'bin');
idspk = idspk(~isnan(ispktr));
ispktr = ispktr(~isnan(ispktr));
posspk = bhv.p_x_d(idspk);

mat = allrmap.fr_s_txu(icond, :, g);
mat = bsxfun(@rdivide, mat, nanmax(mat, [], 2));
pfmat = allpf.ispf_cxu(c, :, g);
meanfr = allrmap.fr_s_cxu(c, :, g);
nb_trial = sum(icond);

% color1 = [63 215 231] / 255;
% color2 = [42 93 229] / 255;
% color_set = [color1 ; color2];
color_set = [65 97 120 ; 65 97 120] / 255;
figure
subplot(3, 1, 1)

fprm = fdt_fprmlist('raster');
fprm.barcol = color_set(1, :);
fprm.xlab = '';
fprm.xt = [];
fct_plot_raster(posspk, ispktr, nb_trial, fprm)


subplot(3, 1, 2)
fprm = fdt_fprmlist('heatmap');
fprm.xlab = '';
% fprm.xt = fprm.xt - 10;
fprm.xt = [];
fprm.xlim = [0 80];

fct_plot_ratemap(mat, fprm);



subplot(3, 1, 3)

fprm = fdt_fprmlist('mean');
fprm.markersize = 5;
fprm.color = color_set(1, :);
fprm.pf_color = fprm.color;
fprm.xt = fprm.xt - 10;
fprm.yt = 0:2:14;
fprm.ylim = [0 max(meanfr+5)];
fprm.xlim = [0 80];
fct_plot_meanrate(meanfr, [], pfmat, [], fprm);
    end
   close all
end 
end