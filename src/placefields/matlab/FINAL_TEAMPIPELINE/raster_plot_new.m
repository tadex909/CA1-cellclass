fpath = 'N:\Vinca\MATLAB\Data\All_data';
main_folder='N:\Vinca\MATLAB\Data\';
h=readtable(strcat(main_folder,'select_sess.xlsx'));
keep_dat=h.session_name(logical(h.doselect));
session_name='VS102_2024-11-05_18-59-31';

load(strcat(fpath,filesep,session_name,'_TrajData.mat'));
load(strcat(fpath,filesep,session_name,'_Ratemap.mat'));
%%
for i = 16

for c=1:max(vertcat(Traj.icondway_tr));
    for g =1:size(allrmap.nbspk_txu,3);

       
% 
% g =1; 
% c = 2;

vel_th = 2;
icondway_tr=[Traj.icondway_tr].';
icond = icondway_tr.'== c;
Speed=vertcat(Traj.XSpeed);

idspk = allcel.itime_spk(allcel.id_spk == allcel.id_cel(g));
idspk = fct_ifreq_swap(idspk, 25000, 1000);
idspk=double(idspk);
idtrack_tr(:,1)=vertcat(Traj.start);
idtrack_tr(:,2)=vertcat(Traj.stop);
idtr =idtrack_tr(icond, :);
ispktr = fct_discretize(idspk, idtr, 'bin');
idspk = idspk(~isnan(ispktr));
ispktr = ispktr(~isnan(ispktr));
ispktr(Speed(idspk)<vel_th)=[]
idspk(Speed(idspk)<vel_th)=[];







p_x_d=vertcat(Traj.VRtraj);
posspk = p_x_d(idspk);

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
fprm.xlim=[20;135]
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