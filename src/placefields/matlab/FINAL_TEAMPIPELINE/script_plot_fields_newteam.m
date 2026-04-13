function script_plot_fields(datalc)
%Load supool et Traj et Ratemap on peut utiliser i=supool.keptsession_idcel pour faire sess par sess
load colormaprom.mat
%% Comparaison all fields avant/après muscimol

allrmap.feat.ispyr_uc = kron(supool.pyr_uc, ones(1, 2));
ispyrpf = allpf.feat.ispf_ucf(:, :, 1) & allrmap.feat.ispyr_uc;

nb_cond = 4;
nb_cel = allcel.nb_cel;


A = cell(nb_cond, 1);
B = cell(nb_cond, 1);
for u = 1:nb_cel
    for c = 1:nb_cond
        if ispyrpf(u, c)  %si c une place cell 
        A{c} = [A{c} ; allrmap.fr_s_cxu(c, : , u)]; %on  met dans A {field condition condw} le FR de condw c par bin pour la cellule n° u (taille 1x80)
        
        if c < 3 % si il s'agit d'une condition différente de la cond 1 (condw 1et2)
        B{c} = [B{c} ; allrmap.fr_s_cxu(c, : , u)]; 
        B{c + 2}  = [B{c + 2} ; allrmap.fr_s_cxu(c + 2, : , u)]; 
        end
        end
    end
end

fprm = fdt_fprmlist('rankfield');
fprm.xt = fprm.xt - 10;
fprm.xlim = fprm.xlim - 10;
fprm.ylab = 'Place Field Number';
edgobj = [15 35 ; 80 90] - 10;

subplot(221)
fprm.colormap = colormapbnl;
fct_plot_rankfield([A{1};A{3};A{2};A{4}], fprm)
%title(datalc.session_name, 'fontsize', 15)
title('Condition ALL PO A/R')
subplot(223)
v = nanmean(fct_matnorm([A{1};A{3};A{2};A{4}]), 1);
plot(v)
hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
hold off


subplot(222)
fprm.colormap = colormapcnl;
fct_plot_rankfield([A{2};A{4}], fprm)
title('Condition 2 PO Retours')


subplot(224)
hold off
v = nanmean(fct_matnorm([A{2};A{4}]), 1);
plot(v)

hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');

fct_fullscreen(gcf)
%fct_save_figure(gcf, [datalc.save_path filesep datalc.session_name 'Maps1.jpg'], 'jpg')

%% Que deviennent les fields avec musci ?

figure
fprm = fdt_fprmlist('rankfield');
fprm.xt = fprm.xt - 10;
fprm.xlim = fprm.xlim - 10;
fprm.ylab = 'Place Field Number';

subplot(221)
fprm.colormap = colormapbnl;
[A_np, ind_perm] = fct_rank_pf([B{1} ; B{2}], 'cont');
fct_plot_rankfield(A_np, fprm)
%title(datalc.session_name, 'fontsize', 15)


subplot(223)
v = nanmean(fct_matnorm([B{1} ; B{2}]), 1);
plot(v)

hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');

subplot(222)
fprm.colormap = colormapcnl;
B_np = [B{3} ; B{4}];
B_np = B_np(ind_perm, :);
B_np = fct_matnorm(B_np, 2);
fct_plot_ratemap(B_np, fprm)


subplot(224)
v = nanmean(fct_matnorm([B{3} ; B{4}]), 1);
plot(v)

hold on
line([edgobj(1, 1) edgobj(1, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(1, 2) edgobj(1, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(2, 1) edgobj(2, 1)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');
line([edgobj(2, 2) edgobj(2, 2)], [min(v) max(v)], 'color', 'r', 'linestyle', '--');

fct_fullscreen(gcf)
fct_save_figure(gcf, [datalc.save_path filesep datalc.session_name 'Maps2.jpg'], 'jpg')








% idxf1 = permute(repmat(idx1, 1, 1, nb_bin), [2 3 1]);
% imagesc(reshape(allrmap.fr_s_cxu(idxf1), nb_bin, nb_pc)')




