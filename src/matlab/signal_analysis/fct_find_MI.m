% -----------------------------
% Written by MARTI Geoffrey and MICHON FX
% 10/16
% 02/16 : vararin -> prm and wrap to find MI in intervals
% -----------------------------


function MI = fct_find_MI(v_x, prm)




is_M = v_x > prm.thresh_speed;



MI.ind_M = fct_find_seq(is_M);
MI.ind_M = fct_merge_seq(MI.ind_M, prm.intermov_dur);
[MI.ind_M, MI.nb_M, MI.length_M] = fct_remove_seq(MI.ind_M, 'minlen', prm.duration);


is_I = true(1, length(v_x));
for k = 1:MI.nb_M % Les périodes de mouvement ne peuvent pas être des périodes d'arrêt
    ind = MI.ind_M(k, 1):MI.ind_M(k, 2);
    is_I(ind) = false;
end

is_I(v_x > prm.thresh_speed) = false;
[MI.ind_I, MI.nb_I, MI.length_I] = fct_find_seq(is_I, 'minlen', prm.duration);
[MI.ind_I, MI.nb_I, MI.length_I] = fct_merge_seq(MI.ind_I, prm.intermov_dur);



% Transition Periods
% crit_length = 50000;
crit_length = prm.duration;


MI.ind_M2I = [];
if (MI.nb_M == 0) || (MI.nb_I == 0)
    MI.nb_M2I = 0;
else
    p = 1;
    for k = 1:MI.nb_M
        tmp = (MI.ind_I(:, 1) - MI.ind_M(k, 2) <= crit_length) & ((MI.ind_I(:, 1) - MI.ind_M(k, 2)) >= 0);
        if any(tmp)
            ind_I_tmp = find(tmp == true, 1, 'first');
            MI.ind_M2I(p, 1:2) = MI.ind_M(k, :);
            MI.ind_M2I(p, 3:4) = MI.ind_I(ind_I_tmp, :);
            p = p + 1;
        end
    end
    MI.nb_M2I = p - 1;
end


MI.ind_I2M = [];
if (MI.nb_M == 0) || (MI.nb_I == 0)
    MI.nb_I2M = 0;
else
    p = 1;
    for k = 1:MI.nb_I
        tmp = (MI.ind_M(:, 1) - MI.ind_I(k, 2) <= crit_length) & ((MI.ind_M(:, 1) - MI.ind_I(k, 2)) >= 0);
        if any(tmp)
            ind_M_tmp = find(tmp == true, 1, 'first');
            MI.ind_I2M(p, 1:2) = MI.ind_I(k, :);
            MI.ind_I2M(p, 3:4) = MI.ind_M(ind_M_tmp, :);
            p = p + 1;
        end
    end
    MI.nb_I2M = p - 1;
end


%% Analysis index of movement & immobility in track

%
% MI.ind_Mperiods = MI.seq_ind_m(find(MI.seq_length_m > prm.duration * freq),:);
% MI.nb_Mperiods = length(MI.ind_Mperiods);
%
% a= MI.seq_ind_m(1,1) - 1;
% for i=2:length(MI.seq_ind_m(:,1))
%     b(i)= MI.seq_ind_m(i,1)- MI.seq_ind_m(i-1,2);
% end
% b(1)=[];
% c =length(v_x) - MI.seq_ind_m (end:2);
% d = [a b c];
%
%
% e = find(d > prm.duration * freq);
% MI.ind_Iperiods = [MI.seq_ind_m(e - 1 ,2) MI.seq_ind_m(e,1)];
% MI.nb_Iperiods=length(MI.ind_Iperiods);










