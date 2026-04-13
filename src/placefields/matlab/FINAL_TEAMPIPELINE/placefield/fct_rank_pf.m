% -----------------------------
% Written by MARTI Geoffrey
% 10/15
% 12/15
% 03/16
% -----------------------------

%% Version avec Cells et Dynamic Fields
% function [A_np] = fct_rank_pf(struct_var, field_rate_map_name, field_name_pf_name, field_stability_name)
% % struct_var = mouse.cell;
% % field_rate_map_name = 'rate_s_nr_cond_bin';
% % field_name_pf_name = 'place_field_nr_cond_bin';
% % field_stability_name = 'stability_index_cond';
%
% % A ./ repmat(max(A_n')', 1, size(A_n, 2)); % La fonction de normalisation
% % par ligne existait déjà !!
%
% for c = 1:size(struct_var(1).(field_rate_map_name), 1)
%     k = 1;
%     for g = 1:numel(struct_var)
%         if any(struct_var(g).(field_name_pf_name)(c, :) == 1) && (struct_var(g).(field_stability_name)(c) >= 0.5)
%             A{c}(k, :) = struct_var(g).(field_rate_map_name)(c, :);
%             k = k + 1;
%         end
%     end
%
%     [A_n{c}, mat_divison] = fct_matnorm(A{c}, 1);
%     [ind_perm, ~] = find(A{c} == mat_divison);
%     A_np{c} = A_n{c}(ind_perm, :);
% end

function [A_np, ind_perm] = fct_rank_pf(rate_map_cell_bin, str_option)


% A = [];
% A_np = [];
% k = 1;
% for g = 1:size(rate_map_cell_bin, 1)
%     if any(boot_pf_cell_bin(g, :) == 1) && (stability_cell(g) >= 0.5)
%         A(k, :) = rate_map_cell_bin(g, :);
%         k = k + 1;
%     end
%
% end
A = rate_map_cell_bin;


switch str_option
    case 'cont' % Real Function
        
        if ~isempty(A) % Il arrive des cas où la matrice A est vide car il n'y a aucune cellule stable
            % A partir de cette étape, on enlève certains PF alors qu'on devrait pas
            % !!! Cela vient du fait qu'on mettait [ind_perm, ~] = find(abs(A ==
            % mat_divison); Le égal est dangereux sous Matlab !!
            [A_n, mat_divison] = fct_matnorm(A, 2);
            crap = abs(A - mat_divison);
            [ind_perm, c] = find(abs(A - mat_divison) < 5*eps);
            ind_perm = unique(ind_perm,'stable');
            %[ind_perm, c] = find(abs(A == mat_divison));
            A_np = A_n(ind_perm, :);
            %ind_perm = unique(ind_perm);
        end
        
        
    case 'bin' % Binary Function (works only for one sequence (one placefield))
        
        %       A = [0 1 1 0 0 ; 1 0 0 0 0 ; 1 1 0 0 0 ; 0 0 0 1 1 ; 0 0 1 1 1 ; 0 1 1 1 0 ];
        pf_size = cell2mat(arrayfun(@(k) sum(A(k, :) == 1), 1:size(A, 1), 'uni', 0));
        [~,pf_order_left] = max(A');
        pf_order_left_pond = pf_order_left*size(A, 2) + pf_size;
        [~, ind_perm] = sort(pf_order_left_pond);
        A_np = A(ind_perm, :);
end



