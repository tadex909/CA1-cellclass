% Author(s): Marti Geoffrey   
% Epsztein Lab 2019

% Based on a vector 'vec' of zeros and ones, this function will look for
% contiguous sequence of 1 and will return the beginning and the end of
% each sequence in a matrix seq_ind (row: sequence number, col:
% start/stop). 'seq_nb' is the total number of sequences and 'seq_length'
% is the length of each sequence.
% As a second input argument, you can put 'minlen' followed by the minimum
% sequence size you are looking for.
% Additionnaly, you can put 'onlyfull' to remove sequences connected to
% the edges of the vector.



function [seq_ind, seq_nb, seq_length] = fct_find_seq(vec, varargin)


mode = 2; % Detection method

if any(strcmp(varargin, 'minlen'))
    minlen = varargin{find(strcmp(varargin, 'minlen')) + 1};
else
    minlen = 1;
end


if any(strcmp(varargin, 'onlyfull'))
    onlyfull = true; % Include/exclude incomplete sequences
else
    onlyfull = false;  
end



if ~any(vec == 0)
    mode = 1;
end

switch mode
    case 1
        
        ind_ones = find(vec == 1);
        space_btw_ones = diff(ind_ones);
        
        
        
        sub_ind_ones_beg_end_seq = find(space_btw_ones > 1);
        sub_ind_ones_beg_end_seq = [0 sub_ind_ones_beg_end_seq length(ind_ones)];
        seq_length = diff(sub_ind_ones_beg_end_seq);
        
        if minlen > 1
            idel = find(seq_length < minlen);
            sub_ind_ones_beg_end_seq(idel) = [];
            seq_length(idel) = [];
        end
        
        
        if seq_length == 0
            seq_nb = 0;
        else
            seq_nb = length(seq_length);
        end
        
        %         if seq_nb > 0
        %            seq_ind = zeros(seq_nb, 2);
        %         else
        %            seq_ind = []; % Juste pour éviter Empty matrix: 0-by-2
        %         end
        
        seq_ind = zeros(seq_nb, 2);
        
        for p = 1:seq_nb
            seq_ind(p, 1) = ind_ones(sub_ind_ones_beg_end_seq(p) + 1);
            seq_ind(p, 2)  = ind_ones(sub_ind_ones_beg_end_seq(p) + seq_length(p));
        end
        
        if onlyfull
            if (seq_nb > 0) && (seq_ind(1, 1) == 1)
                % Exclude first seq if incomplete
                seq_ind(1, :) = [];
                seq_length(1) = [];
                seq_nb = seq_nb - 1;
            end
            
            if (seq_nb > 0) &&  (seq_ind(end, 2) == length(vec))
                % Exclude last seq if incomplete
                seq_ind(end, :) = [];
                seq_length(end) = [];
                seq_nb = seq_nb - 1;
            end
            
        end
        
    case 2
        
        istart_seq = find(diff(vec) > 0);
        istop_seq = find(diff(vec) < 0);
        
        if ~onlyfull
            % Add first seq if incomplete
            if length(istop_seq) == (length(istart_seq) + 1)
                istart_seq = [0 istart_seq];
            end
            
            % Add last seq if incomplete
            if (length(istop_seq) + 1) == length(istart_seq)
                istop_seq = [istop_seq length(vec)];
            end
            
            % Add first and last
            if ~isempty(istart_seq) && (istart_seq(1) > istop_seq(1))
                istart_seq = [0 istart_seq];
                istop_seq = [istop_seq length(vec)];
            end
        else
            % Exclude first seq if incomplete
            if length(istop_seq) == (length(istart_seq) + 1)
                istop_seq(1) = [];
            end
            
            % Exclude last seq if incomplete
            if (length(istop_seq) + 1) == length(istart_seq)
                istart_seq(end) = [];
            end
            
            % Exclude first and last seq if both incomplete
            if ~isempty(istart_seq) && (istart_seq(1) > istop_seq(1))
                istop_seq(1) = [];
                istart_seq(end) = [];
            end
        end
        
        
        if minlen > 1
            idel = find((istop_seq' -(istart_seq + 1)'  + 1) < minlen);
            istart_seq(idel) = [];
            istop_seq(idel) = [];
        end
        
        seq_ind = [(istart_seq + 1)' istop_seq'];
        seq_length = diff(seq_ind') + 1;
        seq_nb = length(seq_length);
        
end


%% Debug
% vec = [0 0 0 1 1 1 0 0];
% vec = [0 0 0 1 1 0 0 0 0 1 0 0 0 1 1 1 0 0 0 0 ];
% vec = [1 1];
% vec = [1 1 0 0];
% vec = [0 0 0 0 1 1 1 1 1];
% vec = [1 1 0 0 1 1 1 1 0 0 0 1 1 1 1 1 0 0 1 1 1];
% vec = [1 0 0 0 1 0 0 0 1];
