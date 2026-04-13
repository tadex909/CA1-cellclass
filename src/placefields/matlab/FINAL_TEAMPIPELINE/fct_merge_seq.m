% Author(s): Marti Geoffrey   
% Epsztein Lab 2019

% This function receives as input a matrix named 'seq' of the indexes of
% start/stop (col) of sequences of contiguous 1 (row = #sequence) and will merge sequences
% separated by less than 'nb_points' fo finally obtain the matrix 'merged_seq'.

function [merged_seq, nb_seq, seq_length] = fct_merge_seq(seq, nb_points)



if isempty(seq)
    merged_seq = [];
    nb_seq = 0;
    seq_length = [];
    return
end


p = 0;
for k = 1:size(seq, 1)
    if (k > 1) && (seq(k, 1) - seq(k-1, 2) <= nb_points)
        merged_seq(p, 2) = seq(k, 2);
    else
        p = p + 1;
        merged_seq(p, :) = seq(k, :);
    end
end
nb_seq = p;
seq_length = diff(merged_seq') + 1;

% Second method
%  isinfpt = seq(2:end, 1) - seq(1:(end-1), 2)< nb_points;
%  isinfpt = isinfpt';
%  nb_merge = length(find(isinfpt == true));
% [imerg, nb_submerge] = fct_sequence_detect(isinfpt);
%
% ind_rem = [];
% for p = 1:nb_submerge
% seq(imerg(p, 1), 2) = seq(imerg(p, 2) + 1, 2);
% ind_rem = [ind_rem (imerg(p, 1) + 1):(imerg(p, 2) + 1)];
% end
% seq(ind_rem, :) = [];
