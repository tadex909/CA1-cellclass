% Author(s): Marti Geoffrey   
% Epsztein Lab 2019

% If you want to remove contiguous sequences of 1 in the matrix seq (col:
% start/stop of sequence, row = #sequence)
% - with a length lower than N1, you can put in arguments: 'minlen' followed by N1
% - with a length greater than N2, you can put in arguments: 'maxlen' followed by N2

function [seq_rem, seq_nb, seq_length] = fct_remove_seq(seq, varargin)



if isempty(seq)
    seq_rem = [];
    seq_nb = 0;
    seq_length = 0;
    return
end

if any(strcmp(varargin, 'minlen'))
    minlen = varargin{find(strcmp(varargin, 'minlen')) + 1};
else
    minlen = 0;
end

if any(strcmp(varargin, 'maxlen'))
    maxlen = varargin{find(strcmp(varargin, 'maxlen')) + 1};
else
    maxlen = max(seq(:));
end



seq_rem = seq;
diffseq = (seq(:, 2) - (seq(:, 1) + 1)  + 1);
idel = (diffseq < minlen) | (diffseq > maxlen);
seq_rem(idel, :) = [];
seq_length = diff(seq_rem') + 1;
seq_nb = length(seq_length);



