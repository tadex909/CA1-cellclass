function [seq_ind, seq_nb, seq_length] = fct_find_seq(vec, varargin)
%FCT_FIND_SEQ Robust sequence finder for logical/numeric vectors.
%
% This compatibility copy shadows the legacy fct_find_seq when theta_s.m
% runs. It normalizes input orientation first, avoiding row/column
% concatenation errors in the original implementation.

vec = vec(:).' ~= 0;

if any(strcmp(varargin, 'minlen'))
    minlen = varargin{find(strcmp(varargin, 'minlen'), 1, 'first') + 1};
else
    minlen = 1;
end

onlyfull = any(strcmp(varargin, 'onlyfull'));

if isempty(vec)
    seq_ind = zeros(0, 2);
    seq_nb = 0;
    seq_length = [];
    return;
end

starts = find(diff([false vec]) == 1);
stops = find(diff([vec false]) == -1);

if onlyfull
    keep = starts > 1 & stops < numel(vec);
    starts = starts(keep);
    stops = stops(keep);
end

seq_length = stops - starts + 1;
keep = seq_length >= minlen;

starts = starts(keep);
stops = stops(keep);
seq_length = seq_length(keep);

seq_ind = [starts(:) stops(:)];
seq_nb = numel(seq_length);
seq_length = seq_length(:).';
end
