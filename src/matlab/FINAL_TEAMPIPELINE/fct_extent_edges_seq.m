% Author(s): Marti Geoffrey   
% Epsztein Lab 2019

% Given a set of intervals 'seq' (matrix with row = #interval and col =
% start/stop index) included in a set of intervals 'seq_e' this 
% function will extend the intervals in 'seq' by 'n' points in each direction 
% if it overlaps with the intervals in 'seq_e'.

function seq_f = fct_extent_edges_seq(seq, seq_e, n)


if isempty(seq)
    seq_f = seq;
    return
end

N = max([seq_e(:) ; seq(:)]);

mtd = 1;

switch mtd
    case 1 % this version tests if seq is fully included in seq_e, otherwise extend seq_e to make it 
        vec = zeros(1, N);
        idx = fct_itv2idx(seq);
        vec(idx) = 1;
        
        vec_e = zeros(1, N);
        idx_e = fct_itv2idx(seq_e);
        vec_e(idx_e) = 1;
        
        vec_c = zeros(1, N);
        idx_c = fct_itv2idx([(seq(:, 1) - n) (seq(:, 2) + n)]);
        idx_c(idx_c < 1 | idx_c > N) = [];
        vec_c(idx_c) = 1;
        
        if sum(vec & vec_e) ~= sum(vec)
            warning('Some intervals of the matrix ''seq'' are not included in ''seq_e''')
            vec_e(vec == 1) = 1;
        end
        seq_f = fct_find_seq(vec_c & vec_e);
        
    case 2
        idx = fct_itv2idx([(seq(:, 1) - n) (seq(:, 2) + n)]);
        idx(idx < 1 | idx > N) = [];
        idx = unique(idx);
        idx_e = fct_itv2idx(seq_e);
        vec = zeros(1, N);
        vec(idx(ismember(idx, idx_e))) = 1;
        seq_f = fct_find_seq(vec);
end





