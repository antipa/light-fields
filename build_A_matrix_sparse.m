function [r,c,v] = build_A_matrix_sparse(gatherer,LF_index)
[r,c,v] = find(gatherer(:));
c = c*LF_index;
end