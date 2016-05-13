function A = make_A_circular(A,P,Q)
% A_out = make_A_circular(A,P,Q)
% A: input matrix
% P: number of theta samples
% Q: number of phi samples
x = linspace(-1,1,P);
[X,Y] = meshgrid(x,x);
R = sqrt(X.^2+Y.^2);
aper = R<=1.05;
good_idx = find(aper(:));

for n = 1:P*Q
    if ~ismember(n,good_idx)
        A(:,n:P*Q:end)=sparse(0);
    end
end