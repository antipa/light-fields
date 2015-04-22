tic
cvx_begin
variable lf_cvx(numel(lf))
minimize(norm(A_sub1*lf_cvx-intens_noisy,2)+.005*norm(lf_cvx,2))
subject to
lf_cvx >= 0
cvx_end
toc
%%
lf_cvx_reshaped = reshape(lf_cvx,[NPhi,NX]);
subplot(2,1,1)
imagesc(lf_cvx_reshaped)