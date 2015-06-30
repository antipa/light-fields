function A_shifted = shift_A_matrix(A_sub1,dx,dy,npx,npy)
%A_shifted = shift_A_matrix(A_sub1,dx,dy,npx,npy)
%
%Shift A matrix projection by dx,dy (positive dx shifts left, positive dy
%shifts down). npx and npy are sensor size in x and y respectively.
dx_px = dx;
dy_px = dy;

cstart = max(dx_px*npy+1+dy_px,1);
rstart = max(-dx_px*npy+1-dy_px,1);
b = mod(rstart-1:npx*npy-cstart,npx)<(npy-dy_px);
a = mod(rstart-1:npx*npy-cstart,npx)>-dy_px;
spvals = a.*b;
shifter = sparse(cstart:npx*npy-rstart+1,rstart:npx*npy-cstart+1,spvals,npx*npy,npx*npy);
%shifter = sparse(rstart:npx*npy-cstart+1,cstart:npx*npy-rstart+1,spvals,npx*npy,npx*npy);
A_shifted = shifter*A_sub1;

return