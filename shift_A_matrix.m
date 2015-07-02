function [A_shifted, shifter]= shift_A_matrix(A_sub1,dx,dy,npx,npy)
%A_shifted = shift_A_matrix(A_sub1,dx,dy,npx,npy)
%
%Shift A matrix projection by dx,dy (positive dx shifts left, positive dy
%shifts down). npx and npy are sensor size in x and y respectively.
dx_px = floor(dx);
dy_px = floor(dy);

cstart = max(dx_px*npy+1+dy_px,1);
rstart = max(-dx_px*npy+1-dy_px,1);
cstartL = floor(cstart);
rstartL = floor(rstart);
b = mod(rstartL-1:npx*npy-cstartL,npx)<(npy-dy_px);
a = mod(rstartL-1:npx*npy-cstartL,npx)>-dy_px;
spvalsL = a.*b;
shifter = sparse(cstartL:npx*npy-rstartL+1,rstartL:npx*npy-cstartL+1,spvalsL,npx*npy,npx*npy);
%shifter = sparse(rstart:npx*npy-cstart+1,cstart:npx*npy-rstart+1,spvals,npx*npy,npx*npy);

%A is now shifted by dx_px and dy_px. Still need to check for and apply
%sub-pixel shifts if necessary.

residual_x = abs(dx-dx_px);
residual_y = abs(dy-dy_px);
if residual_x
    [~, shifter_x] = shift_A_matrix(A_sub1,dx_px+1,dy_px,npx,npy);
    shifter = shifter*(1-residual_x)+residual_x*shifter_x;
    A_shifted = shifter*A_sub1;
end

if residual_y
    [~, shifter_y] = shift_A_matrix(A_sub1,dx,dy_px+1,npx,npy);
    shifter = shifter*(1-residual_y)+residual_y*shifter_y;
end

A_shifted = shifter*A_sub1;


return