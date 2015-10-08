function A_shifted = shift_A_matrix(A_sub1,dx,dy,npx,npy)
%A_shifted = shift_A_matrix(A_sub1,dx,dy,npx,npy)
%
%Shift A matrix projection by dx,dy (positive dx shifts left, positive dy
%shifts down). npx and npy are sensor size in x and y respectively.
dx_px = floor(dx);
dy_px = floor(dy);





residual_x = abs(dx-dx_px);
residual_y = abs(dy-dy_px);

if ~residual_x && ~residual_y
    shifter = make_shifter(dx_px,dy_px,npx,npy);
    A_shifted = shifter*A_sub1;
    return
elseif residual_x && ~residual_y
    %Only x is subpixel
    shifter_x1 = make_shifter(dx_px,dy_px,npx,npy);
    shifter_x2 = make_shifter(dx_px+1,dy_px,npx,npy);
    shifter = shifter_x1*(1-residual_x) + shifter_x2*residual_x;
    A_shifted = shifter*A_sub1;
    return
elseif ~residual_x && residual_y
    %Only y is subpixel
    shifter_y1 = make_shifter(dx_px,dy_px,npx,npy);
    shifter_y2 = make_shifter(dx_px,dy_px+1,npx,npy);
    shifter = shifter_y1*(1-residual_y) + shifter_y2*residual_y;
    A_shifted = shifter*A_sub1;
    return
elseif residual_x && residual_y
    %Both are subpixel
    %Shift X first
    shifter_x1 = make_shifter(dx_px,0,npx,npy);
    shifter_x2 = make_shifter(dx_px+1,0,npx,npy);
    shifter_x = shifter_x1*(1-residual_x) + shifter_x2*residual_x;
    A_shiftedx = shifter_x*A_sub1;
    shifter_y1 = make_shifter(0,dy_px,npx,npy);
    shifter_y2 = make_shifter(0,dy_px+1,npx,npy);
    shifter_y = shifter_y1*(1-residual_y) + shifter_y2*residual_y;
    A_shifted = shifter_y*A_shiftedx;
    return
end

end

function shifter = make_shifter(dx_px,dy_px,npx,npy)

    cstart = max(dx_px*npx+1+dy_px,1);
    rstart = max(-dx_px*npx+1-dy_px,1);

    b = mod(rstart-1:npx*npy-cstart,npx)<(npx-dy_px);
    a = mod(rstart-1:npx*npy-cstart,npx)>=-dy_px;
    spvals = a.*b;
    shifter = sparse(cstart:npx*npy-rstart+1,rstart:npx*npy-cstart+1,spvals,npx*npy,npx*npy);
    
    
end