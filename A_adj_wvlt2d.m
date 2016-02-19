function out = A_adj_wvlt2d(x,At,wvlt,N,order_vec,ntheta,nphi,nx,ny,wvlt_size)
b = At*x;
bord = reshape(b(order_vec),ny,nx,ntheta*nphi);
out = zeros(ntheta*nphi*wvlt_size,1);
for m = 1:ntheta*nphi
    out((m-1)*wvlt_size+1:wvlt_size*m) = wavedec2(bord(:,:,m),N,wvlt)';
end
