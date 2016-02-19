function out = A_wvlt2d(x,S,wvlt,A,nx,ny,ntheta,nphi,deindex_vec)
% wvlt: wavelet strings
% decsize: number of wavelet coefficients for each image
% A : forward matrix
% S : 2d wavelet bookkeeping matrix
% x : vector of wavelet coefficients
% deidex_vec: go from image stack to vector
xp = zeros(nx,ny,ntheta*nphi);
decsize = length(x)/ntheta/nphi;
for m = 1:ntheta*nphi
    vind = (m-1)*decsize+1:m*decsize;
    xp(:,:,m) = waverec2(x(vind),S,wvlt);
end
xm = xp(deindex_vec);
out = A*xm;