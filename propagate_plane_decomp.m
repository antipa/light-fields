function [propagated, varargout] = propagate_plane_decomp(xi,eta,Ui,Z,lambda,varargin)
%apply basic Fresnel propagation
%N. Antipa 9/2/2014 Berkeley-Waller group
%Inputs
%xi : input x-direction vector  size Nx1
%eta : input y-direction vector size Mx1
%Ui : input complex field 2d (MxN)
%Z : propagation distance in same units as xi and eta
%lambda : wavelength
%phase_mask : optional mask to apply prior to propagation
%
%Outputs
%propagated : output field
%X and Y are output x and y grids 3d plaid
if numel(varargin)==0
    phase_mask = ones(length(eta),length(xi));
else
    phase_mask = varargin{1};
end
pad = 1;

    
fx = linspace(-1/2/mean(diff(xi)),1/2/mean(diff(xi)),numel(xi)*pad);
fy = linspace(-1/2/mean(diff(eta)),1/2/mean(diff(eta)),numel(eta)*pad);
[Fx, Fy] = meshgrid(fx,fy);
pad_r1 = floor(.5*(pad*numel(eta)-size(Ui,1)));
pad_r2 = pad*numel(eta)-size(Ui,1)-pad_r1;
pad_c1 = floor(.5*(pad*numel(xi)-size(Ui,2)));
pad_c2 = pad*numel(xi)-size(Ui,2)-pad_c1;
win_pad = 100;
if pad==1;
    %apodize 
    a0 = .355768;
    a1 = .487396;
    a2 = .144232;
    a3 = .012604;
    xp = 1:size(Ui,2);
    xp = xp-mean(xp);
    yp = linspace(1,size(Ui,2),numel(eta));
    yp = yp-mean(yp);
    [Xp, Yp] = meshgrid(xp,yp);
    
    %Rp = sqrt(Xp.^2+Yp.^2);
    Rp = max(cat(3,abs(Xp),abs(Yp)),[],3);
    Rp = Rp-(numel(xi)/2-win_pad);
    Rp(Rp<0) = 0;
    N = 2*max(max(Rp));
    window_2d = a0-a1*cos(2*pi*(Rp-N/2)/(N-1))+a2*cos(4*pi*(Rp-N/2)/(N-1))-a3*cos(6*pi*(Rp-N/2)/(N-1));
    Ui = Ui.*window_2d;
end
Up = padarray(Ui,[pad_r1 pad_c1],0,'pre');
Up = padarray(Up,[pad_r2 pad_c2],0,'post');
Uf = fftshift(fft2(fftshift(Up)));
Hf = exp(1i*2*pi*Z/lambda * sqrt(1-(lambda*Fx).^2 - (lambda*Fy).^2));
Rf = sqrt(Fx.^2 + Fy.^2);
Hf(Rf>1/lambda) = 0;
Uf_prop = Hf.*Uf;
propagated = ifftshift(ifft2(ifftshift(Uf_prop)));
if nargout>1
    
    varargout{1} = Fx*lambda*Z;
    varargout{2} = Fy*lambda*Z;
end



