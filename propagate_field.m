function [propagated, X, Y] = propagate_field(xi,eta,Ui,Z,lambda,varargin)
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

k = 2*pi/lambda;
[XI, ETA] = meshgrid(xi,eta);
kernel = exp(1i*k/2/Z*(XI.^2 + ETA.^2));
fx = linspace(-1/2/mean(diff(xi)),1/2/mean(diff(xi)),numel(xi));
fy = linspace(-1/2/mean(diff(eta)),1/2/mean(diff(eta)),numel(eta));
x = fx*lambda*Z;
y = fy*lambda*Z;
[X, Y] = meshgrid(x,y);
quad_phase = exp(1i*k*Z)./(1i*lambda*Z).*exp(1i*k/2/Z*(X.^2+Y.^2));

propagated = quad_phase.*fftshift(fft2(fftshift(kernel.*Ui.*phase_mask)));




