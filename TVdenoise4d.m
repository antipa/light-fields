function u = TVdenoise4d(f,lambda,iters,ng,angweight,apclip)
%TVDENOISE  Total variation grayscale and color image denoising
% Updated for 4D light fields: 
%   u = TVdenoise4d(f,lambda,iters,ng,angweight,apclip)
%
%Inputs:   
%   f - light field to denoise
%
%   lambda - data weight (i.e. lower this to regularize more)
%
%   iters - passes in denoising (5-8 seems good)
%
%   ng - nonnegativeity constraint (boolean)
%   
%   angweight - weighting for angular dimensions. Defaults is 1
%   
%   apclip - array of same size as f that contains NaN where we don't want
%   to count TV (i.e. transition from inside to outside aperture shouldn't
%   count). Pass in [] to do nothing.
%
%
%   u = TVDENOISE(f,lambda) denoises the input image f.  The smaller
%   the parameter lambda, the stronger the denoising.
%
%   The output u approximately minimizes the Rudin-Osher-Fatemi (ROF)
%   denoising model
%
%       Min  TV(u) + lambda/2 || f - u ||^2_2,
%        u
%
%   where TV(u) is the total variation of u.  If f is a color image (or any
%   array where size(f,3) > 1), the vectorial TV model is used,
%
%       Min  VTV(u) + lambda/2 || f - u ||^2_2.
%        u
%
%   TVDENOISE(...,Tol) specifies the stopping tolerance (default 1e-2).
%
%   The minimization is solved using Chambolle's method,
%      A. Chambolle, "An Algorithm for Total Variation Minimization and
%      Applications," J. Math. Imaging and Vision 20 (1-2): 89-97, 2004.
%   When f is a color image, the minimization is solved by a generalization
%   of Chambolle's method,
%      X. Bresson and T.F. Chan,  "Fast Minimization of the Vectorial Total
%      Variation Norm and Applications to Color Image Processing", UCLA CAM
%      Report 07-25.
%
%   Example:
%   f = double(imread('barbara-color.png'))/255;
%   f = f + randn(size(f))*16/255;
%   u = tvdenoise(f,12);
%   subplot(1,2,1); imshow(f); title Input
%   subplot(1,2,2); imshow(u); title Denoised

% Pascal Getreuer 2007-2008
%  Modified by Jose Bioucas-Dias  & Mario Figueiredo 2010
%  (stopping rule: iters)
%
%
% Last modified by Nick Antipa 5/19/2016
% 
% Added apclip function with NaN outside light field aperture. This
% prevents denoising based on strong transition from aperture to zeros.
% dc: if dc is included in model, don't include in denoising.
% siz: lf size vector [ntheta,nphi,nx,ny]
if lambda < 0
    error('Parameter lambda must be nonnegative.');
end
if isempty(apclip)
    apclip = ones(size(f));
end


tau = .25/2;

N = size(f);
id = [2:N(1),N(1)];
iu = [1,1:N(1)-1];
ir = [2:N(2),N(2)];
il = [1,1:N(2)-1];
ib = [2:N(3),N(3)];
ifr = [1,1:N(3)-1];
i4f = [2:N(4),N(4)];
i4r = [1,1:N(4)-1];


p1 = zeros(size(f));
p2 = zeros(size(f));
p3 = zeros(size(f));
p4 = zeros(size(f));

divp = zeros(size(f));
lastdivp = ones(size(f));

for i=1:iters
    lastdivp = divp;
    
    z = divp - f*lambda;
    z = z.*apclip;
    z1 = z(:,ir,:,:) - z;
    z1(isnan(z1)) = 0;
    z2 = z(id,:,:,:) - z;
    z2(isnan(z2)) = 0;
    z3 = z(:,:,ib,:) - z;
    z3(isnan(z3)) = 0;
    z4 = z(:,:,:,i4f) - z;
    z4(isnan(z4)) = 0;
    
    denom = 1 + tau*sqrt((angweight*z1).^2 + (angweight*z2).^2 + z3.^2 + z4.^2);
    
    p1 = (p1 + tau*z1)./denom;
    p2 = (p2 + tau*z2)./denom;
    p3 = (p3 + tau*z3)./denom;
    p4 = (p4 + tau*z4)./denom;


    divp = angweight*(p1 - p1(:,il,:,:)) + angweight*(p2 - p2(iu,:,:,:))+...
         p3 - p3(:,:,ifr,:) + p4 - p4(:,:,:,i4r); % divergence
end

u = f - divp/lambda;

if ng
    u = max(0,u);
end

% threeslice(u,88);

end