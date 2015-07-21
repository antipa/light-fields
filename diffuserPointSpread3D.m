%voxel dimensions
vx = 2; %in x-direction
vy = 4; %in y-direction
vz = 6; %in z-direction

%total rays traced from all z-planes of voxel
rays = 1000;

%angle spread of rays, in degrees
thetaSpread = 2;
phiSpread = 2;

%define number of planes in z-direction
k = 5;
dz = vz/k;

%distance from closest plane of voxel to diffuser
z0 = 100;

%distance from diffuser to sensor
z1 = 150;

%sensor dimensions
npx = 5; %number of pixels in x-direction
npy = 5; %number of pixels in y-direction
x_range = 5; %total range of sensor in x-direction, in microns
y_range = 5; %total range of sensor in y-direction, in microns

%rays traced from each plane
raysPerPixel = round(rays./k);

%diffuser properties
strength = 50;
indexEnv = 1;
indexDiff = 1.5;
diff_upsample = false;
in = load('../Output/diffuser.mat');
diffuser_in = in.filtered * strength;

%coordinate system in physical units from the diffuser file
x = in.x;

%first diffuser pixel starts at 0 and goes to 6000
x = x - min(x);

if diff_upsample
    diffuser = imresize(diffuser_in,diff_upsample,'bicubic');
    x = linspace(min(x),max(x),numel(x)*diff_upsample);
else
    diffuser = diffuser_in;
end

y = x;
px = mean(diff(x)); %diffuser "pixel" size in um/pixel, physical units

[Fx, Fy] = gradient(diffuser);

%initialize the histogram
gatherer = zeros(npy,npx);

%repeat for each plane in z-direction
for j = 1:k
    
%random x, y, theta, and phi column vectors   
xr = rand(raysPerPixel,1) * vx;
yr = rand(raysPerPixel,1) * vy;
th = (rand(raysPerPixel,1) - 0.5) * thetaSpread;
ph = (rand(raysPerPixel,1) - 0.5) * phiSpread;
 
%distance from z-plane to diffuser
z = z0 + dz * (j-1);

%propagate rays to the diffuser
xo = z * tand(th) + xr;
yo = z * tand(ph) + yr;

Fyr = interp2(x,y,Fy,xo,yo);
Fxr = interp2(x,y,Fx,xo,yo);

%throwing out points that did not hit the diffuser and could
%not be interpolated
good = ~isnan(Fxr) & ~isnan(Fyr);
Fxr = Fxr(good);
Fyr = Fyr(good);
th = th(good);
ph = ph(good);
xo = xo(good);
yo = yo(good);
            
%outputs direction cosines after refraction
[uxp, uyp, uzp] = refraction(Fxr, Fyr, th, ph, indexDiff, indexEnv,'angles');

%propagate from diffuser to the sensor
[yo, xo] = propagation(uyp, uzp, z1, yo, uxp, xo);

%create image at the sensor using 2D histogramming
dpx = x_range/npx;
dpy = y_range/npy;
% clf
% figure(1)
% imagesc(gather_rays_nohist(xo,yo,npx,npy,dpx,dpy,0,0,px));
gatherer = gatherer + gather_rays_nohist(xo,yo,npx,npy,dpx,dpy,0,0,px);
end

imagesc(gatherer);