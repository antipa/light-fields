%object dimensions in voxels
vx = 1; %in x-direction
vy = 1; %in y-direction
vz = 1; %in z-direction

%voxel dimensions in microns
voxX = 10;
voxY = 10;
voxZ = 100;

%total rays traced from entire object
rays = 10000000;

%angle spread of rays, in degrees
thetaSpread = 5;
phiSpread = 5;

%define number of planes in z-direction for each voxel
k = 1;

%assuming voxels are size 1
dz = voxZ/k;

%distance from closest plane of voxel to diffuser, in microns
z0 = 100;

%distance from diffuser to sensor, in microns
z1 = 200;

%sensor dimensions
npx = 50; %number of pixels in x-direction
npy = 50; %number of pixels in y-direction
x_range = 150; %total range of sensor in x-direction, in microns
y_range = 150; %total range of sensor in y-direction, in microns

%rays traced from each voxel
raysPerVoxel = round(rays./(vx*vy*vz));

%diffuser properties
strength = 0;
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

%initialize the measurement matrix
hMatrix = zeros(npx*npy,vx*vy*vz);

count = 0;
for a = 1:vx
    for b = 1:vy
        %repeat for each plane in z-direction
        for c = 1:vz
            gatherer = zeros(npy,npx);
            for j = 1:k
                %random x, y, theta, and phi column vectors
                xr = (rand(raysPerVoxel,1) + (a-1)) * voxX;
                yr = (rand(raysPerVoxel,1) + (b-1)) * voxY;
                th = (rand(raysPerVoxel,1) - 0.5) * thetaSpread;
                ph = (rand(raysPerVoxel,1) - 0.5) * phiSpread;
                
                %distance from z-plane to diffuser
                z = z0 + voxZ * (c-1) + dz * (j-1);
                
                %propagate rays to the diffuser
                xo = z * tand(th) + xr + 25;
                yo = z * tand(ph) + yr + 25;
                
                Fyr = interp2(x(1:ceil(max(xo))+5),y(1:ceil(max(yo))+5),Fy((1:ceil(max(yo))+5),(1:ceil(max(xo))+5)),xo,yo);
                Fxr = interp2(x(1:ceil(max(xo))+5),y(1:ceil(max(yo))+5),Fx((1:ceil(max(yo))+5),(1:ceil(max(xo))+5)),xo,yo);
                
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
                gatherer = gatherer + hist4(xo,yo,npx,npy,dpx,dpy,0,0,px);
                %imagesc(gatherer);
                %pause(1/24);
            end
            %adding a column of the H matrix
            hMatrix(:,(sub2ind([vy,vx,vz],b,a,c))) = gatherer(:);
            
            %count to display progress to user
            count = count + 1;
            [num2str(count) '/' num2str(vx*vy*vz)]
        end
    end
end