pixelSize = 1; %pixel size in physical units
sensorSizeX = 10; %given in pixels
sensorSizeY = 0;
thetaSpread = 10; %in degrees
phiSpread = 0;
nrays = 100;

%generate random (x,y) positions within each pixel
%index (0,0) is in bottom left corner of the sensor
raysPerPixel = round(nrays./sensorSizeX);
xr = [];
yr = [];
currentPixel = 1;
pixelIndex = [];

%should preallocate xr, yr, and pixelIndex for more efficiency
for i = 0:sensorSizeX - 1
%     for j = 0:sensorSizeY - 1
        xr = [xr; (rand(raysPerPixel,1) + i)*pixelSize];
%         yr = [yr; (rand(raysPerPixel,1) + j)*pixelSize];
        pixelIndex = [pixelIndex; zeros(raysPerPixel,1) + currentPixel];
        currentPixel = currentPixel + 1;
%     end
end

%generate random angles within the angular spread
%negative angle means below the horizontal
th = rand(raysPerPixel * sensorSizeX,1)*thetaSpread - thetaSpread./2;

%generate random phi angles within the spread
%negative angle means below the horizontal
ph = rand(raysPerPixel * sensorSizeX,1)*phiSpread - phiSpread./2;

%distance to diffuser
z = 100;

%propagate to the diffuser by a distance z
%angle at which diffuser is hit stays the same
xo = z * tand(th) + xr;

%looking at whole diffuser
in = load('../Output/diffuser.mat');

%diffuser strength
strength = 50;
diffuser_in = in.filtered * strength;

%coordinate system in physical units from the diffuser file
x = in.x;
diff_upsample = false;
if diff_upsample
    diffuser = imresize(diffuser_in,diff_upsample,'bicubic');
    x = linspace(min(x),max(x),numel(x)*diff_upsample);
else
    diffuser = diffuser_in;
end
y = x;
px = mean(diff(x)); %diffuser "pixel" size in um/pixel
%physical units

%Setup gradients
% if dy_idx == 0 %flat in the y-direction
    Fx = gradient(diffuser(1,:));
    Fy = zeros(size(Fx));
% elseif dx_idx == 0 %flat in the x-direction
%     Fy = gradient(diffuser);
%     Fx = zeros(size(Fy));
% else %sloped in both directions
%     [Fx, Fy] = gradient(diffuser);
%  end

%Calculate surface norm at random positions by
%interpolation within gradient
%--------------------------------------------------------------
% %Get z and surface normal at each random (xr(i),yr(i)) pair
% if dy_idx==0
    Fxr = interp1(x,Fx,xo);  %Interpolate x gradient
    Fyr = zeros(size(Fxr));
%     %zr = interp1(xg,diff_crop,xr);   %Interpolate surface
% elseif dx_idx==0
%     Fyr = interp1(yg,Fy_crop,yr);  %Interpolate x gradient
%     Fxr = zeros(size(Fyr));
%     %zr = interp1(yg,diff_crop,yr);   %Interpolate surface
% else
%     %zr = interp2(xg,yg,diff_crop,xr,yr);   %Interpolate surface
%     Fyr = interp2(xg,yg,Fx_crop',xr,yr);  %Interpolate x gradient
%     Fxr = interp2(xg,yg,Fy_crop',xr,yr);  %Interpolate y gradiet
% end
%-----------------------------------------------------------
%index of refraction
index = 1.5;

%refraction
[uxp, ~, ~] = refraction(Fxr, Fyr, th, 0, index);
%outputs the new angle after refraction in the x-direction

uxp = 90 - acosd(uxp);

%constructs cell array for output
output = cell(2,3);
output{1,1} = 'Pixel';
output{2,1} = pixelIndex;
output{1,2} = 'x';
output{2,2} = xo;
output{1,3} = 'theta';
output{2,3} = uxp;

%should not be getting negative output x-value because the
%sensor and the diffuser are the same size and the origin is in
%the bottom left corner of the sensor