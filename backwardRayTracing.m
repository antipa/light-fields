try
    tracerays_settings;
catch exception
    display('Error loading settings. Read instructions in tracerays_settings.m.template.');
end

if sensorSizeY == 0 && phiSpread == 0
    twoD = true;
else
    twoD = false;
end

%generate random (x,y) positions within each pixel
if twoD
    raysPerPixel = round(rays./sensorSizeX);
else
    raysPerPixel = round(rays./(sensorSizeX .* sensorSizeY));
end

xr = [];
yr = [];
currentPixel = 1;
pixelIndex = [];

if twoD
    %should preallocate xr, yr, and pixelIndex for more efficiency
    for i = 0:sensorSizeX - 1
        xr = [xr; (rand(raysPerPixel,1) + i) * pixelSize];
        pixelIndex = [pixelIndex; zeros(raysPerPixel,1) + currentPixel];
        currentPixel = currentPixel + 1;
    end
else
    for i = 0:sensorSizeX - 1
        for j = 0:sensorSizeY - 1
            xr = [xr; (rand(raysPerPixel,1) + i)*pixelSize];
            yr = [yr; (rand(raysPerPixel,1) + j)*pixelSize];
            pixelIndex = [pixelIndex; zeros(raysPerPixel,1) + currentPixel];
            currentPixel = currentPixel + 1;
        end
    end
end

%generate random angles within the angular spread
%negative angle means below the horizontal
if twoD
    th = rand(raysPerPixel * sensorSizeX,1)*thetaSpread - thetaSpread./2;
    ph = zeros(raysPerPixel * sensorSizeX,1);
else
    %generate random phi angles within the spread
    %negative angle means below the horizontal
    th = rand(raysPerPixel * sensorSizeX * sensorSizeY,1)*thetaSpread - thetaSpread./2;
    ph = rand(raysPerPixel * sensorSizeY * sensorSizeX,1)*phiSpread - phiSpread./2;
end

%propagate to the diffuser by a distance z
%angle at which diffuser is hit stays the same
xo = z * tand(th) + xr;

if twoD
    yo = zeros(length(xr),1);
else
    yo = z * tand(ph) + yr;
end

%looking at whole diffuser
in = load('../Output/diffuser.mat');

diffuser_in = in.filtered * strengthB;

%coordinate system in physical units from the diffuser file
x = in.x;
%first pixel starts at 0
x = x - min(x);
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


if twoD
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
    %
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
else
    [Fx, Fy] = gradient(diffuser);
end


if twoD
    Fxr = interp1(x,Fx,xo); %Interpolate x gradient
    Fyr = zeros(size(Fxr));
else
    Fyr = interp2(x,y,Fx,xo,yo);
    Fxr = interp2(x,y,Fy,xo,yo);
end

good = ~isnan(Fxr) & ~isnan(Fyr);
Fxr = Fxr(good);
Fyr = Fyr(good);
th = th(good);
ph = ph(good);
xo = xo(good);
yo = yo(good);
pixelIndex = pixelIndex(good);


%refraction
[uxp, uyp,~] = refraction(Fxr, Fyr, th, ph, indexEnv, indexDiff);
%outputs the new angle after refraction in the x-direction
uxp = 90 - acosd(uxp);
uyp = 90 - acosd(uyp);

%constructs cell array for output
output = cell(2,5);
output{1,1} = 'Pixel';
output{2,1} = pixelIndex;
output{1,2} = 'x';
output{2,2} = xo;
output{1,3} = 'theta';
output{2,3} = uxp;
output{1,4} = 'y';
output{2,4} = yo;
output{1,5} = 'phi';
output{2,5} = uyp;

stepX = (max(xRange) - min(xRange)) ./ gridX;
stepT = (tRange(2) - tRange(1)) ./ gridT;
xValues = [];
for l = xRange(1):stepX:xRange(2)
    for q = 1:gridX
        xValues = [xValues l];
    end
end
yValues = [tRange(2):-stepT:tRange(1)];

figure(1);
hold on
grid on
aMatrix = zeros(sensorSizeX, gridX .* gridT);
for i = 1:sensorSizeX
    list = find(output{2,1} == i);
%     scatter(output{2,2}(list),output{2,3}(list) * indexDiff);
    for j = 1: gridX * gridT
        xmin = xValues(j);
        xmax = xValues(j+gridX);
        ymin = yValues(mod(j-1,length(yValues)-1) + 1);
        ymax = yValues(mod(j-1,length(yValues)-1) + 2);
        count = 0;
        a = output{2,2}(list) > xmin & output{2,2}(list) <= xmax & ...
            output{2,3}(list) > ymin & output{2,3}(list) <= ymax;
        b = output{2,2}(list) > xRange(1) & output{2,2}(list) <= xRange(2) & ...
            output{2,3}(list) > tRange(2) & output{2,3}(list) <= tRange(1);
        if sum(b)
        aMatrix(i,j) = sum(a) ./ sum(b);
        end
    end
end
hold off
