profile on;
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

in = load('../Output/diffuser.mat');
diffuser_in = in.filtered * strengthB;

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

if twoD
    raysPerPixel = round(rays./sensorSizeX);
    %2D gradient generated using first row of the diffuser
    Fx = gradient(diffuser(1,:));
    Fy = zeros(size(Fx));
else
    raysPerPixel = round(rays./(sensorSizeX * sensorSizeY));
    [Fx, Fy] = gradient(diffuser);
end

%size of each bin
stepX = (max(xRange) - min(xRange)) / gridX;
stepT = (max(tRange) - min(tRange)) / gridT;

if twoD
    xValues = [];
    for l = min(xRange):stepX:max(xRange)
        for q = 1:gridT
            xValues = [xValues l];
        end
    end
    
    tValues = min(tRange):stepT:max(tRange);
    aMatrix = zeros(sensorSizeX,(gridX .* gridT));
    %pixel indexing goes from bottom to top, left to right
    figure(1);
    hold on
    grid on
    for i = 0:sensorSizeX - 1
        xr = (rand(raysPerPixel,1) + i) * pixelSize;
        th = rand(raysPerPixel,1)*thetaSpread - thetaSpread./2;
        ph = zeros(raysPerPixel,1);
        %propagate to the diffuser by a distance z
        %angle at which diffuser is hit stays the same
        xo = z * tand(th) + xr;
        yo = zeros(length(xr),1);
        Fxr = interp1(x,Fx,xo); %Interpolate x gradient
        Fyr = zeros(size(Fxr));
        %throwing out points that did not hit the diffuser and could
        %not be interpolated
        good = ~isnan(Fxr);
        Fxr = Fxr(good);
        Fyr = Fyr(good);
        th = th(good);
        ph = ph(good);
        xo = xo(good);
        yo = yo(good);
        %refraction
        [uxp,~,~] = refraction(Fxr, Fyr, th, ph, indexEnv, indexDiff);
        %outputs the new angle after refraction in the x-direction
        uxp = 90 - acosd(uxp);
        scatter(xo,uxp * indexDiff);
        %constructing a row of the A matrix
        for j = 1: gridX * gridT
            xmin = xValues(j);
            xmax = xValues(j+gridT);
            tmin = tValues(mod(j-1,length(tValues)-1) + 1);
            tmax = tValues(mod(j-1,length(tValues)-1) + 2);
            a = xo > xmin & xo <= xmax & uxp > tmin & uxp <= tmax;
            aMatrix(i+1,j) = sum(a);
        end
    end
    aMatrix = aMatrix ./ repmat(sum(aMatrix),[size(aMatrix,1),1]);
    %hold off;
else
    %size of each bin
    stepY = (max(yRange) - min(yRange)) / gridY;
    stepP = (max(pRange) - min(pRange)) / gridP;
    
    %rng(0) %for testing, seed the random number generator
    
    %left most edge of the grid (not included)
    xMinFirst = min(xRange);
    yMinFirst = min(yRange);
    tMinFirst = min(tRange);
    pMinFirst = min(pRange);
    
    % Collect the histogram of rays per sensor-pixel * diffuser bin
    % 1st dim = sensor pixel
    % 2nd dim = diffuser bin
    % value   = # rays
    %   tic
    if useEfficient
        gatherer = zeros(sensorSizeX*sensorSizeY, gridX*gridY*gridT*gridP);
    end
    if useSimple
        lightField = zeros(sensorSizeX*sensorSizeY, gridX*gridY*gridT*gridP);
    end
    for i = 0:sensorSizeX - 1
        for j = 0:sensorSizeY - 1
            %random x and y positions on the sensor
            xr = (rand(raysPerPixel,1) + i)*pixelSize;
            yr = (rand(raysPerPixel,1) + j)*pixelSize;
            
            %generate random theta and phi angles within
            %the (- 0.5 * spread, + 0.5 * spread)
            %negative angle means below the horizontal
            th = (rand(raysPerPixel,1) - 0.5) * thetaSpread;
            ph = (rand(raysPerPixel,1) - 0.5) * phiSpread;
            
            %x and y positions on the diffuser as a result of propagating
            %by distance z
            xo = z * tand(th) + xr;
            yo = z * tand(ph) + yr;
            
%             good = xo > xMinFirst & xo <= 6001 & yo > yMinFirst & yo <= 6001;
%             xo = xo(good);
%             yo = yo(good);
%             th = th(good);
%             ph = ph(good);
            
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
            
            %outputs new theta and phi angles after refraction
            [uxp, uyp, uzp] = refraction(Fxr, Fyr, th, ph, indexEnv, indexDiff,'angles');
            
            %propagate a second distance from the diffuser to the light
            %source
            [yo, xo] = propagation(uyp, uzp, z0, yo, uxp, xo);
            
            
            if useEfficient
                gatherer = gathererEfficient(gatherer,gridX,gridY,gridT,gridP,xMinFirst,yMinFirst,tMinFirst,...
                    pMinFirst,stepX,stepY,stepT,stepP,i,j,sensorSizeX,sensorSizeY,xo,yo,uxp,uyp);
            end
            if useSimple
                lightField = gathererSimple(lightField,gridX,gridY,gridT,gridP,xMinFirst,yMinFirst,tMinFirst,...
                    pMinFirst,stepX,stepY,stepT,stepP,i,j,sensorSizeX,sensorSizeY,xo,yo,uxp,uyp);
            end
        end
    end
    
    %for testing, check to see if outputs from the gatherer functions agree
    if useSimple && useEfficient
        assert(isequal(lightField,gatherer));
        fprintf('Total rays in gatherer: %d\n', sum(sum(gatherer)));
    end
end

%fprintf('%d\n',toc);
profile viewer;