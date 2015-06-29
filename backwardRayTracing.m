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

%looking at whole diffuser
in = load('../Output/diffuser.mat');
diffuser_in = in.filtered * strengthB;

%coordinate system in physical units from the diffuser file
x = in.x;
%first pixel starts at 0 and goes to 6000
x = x - min(x);
diff_upsample = false;
if diff_upsample
    diffuser = imresize(diffuser_in,diff_upsample,'bicubic');
    x = linspace(min(x),max(x),numel(x)*diff_upsample);
else
    diffuser = diffuser_in;
end
y = x;
px = mean(diff(x)); %diffuser "pixel" size in um/pixel, physical units

%generate random (x,y) positions within each pixel
%setup gradients using first row of the diffuser
if twoD
    raysPerPixel = round(rays./sensorSizeX);
    Fx = gradient(diffuser(1,:));
    Fy = zeros(size(Fx));
else
    raysPerPixel = round(rays./(sensorSizeX * sensorSizeY));
    [Fx, Fy] = gradient(diffuser);
end

%constructing A matrix by dividing the diffuser into a grid
stepX = (max(xRange) - min(xRange)) / gridX;
stepT = (max(tRange) - min(tRange)) / gridT;

if twoD
    xValues = [];
    for l = min(xRange):stepX:max(xRange)
        for q = 1:gridT
            xValues = [xValues l];
        end
    end
    
    %start with the bottom left corner
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
    %construct a 4d matrix
    gridY = 3;
    yRange = [0 10];
    stepY = (max(yRange) - min(yRange)) / gridY;
    gridP = 2;
    pRange = [-2 2];
    stepP = (max(pRange) - min(pRange)) / gridP;
    %start with the bottom left corner
    pValues = min(pRange):stepP:max(pRange);
    lightField = zeros(sensorSizeX*sensorSizeY,gridX*gridY*gridT*gridP);
    for i = 0:sensorSizeX - 1
        for j = 0:sensorSizeY - 1
            xr = (rand(raysPerPixel,1) + i)*pixelSize;
            yr = (rand(raysPerPixel,1) + j)*pixelSize;
            %generate random phi angles within the spread
            %negative angle means below the horizontal
            th = rand(raysPerPixel,1)*thetaSpread - thetaSpread./2;
            ph = rand(raysPerPixel,1)*phiSpread - phiSpread./2;
            %propagate to the diffuser by a distance z
            %angle at which diffuser is hit stays the same
            xo = z * tand(th) + xr;
            yo = z * tand(ph) + yr;
            Fyr = interp2(x,y,Fx,xo,yo);
            Fxr = interp2(x,y,Fy,xo,yo);
            %throwing out points that did not hit the diffuser and could
            %not be interpolated
            good = ~isnan(Fxr) & ~isnan(Fyr);
            Fxr = Fxr(good);
            Fyr = Fyr(good);
            th = th(good);
            ph = ph(good);
            xo = xo(good);
            yo = yo(good);
            %refraction
            [uxp, uyp,~] = refraction(Fxr, Fyr, th, ph, indexEnv, indexDiff);
            %outputs the new angle after refraction in the x-direction
            uxp = 90 - acosd(uxp);
            uyp = 90 - acosd(uyp);
            tic
            for a = 1:gridX
                xmin = a * stepX;
                xmax = (a+1) * stepX;
                for b = 1:gridY
                    ymin = b * stepY;
                    ymax = (b+1) * stepY;
                    for c = 1:gridT
                        tmin = c * stepT;
                        tmax = (c+1) * stepT;
                        for d = 1:gridP
                            pmin = d * stepP;
                            pmax = (d+1) * stepP;
                            inBox = xo > xmin & xo <= xmax & uxp > tmin & uxp <= tmax & yo > ymin ...
                                & yo <= ymax & uyp > pmin & uyp <= pmax;
                            lightField(sub2ind(size(zeros(sensorSizeX,sensorSizeY)),i+1,j+1),sub2ind(size(zeros(gridX,gridY,gridT,gridP)),a,b,c,d)) = sum(inBox);
                        end
                    end
                end
            end
             toc
        end
    end
end