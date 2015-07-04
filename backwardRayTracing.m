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

%using whole diffuser
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
   
    stepY = (max(yRange) - min(yRange)) / gridY;

    stepP = (max(pRange) - min(pRange)) / gridP;
    
    % Collect the histogram of rays per sensor-pixel * diffuser bin
    % 1st dim = sensor pixel
    % 2nd dim = diffuser bin
    % value   = #rays
    
    %     rng(0) % SEED the random number generator.
    
    gatherer = zeros(sensorSizeX*sensorSizeY,gridX*gridY*gridT*gridP);
    lightField = zeros(sensorSizeX*sensorSizeY,gridX*gridY*gridT*gridP);
    for i = 0:sensorSizeX - 1
        for j = 0:sensorSizeY - 1
            % xr: array of "raysPerPixel" random value for x value of ray
            xr = (rand(raysPerPixel,1) + i)*pixelSize;
            yr = (rand(raysPerPixel,1) + j)*pixelSize;
            
            %generate random phi angles within
            %   the (- 0.5 * spread, + 0.5 * spread)
            %negative angle means below the horizontal
            th = (rand(raysPerPixel,1) - 0.5) * thetaSpread;
            ph = (rand(raysPerPixel,1) - 0.5) * phiSpread;
            
            %propagate to the diffuser by a distance z
            %angle at which diffuser is hit stays the same
            % xo, yo : position at the diffuser
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
            %refraction
            [uxp, uyp, ~] = refraction(Fxr, Fyr, th, ph, indexEnv, indexDiff);
            %outputs the new angle after refraction in the x-direction
            uxp = 90 - acosd(uxp);
            uyp = 90 - acosd(uyp);
            
            %%% gatherer computation starts here
            % x bins:
            %   min(xRange) is the x-center of the first bin
            %   gridX is the number of bins on x axis, stepX is the size of
            %   a bin on x axis
            % A bin includes the right-end but does not include the
            %   left-end
            %
            xMinFirst = min(xRange) - 0.5 * stepX;
            yMinFirst = min(yRange) - 0.5 * stepY;
            tMinFirst = min(tRange) - 0.5 * stepT;
            pMinFirst = min(pRange) - 0.5 * stepP;
            
            %bin numbers for each ray
            xo_r_sub = ceil((xo - xMinFirst) / stepX);
            yo_r_sub = ceil((yo - yMinFirst) / stepY);
            uxp_r_sub = ceil((uxp - tMinFirst) / stepT);
            uyp_r_sub = ceil((uyp - pMinFirst) / stepP);
            
            %check to see if the ray is in the grid
            good = xo_r_sub>0&yo_r_sub>0&xo_r_sub<=gridX&yo_r_sub<=gridY...
                & uxp_r_sub>0 &uxp_r_sub <= gridT & uyp_r_sub >0 & uyp_r_sub <= gridP;
            xo_r_sub = xo_r_sub(good);
            yo_r_sub = yo_r_sub(good);
            uxp_r_sub = uxp_r_sub(good);
            uyp_r_sub = uyp_r_sub(good);
            
            %gives each ray a 4D bin index
            ray_ind = sort(sub2ind([gridX,gridY,gridT,gridP],xo_r_sub,yo_r_sub,uxp_r_sub,uyp_r_sub));
            %gives index of each bin that contains at least 1 ray
            [ray_unique, ray_entry] = unique(ray_ind);
            
            %if rays hit the grid
            if nnz(good)
                if length(ray_unique) > 1
                    %number of rays in each bin
                    ray_hist = diff(ray_entry);
                    %number of rays in last bin
                    ray_hist(end+1) = length(ray_ind)-ray_entry(end)+1;
                    gatherer(sub2ind([sensorSizeY,sensorSizeX],j+1,i+1),ray_unique) = ray_hist;
                elseif length(ray_unique) == 1
                    %if 1 bin want total number of rays that hit the grid
                    %or 1 ray that hit the grid
                    gatherer(sub2ind([sensorSizeY,sensorSizeX],j+1,i+1),ray_unique) = sum(good);
                end
            end
            %if rays didn't hit the grid, leave the zeros
            
            % lightfield computation
            % x bins:
            %   min(xRange) is the x-center of the first bin
            %   gridX is the number of bins on x axis, stepX is the size of
            %   a bin on x axis
            % A bin includes the right-end but does not include the
            %   left-end
            for a = 1:gridX
                xmin = xMinFirst + (a - 1) * stepX; % x min for this bin
                xmax = xmin + stepX;
                for b = 1:gridY
                    ymin = yMinFirst + (b - 1) * stepY;
                    ymax = ymin + stepY;
                    for c = 1:gridT
                        tmin = tMinFirst + (c - 1) * stepT;
                        tmax = tmin + stepT;
                        for d = 1:gridP
                            pmin = pMinFirst + (d - 1) * stepP;
                            pmax = pmin + stepP;
                            
                            inBox = xo > xmin & xo <= xmax & ...
                                uxp > tmin & uxp <= tmax & ...
                                yo > ymin & yo <= ymax & ...
                                uyp > pmin & uyp <= pmax;
                            lightField(sub2ind([sensorSizeY,sensorSizeX],j+1,i+1),...
                                sub2ind([gridX,gridY,gridT,gridP],a,b,c,d)) = sum(inBox);
                        end
                    end
                end
            end
        end
    end
end