function lightField = gathererSimple(lightField,gridX,gridY,gridT,gridP,xMinFirst,...
    yMinFirst,tMinFirst,pMinFirst,stepX,stepY,stepT,stepP,i,j,sensorSizeX,sensorSizeY,xo,yo,uxp,uyp)
%This function loops over every bin on the diffuser and counts how many
%rays from a single sensor pixel hit each bin. See gathererEfficient.m for
%a more efficient implementation.

%bins include the right end but do not include the left end
for a = 1:gridX
    xmin = xMinFirst + (a - 1) * stepX; % left end of the bin
    xmax = xmin + stepX; %right end of the bin
    for b = 1:gridY
        ymin = yMinFirst + (b - 1) * stepY;
        ymax = ymin + stepY;
        for c = 1:gridT
            tmin = tMinFirst + (c - 1) * stepT;
            tmax = tmin + stepT;
            for d = 1:gridP
                pmin = pMinFirst + (d - 1) * stepP;
                pmax = pmin + stepP;
                
                %checks to see which rays are in the bin you are
                %considering
                inBox = xo > xmin & xo <= xmax & ...
                    uxp > tmin & uxp <= tmax & ...
                    yo > ymin & yo <= ymax & ...
                    uyp > pmin & uyp <= pmax;
                
                %assigns the number of rays in the bin to the sensor pixel
                %the rays came from and the linear index of the bin
                lightField(sub2ind([sensorSizeY,sensorSizeX],j+1,i+1),...
                    sub2ind([gridX,gridY,gridT,gridP],a,b,c,d)) = sum(inBox);
            end
        end
    end
end
end