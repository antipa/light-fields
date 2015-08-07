function lightField = gathererSimple(lightField,gridX,gridY,gridT,gridP,xMinFirst,...
    yMinFirst,tMinFirst,pMinFirst,stepX,stepY,stepT,stepP,i,j,sensorSizeX,sensorSizeY,xo,yo,uxp,uyp)
%This function loops over every bin in 4D space and counts how many
%rays from a single sensor pixel hit each bin. See gathererEfficient.m for
%a more efficient implementation.
%------Inputs--------
%lightField: (2D array) the preallocated histogram array that will be outputted.
%gridX: (scalar) number of bins in x-direction of histogram.
%gridY: (scalar) number of bins in y-direction of histogram.
%gridT: (scalar) number of bins in t-direction of histogram.
%gridP: (scalar) number of bins in p-direction of histogram.
%xMinFirst: (scalar) left most edge of the histogram grid in the x-direction.
%yMinFirst: (scalar) left most edge of the histogram grid in the y-direction.
%tMinFirst: (scalar) left most edge of the histogram grid in the t-direction.
%pMinFirst: (scalar) left most edge of the histogram grid in the p-direction.
%stepX: (scalar) the width of each histogram bin in the x-direction.
%stepY: (scalar) the width of each histogram bin in the y-direction.
%stepT: (scalar) the width of each histogram bin in the theta-direction.
%stepP: (scalar) the width of each histogram bin in the phi-direction.
%i: (scalar) index in the x-direction of the sensor pixel you are currently looking at.
%j: (scalar) index in the y-direction of the sensor pixel you are currently looking at.
%sensorSizeX: (scalar) the number of pixels on the sensor in the
%x-direction.
%sensorSizeY: (scalar) the number of pixels on the sensor in the
%y-direction.
%xo: (column vector) the x-position of the light rays.
%yo: (column vector) the y-position of the light rays.
%uxp: (column vector) angle in degrees in theta-direction of the light rays.
%uyp: (column vector) angle in degrees in phi-direction of the light rays.
%------Outputs--------
%lightfield: (2D array) the histogram where every row is the linear index of a sensor 
%pixel and every column is the linear index of a bin in 4D space. The value in the
%array with index (r,c) represents the number of light rays coming from a single
%sensor pixel r that landed in bin c. 

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