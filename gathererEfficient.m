function lightField = gathererEfficient(lightField,gridX,gridY,gridT,gridP,xMinFirst,...
    yMinFirst,tMinFirst,pMinFirst,stepX,stepY,stepT,stepP,i,j,sensorSizeX,sensorSizeY,xo,yo,uxp,uyp)
%This function calculates the number of rays from a single sensor pixel that
%hit each bin on the diffuser. It is more efficient because it bins all of
%the rays from each sensor pixel simultaneously. See gathererSimple.m for 
%the simple implementation. 

%bins include the right end but do not include the left end
%calculates bin number for the rays in each dimension
xo_r_sub = ceil((xo - xMinFirst) / stepX);
yo_r_sub = ceil((yo - yMinFirst) / stepY);
uxp_r_sub = ceil((uxp - tMinFirst) / stepT);
uyp_r_sub = ceil((uyp - pMinFirst) / stepP);

%checks to see if the bin is in the grid
good = xo_r_sub > 0 & xo_r_sub <= gridX &...
    yo_r_sub > 0 & yo_r_sub <= gridY &...
    uxp_r_sub > 0 & uxp_r_sub <= gridT &...
    uyp_r_sub > 0 & uyp_r_sub <= gridP;

%only keeps bins that are in the grid
xo_r_sub = xo_r_sub(good);
yo_r_sub = yo_r_sub(good);
uxp_r_sub = uxp_r_sub(good);
uyp_r_sub = uyp_r_sub(good);

%gives each ray a linear bin index based on the 4D bin coordinates
ray_ind = sort(sub2ind([gridX,gridY,gridT,gridP],xo_r_sub,yo_r_sub,uxp_r_sub,uyp_r_sub));

%gives bin number and linear index of each bin that contains at least 1 ray
[ray_unique, ray_entry] = unique(ray_ind);

%if rays hit the grid
if nnz(good)
    if length(ray_unique) > 1
        %calculates number of rays in each bin, except the last one
        ray_hist = diff(ray_entry);
        %calculates number of rays in last bin
        ray_hist(end+1) = length(ray_ind)-ray_entry(end)+1;
        %assigns the number of rays that hit each bin to the linear index of the bin 
        lightField(sub2ind([sensorSizeY,sensorSizeX],j+1,i+1),ray_unique) = ray_hist;
    elseif length(ray_unique) == 1
        %if 1 bin want total number of rays that hit the grid
        %or 1 ray that hit the grid
        lightField(sub2ind([sensorSizeY,sensorSizeX],j+1,i+1),ray_unique) = sum(good);
    end
end
%if rays didn't hit the grid, leave the zeros that are used in
%preallocation