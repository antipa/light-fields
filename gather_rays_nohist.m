function gatherer = gather_rays_nohist(xo,yo,npx,npy,dpx,dpy,offsetx,offsety,px)
%left px as input to be consistent with hist4
%Initial offset to go from local to global sensor coordinates. 
nidx = offsetx*px/dpx;   %Multiply by px/dpx to get nidx in sensor coords
midx = offsety*px/dpy;

%Convert ray physical units to pixel subscript (multiple rays not have same
%coordinate at their local pixel center)
xo_r_sub = round(xo/dpx+nidx)-1;
yo_r_sub = round(yo/dpy+midx)-1;

%array that encodes which rays hit the sensor as 1s and those that miss as
%0s
good = xo_r_sub>0&yo_r_sub>0&xo_r_sub<=npx&yo_r_sub<=npy;               

%Throw out rays that don't strike sensor
xo_r_sub = xo_r_sub(good);
yo_r_sub = yo_r_sub(good);

%Sort rays in linear index
ray_ind = sort(sub2ind([npy,npx],yo_r_sub,xo_r_sub));

%Figure out the unique ray indices. ray_entry has the indices within
%ray_ind where each ray occurs. 
[ray_unique, ray_entry] = unique(ray_ind);

%Only do anything if rays hit sensor
if nnz(good)
    if length(ray_unique)>1
        ray_hist = diff(ray_entry);
        ray_hist(end+1) = length(ray_ind)-ray_entry(end)+1;
        [rr, cc] = ind2sub([npy,npx],ray_unique);
        gatherer = sparse(rr,cc,ray_hist,npy,npx);

    elseif length(ray_unique)==1
        gatherer = zeros(npy,npx);
        %The hist function doesn't work for 1 bin, so just
        %count the number of rays that are
        gatherer(ray_unique) = nnz(good);
        gatherer = sparse(gatherer);
    end
else  %If no rays hit sensor, return sparse zeros
    gatherer = sparse(npy,npx);
end
return