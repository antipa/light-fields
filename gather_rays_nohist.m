function gatherer = gather_rays_nohist(xo,yo,npx,npy,dpx,dpy,offsetx,offsety,px)
nidx = offsetx;
midx = offsety;
%bin outputs at sensor pixels in local physical units
xo_r = (xo+dpx/2-mod(xo+dpx/2,dpx));
yo_r = (yo+dpy/2-mod(yo+dpy/2,dpy));

%Shift physically located rays to global coordinates
xo_r_global = xo_r + (nidx(1)-1)*px;
yo_r_global = yo_r + (midx(1)-1)*px;     

%xo_r_sub = 1+round((xo_r-min(xo_r))/dpx);
%yo_r_sub = 1+round((yo_r-min(yo_r))/dpy);
xo_r_sub = 1+round((xo_r_global)/dpx);
yo_r_sub = 1+round((yo_r_global)/dpy);
good = xo_r_sub>0&yo_r_sub>0&xo_r_sub<=npx&yo_r_sub<=npy;               
xo_r_sub = xo_r_sub(good);
yo_r_sub = yo_r_sub(good);
%gatherer = zeros(npy,npx);

ray_ind = sort(sub2ind([npy,npx],yo_r_sub,xo_r_sub));
[ray_unique, ray_entry] = unique(ray_ind);

if nnz(good)
    if length(ray_unique)>1
        ray_hist = diff(ray_entry);  
        ray_hist(end+1) = length(ray_ind)-ray_entry(end);
        [rr, cc] = ind2sub([npy,npx],ray_unique);
        gatherer = sparse(rr,cc,ray_hist,npy,npx);


    elseif length(ray_unique)==1
        gatherer = zeros(npy,npx);
        %The hist function doesn't work for 1 bin, so just
        %count the number of rays that are
        gatherer(ray_unique) = nnz(good);
        gatherer = sparse(gatherer);
    end
else
    gatherer = sparse(npy,npx);
end
return