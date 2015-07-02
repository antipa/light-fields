function gatherer = gather_rays_nohist(xo,yo,npx,npy,dpx,dpy,offsetx,offsety,px)
%left px as input to be consistent with hist4
nidx = offsetx;
midx = offsety;   
xo_r_sub = round(xo/dpx)+nidx(1);
yo_r_sub = round(yo/dpy)+midx(1);
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