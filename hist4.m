function gatherer = hist4(xo,yo,npx,npy,dpx,dpy,offsetx,offsety,px)
nidx = offsetx;
midx = offsety;
centers = {};
centers{1} = 0:dpx:dpx*npx - dpx;
centers{2} = 0:dpy:dpy*npy - dpy;
globalx = xo + (nidx(1)-1)*px;
globaly = yo + (midx(1)-1)*px;
good = globalx >= - dpx/2 & globalx <= dpx*npx & globaly >= -dpy/2 & globaly <= dpx*npx;
xo1 = globalx(good);
yo1 = globaly(good);
matrix = [yo1 xo1];
gatherer = hist3(matrix,centers);
end