function gatherer = hist4(xo,yo,npx,npy,dpx,dpy,offsetx,offsety,px)
centers = {};
centers{1} = 0:dpx:x_range - dpx;
centers{2} = 0:dpy:y_range - dpy;
globalx = xo + (nidx(1)-1)*px;
globaly = yo + (midx(1)-1)*px;
good = globalx >= - dpx/2 & globalx <= x_range & globaly >= -dpy/2 & globaly <= y_range;
xo1 = globalx(good);
yo1 = globaly(good);
matrix = [yo1 xo1];
gatherer = hist3(matrix,centers);
end