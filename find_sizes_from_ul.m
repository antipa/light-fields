function [M_f, N_f] = find_sizes_from_ul(C0_updated,im_w,im_h,dpx_p,dpy_p,uxp,uyp)
d = [im_w; im_h];

cvx_begin quiet
variable N(1)
variable M(1)

maximize dpx_p*N+dpy_p*M
subject to
C = C0_updated;
C+N*dpx_p*uxp<=d;
C+N*dpx_p*uxp>=1;
C+M*dpy_p*uyp<=d;
C+M*dpy_p*uyp>=1;
C+N*dpx_p*uxp+M*dpy_p*uyp<=d;
C+N*dpx_p*uxp+M*dpy_p*uyp>=1;
cvx_end

M_f = floor(M);
N_f = floor(N);

end