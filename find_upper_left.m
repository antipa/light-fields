function c = find_upper_left(xgm,im_w,im_h,dpx_p,dpy_p,uxp,uyp,ulx,uly)
    d = [im_w; im_h];
    C0 = xgm';
    
    cvx_begin quiet
    variable N(1)
    variable M(1) 
    variable P(1)
    variable Q(1)
    maximize dpx_p*N+dpy_p*M
    subject to
    C = C0+P*dpx_p*uxp+Q*dpy_p*uyp;
    C>=[ulx;uly];
    C+N*dpx_p*uxp<=d;
    C+N*dpx_p*uxp>=[ulx;uly];
    C+M*dpy_p*uyp<=d;
    C+M*dpy_p*uyp>=[ulx;uly];
    C+N*dpx_p*uxp+M*dpy_p*uyp<=d;
    C+N*dpx_p*uxp+M*dpy_p*uyp>=[ulx;uly];
    cvx_end
%     c = zeros(2,4);
%     c(:,1) = C0 + floor(P)*dpx_p*uxp + floor(Q)*dpy_p*uyp;
%     c(:,2)= C0 + ceil(P)*dpx_p*uxp + floor(Q)*dpy_p*uyp;
%     c(:,3)= C0 + floor(P)*dpx_p*uxp + ceil(Q)*dpy_p*uyp;
    c = C0 + ceil(P)*dpx_p*uxp + ceil(Q)*dpy_p*uyp;
    
 return   

end