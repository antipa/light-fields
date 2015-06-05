function [yo, xo] = propagation(uyp, uzp, z0, yr, uxp, xr)
    yo = uyp./uzp*z0+yr;
    xo = uxp./uzp*z0+xr;
end