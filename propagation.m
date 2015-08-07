function [yo, xo] = propagation(uyp, uzp, z0, yr, uxp, xr)
%This function takes in direction cosines and positions of light rays 
%as well as a propagation distance and calculates the new positions of the
%light rays after they are propagated.
%------Inputs--------
%uyp: (column vector) y-direction cosines of the light rays.
%uzp: (column vector) z-direction cosines of the light rays.
%z0: (scalar) the propagation distance of the light rays.
%yr: (column vector) the y-position of the light rays.
%uxp: (column vector) x-direction cosines of the light rays.
%xr: (column vector) the x-position of the light rays.
%------Outputs--------
%yo: (column vector) the new y-position of the light rays after propagation.
%xo: (column vector) the new x-position of the light rays after propagation.
    yo = uyp./uzp*z0+yr;
    xo = uxp./uzp*z0+xr;
end