function [uxp, uyp, uzp] = refraction(Fxr, Fyr, th, ph, index, index_p, varargin)
%This function is for refracting light rays through a surface and getting
%the output direction cosines.
%------Inputs--------
%Fxr: (column vector) gradient of the surface in x-direction.
%Fyr: (column vector) gradient of the surface in y-direction.
%th: (column vector) angle in degrees in theta direction of the light rays (see optional inputs).
%ph: (column vector) angle in degrees in phi direction of the light rays (see optional inputs).
%index: (scalar) index of refraction of the medium the rays start in.
%index_p: (scalar) index of refraction of the medium the rays end in.
%------Optional Inputs------
%First optional input is either 'angles' or 'cosines'.
%If using 'cosines', th and ph inputs are cosine direction column vectors in the
%x- and y-directions respectively. Must also have a second optional input, 
%which is a z-direction cosine vector (column vector).
%If no optional inputs, assuming th and ph inputs are angles in degrees.
%-----Outputs---------
%uxp: (column vector) x-direction cosines after refraction.
%uyp: (column vector) y-direction cosines after refraction.
%uzp: (column vector) z-direction cosines after refraction.

if strcmp(varargin{1},'cosines')
  uxn = th;
  uyn = ph;
  if isMatrix(varargin{2})
  uzn = varargin{2};
  elseif isMatrix(varargin{1})
      uzn = varargin{1};
  end
  
elseif strcmp(varargin{1},'angles') || nargin == 6
%Normal vectors. ith row is [x,y,z] normal at (xr(i),yr(i),zr(i)
normals_norm = sqrt(Fxr.^2+Fyr.^2+1);   %Length of each vector
normals = [-Fxr./normals_norm,-Fyr./normals_norm,ones(size(Fxr))./normals_norm];

%Convert theta and phi from degrees into vector representation
ux = tand(th);
uy = tand(ph);
uz = ones(size(ux));
norms = sqrt(ux.^2+uy.^2+1);

%Normalize to get direction
%cosines
uxn = ux./norms;
uyn = uy./norms;
uzn = uz./norms;

end

%Calculate magnitude of incident angle I

I = acos(sum(normals.*[uxn, uyn, uzn],2));
%Use snell's law to calculate Ip
Ip = asin(index/index_p*sin(I));
%define gamma = n'cosI'-ncosI
Gamma = index_p*cos(Ip)-index*cos(I);

%Calculate new direction cosines
uxp = 1/index_p * (index*uxn+Gamma.*normals(:,1));
uyp = 1/index_p * (index*uyn+Gamma.*normals(:,2));
uzp = 1/index_p * (index*uzn+Gamma.*normals(:,3));

end
