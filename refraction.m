function [uxp, uyp, uzp] = refraction(Fxr, Fyr, th, ph, index)
%Normal vectors. ith row is [x,y,z] normal at (xr(i),yr(i),zr(i)
normals_norm = sqrt(Fxr.^2+Fyr.^2+1);   %Length of each vector
normals = [-Fxr./normals_norm,-Fyr./normals_norm,ones(size(Fxr))./normals_norm];

%Convert theta and phi from degrees into vector representation
ux = tand(th);
uy = tand(ph);
uz = ones(size(ux));
norms = sqrt(ux.^2+uy.^2+1);

%Normalize (probably not necessary?) to get direction
%cosines
uxn = ux./norms;
uyn = uy./norms;
uzn = uz./norms;


%Calculate magnitude of incident angle I

I = acos(sum(normals.*[uxn, uyn, uzn],2));
%Use snell's law to calculate Ip
index_p = 1;
Ip = asin(index/index_p*sin(I));
%define gamma = n'cosI'-ncosI
Gamma = index_p*cos(Ip)-index*cos(I);

%Calculate new direction cosines
uxp = 1/index_p * (index*uxn+Gamma.*normals(:,1));
uyp = 1/index_p * (index*uyn+Gamma.*normals(:,2));
uzp = 1/index_p * (index*uzn+Gamma.*normals(:,3));

end
