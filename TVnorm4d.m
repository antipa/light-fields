function y = TVnorm4d(x,angweight,apfunc)
%4D TV norm with adjustible angular 
% Assumes dimensions 1 and 2 are angular.
% angweight: amplifies angular contribution by amount angweight
% apfunc: 4D array with NaN where the aperture doesn't exist. Any
% difference between a value in the aperture and NaN will be zeroed out.

if isempty(apfunc)
    apfunc = ones(size(x));
end
ss = size(x);
d1 = zeros(ss);
d2 = zeros(ss);
d3 = zeros(ss);
d4 = zeros(ss);
d1(1:end-1,:,:,:) = diff(x.*apfunc,1,1)*angweight;
d1(isnan(d1)) = 0;
d2(:,1:end-1,:,:) = diff(x.*apfunc,1,2)*angweight;
d2(isnan(d2)) = 0;
d3(:,:,1:end-1,:) = diff(x.*apfunc,1,3);
d3(isnan(d3)) = 0;
d4(:,:,:,1:end-1) = diff(x.*apfunc,1,4);
d4(isnan(d4)) = 0;
y = sum(sum(sum(sum(sqrt(d1.^2+d2.^2+d3.^2+d4.^2)))));

return
