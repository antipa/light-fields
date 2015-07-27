%load H matrix generated from diffuserPointSpread3D.m
load('H2.mat');

%object size in x, y, and z from diffuserPointSpread3D.m
sizeX = 10;
sizeY = 10;
sizeZ = 10;

%sensor size from diffuserPointSpread3D.m
sensorSizeX = 50;
sensorSizeY = 50;

%read in image
in = imread('Tanzania-0741.jpg','jpg');
% r = double(in(:,:,1));
% g = double(in(:,:,2));
% b = double(in(:,:,3));
% avg = (r + g + b) ./ 3;
% bw = uint8(avg);
bw = rgb2gray(in);

%shrink or stretch image to size of object when H matrix was created
[r, c] = size(bw);
rowVector = linspace(1,r,sizeY);
rowVector = round(rowVector);
columnVector = linspace(1,c,sizeX);
columnVector = round(columnVector);
bw = bw(rowVector,columnVector);
clf
figure(1);
imshow(bw);

%crop image to size of object when H matrix was created
% bw = bw(1000:1009,700:709);
% clf 
% figure(1);
% imshow(bw);

%all zeros in 3rd dimension because picture is 2D
bw = im2double(bw);
a = zeros(sizeY,sizeX,sizeZ);
a(:,:,1) = bw;
a = reshape(a,[sizeY*sizeX*sizeZ,1]);

sensorImage = hMatrix * a;
sensorImage = reshape(sensorImage, [sensorSizeY, sensorSizeX]);
figure(2);
colormap default;
imagesc(sensorImage);
