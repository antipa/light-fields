%load H matrix generated from diffuserPointSpread3D.m
in1 = load('H_5deg.mat');
hMatrix = in1.hMatrix;
in2 = load('H_5deg_2.mat');
hMatrix2 = in2.hMatrix;

%%
%object size in x, y, and z from diffuserPointSpread3D.m
sizeX = 10;
sizeY = 10;
sizeZ = 10;

%sensor size from diffuserPointSpread3D.m
sensorSizeX = 50;
sensorSizeY = 50;

%noise added to sensor
noise = 0.000;

%regularization
lambda = 0.0005;

%read in image
%in = imread('Tanzania-0741.jpg','jpg');
% r = double(in(:,:,1));
% g = double(in(:,:,2));
% b = double(in(:,:,3));
% avg = (r + g + b) ./ 3;
% bw = uint8(avg);
%bw = rgb2gray(in);

%shrink or stretch image to size of object when H matrix was created
bw = zeros(10,10,10)
[r, c] = size(bw);
rowVector = linspace(1,r,sizeY);
rowVector = round(rowVector);
columnVector = linspace(1,c,sizeX);
columnVector = round(columnVector);
bw = bw(rowVector,columnVector);
bw(:,:) = 0;
bw(3,3) = 255;
bw(3,7) = 255;
bw(6,3) = 255;
bw(7,4:6) = 255;
bw(6,7) = 255;

clf
figure(1);
imshow(bw);

%crop image to size of object when H matrix was created
% bw = bw(1000:1009,700:709);
% clf 
% figure(1);
% imshow(bw);

%all ones(white) in 3rd dimension because picture is 2D
bw = im2double(bw);
a = zeros(sizeY,sizeX,sizeZ);
a(:,:,1) = bw;
a = reshape(a,[sizeY*sizeX*sizeZ,1]);

sensorImage = hMatrix2 * a;
sensorImage_reshaped = reshape(sensorImage, [sensorSizeY, sensorSizeX]);
figure(2);
colormap default;
imagesc(sensorImage_reshaped);

%add noise to sensor image and invert and reshape
sensorImage_noisy = sensorImage + abs(noise*max(sensorImage)*randn(size(sensorImage)));

%%
cvx_begin
variable recovered(size(hMatrix,2),1)
minimize norm(hMatrix*recovered - sensorImage_noisy) + 10*norm(recovered,1)
subject to
recovered>=0
cvx_end
%%

%recovered = (hMatrix' * hMatrix + )\(hMatrix' * sensorImage_noisy);
recovered_reshaped = reshape(recovered,[sizeY,sizeX,sizeZ]);
figure(3);
colormap gray;
imagesc(recovered_reshaped(:,:,1));