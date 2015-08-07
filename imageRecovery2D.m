%object size in x, y, and z from diffuserPointSpread3D.m
sizeX = 10;
sizeY = 10;
sizeZ = 10;

%sensor size from diffuserPointSpread3D.m
sensorSizeX = 50;
sensorSizeY = 50;

%if want to save recovery as .gif
save_file = true;

%noise added to sensor
noise = 0.001;

%regularization
lambda = 500000;

 %read in image
% in = imread('Tanzania-0741.jpg','jpg');
% r = double(in(:,:,1));
% g = double(in(:,:,2));
% b = double(in(:,:,3));
% avg = (r + g + b) ./ 3;
% bw = uint8(avg);
% bw = rgb2gray(in);

%shrink or stretch image to size of object when H matrix was created
% [r, c] = size(bw);
% rowVector = linspace(1,r,sizeY);
% rowVector = round(rowVector);
% columnVector = linspace(1,c,sizeX);
% columnVector = round(columnVector);
% bw = bw(rowVector,columnVector);


%crop image to size of object when H matrix was created
% bw = bw(1000:1009,700:709);
% clf
% figure(1);
% imshow(bw);

%creating the z-planes of the object
%bw is the plane with the smiley face
bw(sizeY,sizeX) = 0;
bw2 = bw;
bw(3,3) = 255;
bw(3,7) = 255;
bw(6,3) = 255;
bw(7,4:6) = 255;
bw(6,7) = 255;
bw = im2double(bw);

%bw2 is the plane with the horizontal line
bw2(6,3:8) = 255;
bw2 = im2double(bw2);

% clf
% figure(1);
% imshow(bw);

%a is the 3D object
a = zeros(sizeY,sizeX,sizeZ);
a(:,:,10) = bw;
a(:,:,5) = bw2;
a = reshape(a,[sizeY*sizeX*sizeZ,1]);

%resulting sensor image
sensorImage = hMatrix2 * a;
sensorImage_reshaped = reshape(sensorImage, [sensorSizeY, sensorSizeX]);
figure(2);
colormap default;
imagesc(sensorImage_reshaped);

%add noise to sensor image
sensorImage_noisy = sensorImage + noise*max(sensorImage)*randn(size(sensorImage));

%invert and recover the object as a column vector
recovered = (hMatrix' * hMatrix + lambda * eye(size(hMatrix' * hMatrix)))\(hMatrix' * sensorImage_noisy);
recovered_reshaped = reshape(recovered,[sizeY,sizeX,sizeZ]);
h3 = figure(3);
colormap gray;

set(0,'CurrentFigure',h3)

if save_file
    str = input('Output file name','s');
    filename = [str '.gif'];
end

%cycle through the z-planes of the object to create gif 
for v = 1:sizeZ
    maxRecovered = max(recovered);
    imagesc(recovered_reshaped(:,:,v));
    caxis([0 maxRecovered]);
    im = recovered_reshaped(:,:,v);
     if save_file
        if v == 1;
             %adjust z-planes to have same colormap scale
            imwrite(im*255/maxRecovered,filename,'gif', 'Loopcount',inf);
        else
             imwrite(im*255/maxRecovered,filename,'gif','WriteMode','append');
        end
     end
    pause(1/2);
end