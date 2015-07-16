
%in = load('./dice_137_50-350_color.mat');
%in = load('./dragon_125_200-500_color.mat');

monochrome = 0;   %1 for red, 2 for green, 3 for blue, 0 for color

%A_in = load('../Output/A_sub_100_100_5_5_5e4_CORRECT3.mat');
A_in = load('./A_sub_50_50_5_5_strong.mat')
save_file = 0;
nsx = 500;    %Number of sensor pixels in x
nsy = 500;  %Number of sensor pixels in y
nx = 50;   %Number of light field x bins
ny = 50;   %Number of light field y bins
ntheta = 5; %Number of theta bins
nphi = 5; %Number of phi bins
A_sub1 = A_in.A_sub; 
%A_sub_inv = A_in_inv.A_sub;
clear A_in
clear A_in_inv
A_sub1 = A_sub1./mean(sum(A_sub1,1));   %Normalize A matrix
%A_sub_inv = A_sub_inv./mean(sum(A_sub_inv,1));
%%
lambda = .00001;    %regularization
noise = .001;  %Sensor noise2
col = {'red','green','blue'};   %Color index
h5 = figure(5);
clf
h6 = figure(6);
clf
h7 = figure(7);
clf
if monochrome
    lf_final = zeros(ny*nphi,nx*ntheta,1);
    in = load('../Output/dragon_bunny_50x50.mat');
    if monochrome == 1
        lf = in.lfr;
    elseif monochrome == 2
        lf = in.lfg;
    else 
        lf = in.lgb;   
    end
    clear in
    sensor = A_sub1 * lf(:);
    set(0,'CurrentFigure',h5)
    sensor_reshaped = reshape(sensor,[nsx,nsy]);
    imagesc(sensor_reshaped);
    intens_noisy = sensor + abs(noise*max(sensor)*randn(size(sensor)));
    AtA = A_sub1'*A_sub1;
    AtA_r = (AtA+lambda*speye(size(AtA)));
    tic
    recovered = AtA_r\(A_sub1'*(intens_noisy));
    toc
    recovered_reshaped = reshape(recovered,[nphi,ntheta,nx,ny]);
else  %Color case
    sensor = zeros(nsx,nsy,3);    %Preallocate sensor
    recovered_reshaped = cell(1,3); %Preallocate space for recovered light field
    lf_final = zeros(ny*nphi,nx*ntheta,3);     %Preallocate final light field
    in = load('../Output/dragon_bunny_50x50.mat');   %Load light field
    set(0,'CurrentFigure',h5)
    AtA = A_sub1'*A_sub1;   %Computer A'A
    AtA_r = (AtA+lambda*speye(size(AtA)));   %Add regularization
    lf_cell = cell(3,1);
    lf_cell{1} = in.lfr;
    lf_cell{2} = in.lfg;
    lf_cell{3} = in.lfb;
    tic
    for n = 1:3   %Loop over all three color channels

        
        sensor_mono = A_sub1 * lf_cell{n}(:);   %Compute image using forward model
        
        intens_noisy = sensor_mono + abs(noise*max(sensor_mono)*randn(size(sensor_mono)));   %Add noise
        sensor(:,:,n) = reshape(intens_noisy,[nsx,nsy]);    %Reshape to image
        
        recovered = AtA_r\(A_sub1'*(intens_noisy));  %Invert
        
        recovered_reshaped{n} = reshape(recovered,[nphi,ntheta,nx,ny]); %Reshape inversion to light field (4d)
    end
    toc
    imshow((sensor/prctile(sensor(:),99)));  %Display sensor image
    title('Sensor image')
    axis image
end

%%

%Make lenslet-style image
for n = 1:nx
    for m = 1:ny
        if monochrome
            lf_final((n-1)*ntheta+1:n*ntheta,(m-1)*nphi+1:m*nphi) = (recovered_reshaped(:,:,m,n));
        else
            for ncol = 1:3
                lf_final((n-1)*ntheta+1:n*ntheta,(m-1)*nphi+1:m*nphi,ncol) = uint8(recovered_reshaped{ncol}(:,:,m,n));            
            end
        end
    end
end
set(0,'CurrentFigure',h6)
imagesc(uint8(lf_final))
axis image

title('full recovered light field')


%%
count = 0;

set(0,'CurrentFigure',h7)
if save_file
    str = input('Output file name','s')
    filename = ['./Output/',str,'.gif'] 
end

%Loop over theta and phi, display LF(x,y) at each (theta,phi)
for n = 1:nphi
    for m = 1:ntheta
        count = count+1;
        imagesc(uint8(lf_final(n:5:end,m:5:end,:)))
        if save_file   %Only used if saving output as gif
            if ~monochrome
                [imind,cm] = rgb2ind(uint8(lf_final(n:5:end,m:5:end,:)),256);
                if count == 1;
                  imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
                else
                  imwrite(imind,cm,filename,'gif','WriteMode','append');
                end
            else
                if count == 1;
                  imwrite(uint8(lf_final(n:5:end,m:5:end)),filename,'gif', 'Loopcount',inf);
                else
                  imwrite(uint8(lf_final(n:5:end,m:5:end)),filename,'gif','WriteMode','append');
                end
            end
        end
        axis image
        pause(1/5)
    end
end


