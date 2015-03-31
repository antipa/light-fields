NPhi = 16;
NX = 300;
npix = 5000;
%in = load('./dice_137_50-350_color.mat');
%in = load('./dragon_125_200-500_color.mat');
in = load('./dragon_bunny_230_510-810_color.mat');
lfr = imresize(in.lfr,[NPhi, NX],'bicubic');
%lfr = 255*ones(size(lfr));
monochrome = 0;   %1 for red, 2 for green, 3 for blue, 0 for color
%A_in = load('./A_sub_no_diffuser_16x300_10um_z0_228um.mat');
A_in = load('./A_sub_16x300_10um_z0_60um.mat');
%A_in = load('./A_sub_1_degree_16x300_10um_z0_228um.mat');
%A_in = load('A_sub_no_diffuser_16x300_10um_z0_228um.mat');
%A_in = load('A_sub_no_diffuser_16x300_10um_z0_228um_2e6.mat');
%A_in = load('A_sub_16x300_10um_z0_60um_1e5_image_wide_gradient');
A_sub1 = A_in.A_sub;
A_sub1 = A_sub1./mean(sum(A_sub1,1));   %Normalize A matrix
lambda = .005;    %regularization
noise = .01;  %Sensor noise


if monochrome
    lf_m = lfr(:,:,monochrome);
    sensor = A_sub1 * lf_m(:);
    figure(5),clf,plot(sensor)
else
    sensor = zeros(npix,3);
    recovered_reshaped = zeros(NPhi,NX,3);
    h5 = figure(5);
    clf,hold on
    h6 = figure(6);
    clf,hold on
    for m = 1:3
        lf_m = lfr(:,:,m);
        sensor(:,m) = A_sub1 * lf_m(:);

        intens_noisy = sensor(:,m)+noise*max(sensor(:,m))*randn(size(sensor(:,m)));
        if m ==1
            set(0,'CurrentFigure',h5)            
            stairs(sensor(:,1),'r-')
            title('sensor output')
            xlabel('pixel')
            ylabel('intensity')
        elseif m == 2
            set(0,'CurrentFigure',h5) 
            stairs(sensor(:,2),'g-')
        else
            set(0,'CurrentFigure',h5) 
            stairs(sensor(:,3),'b-')
        end
        
        recovered = (A_sub1'*A_sub1+lambda*eye(size(A_sub1,2),size(A_sub1,2)))\(A_sub1'*(intens_noisy));
        recovered_reshaped(:,:,m) = reshape(recovered,NPhi,NX);        
    end
end
%%
set(0,'CurrentFigure',h6)
subplot(2,1,1)
imagesc(uint8(lfr+noise*randn(NPhi,NX,3)*255));
title('Original light field with noise added')
subplot(2,1,2)
imagesc(uint8(recovered_reshaped))
title(['Recovered light field with ',num2str(noise*100),'% sensor noise'])

cmax = .2;
figure(7),imagesc(uint8(255*imadjust((abs(lfr/255-recovered_reshaped/255)),[0;cmax],[0 1])))
title(['abs(image-recovered) on 0 to ', num2str(uint8(255*cmax)),' bit scale'])

figure(h6)