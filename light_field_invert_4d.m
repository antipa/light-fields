
%in = load('./dice_137_50-350_color.mat');
%in = load('./dragon_125_200-500_color.mat');

monochrome = 0;   %1 for red, 2 for green, 3 for blue, 0 for color

%A_in = load('../Output/A_sub_100_100_5_5_5e4_CORRECT3.mat');
%A_in = load('../Output/A_under_50x50x5x5_125x125_far_z.mat');
solver = 'lsmr';
Phi_func = @(x)sum(diff(x));
save_file = 0;
nsx = A_in.settings.npx;    %Number of sensor pixels in x
nsy = A_in.settings.npy;  %Number of sensor pixels in y
nx = A_in.settings.N;   %Number of light field x bins
ny = A_in.settings.M;   %Number of light field y bins
ntheta = A_in.settings.P; %Number of theta bins
nphi = A_in.settings.Q; %Number of phi bins
A_sub1 = A_in.A_sub; 
%A_sub_inv = A_in_inv.A_sub;
clear A_in
clear A_in_inv
A_sub1 = A_sub1./mean(sum(A_sub1,1))*1e0;   %Normalize A matrix


%A_sub_inv = A_sub_inv./mean(sum(A_sub_inv,1));
%%
lambda = .0001;    %regularization
noise = .001;  %Sensor noise2
col = {'red','green','blue'};   %Color index
h5 = figure(5);
clf
h6 = figure(6);
clf
h7 = figure(7);
clf
monochrome = 1;

if monochrome
    lf_final = zeros(ny*nphi,nx*ntheta,1);
    in = load('../Output/test.mat')
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
    in = load('../Output/dragon_tooth_50x50x5x5.mat');
    set(0,'CurrentFigure',h5)
    
    lf_cell = cell(3,1);
    lf_cell{1} = in.lfr;
    lf_cell{2} = in.lfg;
    lf_cell{3} = in.lfb;
    tic
    for n = 1:3   %Loop over all three color channels

        [C,L] = wavedec(lf_cell{n}(:),3,'bior4.4');
        [C_t,L_t] = wavedec(lf_cell{n}(:),5,'rbio4.4');
        %A_func = @(x)A_sub1*waverec(padarray(x',length(C)-size(A_sub1,2),0,'post'),L,'bior4.4');
        %AT_func = @(y)transpose(waverec(padarray(A_sub1'*y,length(C)-size(A_sub1,2),0,'post'),L,'rbio4.4'));
        sensor_mono = A_sub1 * lf_cell{n}(:);   %Compute image using forward model
        
        intens_noisy = sensor_mono + abs(noise*max(sensor_mono)*randn(size(sensor_mono)));   %Add noise
        sensor(:,:,n) = reshape(intens_noisy,[nsx,nsy]);    %Reshape to image
        switch lower(solver)
            case('twist')
                fprintf('using twist\n')
                recovered =  TwIST(intens_noisy,A_sub1,5e-2,'Phi',Phi_func,'StopCriterion',1,'ToleranceA',.00001);
            
            case('back_divide')
                AtA = A_sub1'*A_sub1;   %Computer A'A
                AtA_r = (AtA+lambda*speye(size(AtA)));   %Add regularization
                recovered = AtA_r\(A_sub1'*(intens_noisy));  %Invert
                rec_structure = 'vector';
            case('lsmr')
                recovered = lsmr(A,b,1e-5,[],[],[],[],[],1);
                rec_structure = 'vector';
        end
        switch lower(rec_structure)
            case('vector')
                recovered_reshaped{n} = reshape(recovered,[nphi,ntheta,nx,ny]); %Reshape inversion to light field (4d)
            case('stack')
                recovered_reshaped{n} = 
            case('wavelet')
        end
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
    filename = ['../Output/',str,'.gif'] 
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


