% load('/Users/nick.antipa/Documents/MATLAB/Output/1deg_20151121_A_sub_128x128x5x5.mat')
% AtA = A_sub'*A_sub;
%%
crop = 1;
dist = 1;
circ_aper = 0;
if crop
    crop_c = 1280;
    crop_r = 1080;
    crop_half_size = 2^9;
end
nsx = settings.npx;    %Number of sensor pixels in x
nsy = settings.npy;  %Number of sensor pixels in y
nx = settings.N;   %Number of light field x bins
ny = settings.M;   %Number of light field y bins
ntheta = settings.P; %Number of theta bins
save_file= 0;
nphi = settings.Q; %Number of phi bins
monochrome = 1;

im_in = imread('/Users/nick.antipa/Documents/Diffusers/Images_1deg_4f_20151202/flowerReal_z362_0005.tif');

%in = load('../Output/1deg_20151129_A_sub_128_5_z260_gptiesensor_sim.mat');
%im_in = in.s_resh;


figure(9)

if dist
    OA_c = crop_c-crop_half_size+1;
    OA_r = crop_r-crop_half_size+1;
    c = 1/.9955;
    r_off = -2.1;
    c_off = -.58;
    im_x = (1:size(im_in,2))-OA_c;
    im_y = (1:size(im_in,1))-OA_r;
    [im_X,im_Y] = meshgrid(im_x/crop_half_size,im_y/crop_half_size);
    im_rho = sqrt(im_X.^2+im_Y.^2);
    im_th = atan2(im_Y,im_X);
    im_rhop = im_rho*c;
    X_p = im_rhop.*cos(im_th);
    Y_p = im_rhop.*sin(im_th);
    im_in_dist = zeros(size(im_crop));
    for k = 1
        im_in_dist = interp2(im_X*crop_half_size-c_off,im_Y*crop_half_size-r_off,double(im_in),X_p*crop_half_size,Y_p*crop_half_size);
        im_dist = im_in_dist;
    end
    im_dist(isnan(im_dist))=0;
    
end

if crop
    im_crop = double(im_dist(crop_r-crop_half_size+1:crop_r+crop_half_size,crop_c-crop_half_size+1:crop_c+crop_half_size));
else
    im_crop = im_dist;
end
im_crop = im_crop-min(im_crop(:));
%im_crop = im_crop-min(min(im_crop));

imagesc_pctl(im_crop,5,99.9);
axis image
a9 = gca
%im_crop(s_resh==0)=0;
linkaxes([a1,a9],'xy')

% corred = corr2_fft(im_crop,s_resh);
% figure(99),
% contour(corred)
% axis image
% axis([1020 1028 1020 1028])
% grid on
count= 0;
order_vec = [];
for m = 1:nphi
    for n = 1:ntheta
        count = count+1;
        order_vec = cat(2,order_vec,count:ntheta*nphi:nx*ny*ntheta*nphi);
    end
end
%%
im_init = imresize(im_crop,[settings.M,settings.N]);
im_stack = repmat(im_init,[1,1,settings.P,settings.Q]);
im_stack = permute(im_stack,[4,3,2,1]);
vec_init = im_stack(:);
clear im_stack;
clear im_init;
if circ_aper
    for n = 1:settings.Q*settings.P
        if ~ismember(n,good_idx)
            vec_init(n:settings.P*settings.Q:end)=0;
        end
    end
end
    
%im_crop(s_resh==0)=0;

c = 6500;
lamb = 0;
write_gif = 0;
use_twist = 1;
if ~use_twist
    recovered = (AtA+lamb*speye(size(AtA)))\(A_sub'*double(im_crop(:)));
elseif use_twist
    %Phi_func = @(x)(1e8*x'*-(x<0));
    Phi_func = @(x) norm(x,2);
    %Phi_func = @(x) TVnorm_lf(x,ntheta,nphi,ny,nx)
    %Phi_func = @(x) TVnorm_lf(x,order_vec);
    Psi_func = @(x,tau) soft(x,tau);
    %Phi_func = 6
    recovered =  TwIST(double(im_crop(:))/c,A_sub/c,300/c,...
        'StopCriterion',1,...
        'Psi',Psi_func,...
        'Phi',Phi_func,...
        'Debias',1,...
        'MaxiterA',1000,...
        'Monotone',1,...
        'ToleranceA',.000001);
end
  
%% display

h7 = figure(7);
set(0,'CurrentFigure',h7)

%Loop over theta and phi, display LF(x,y) at each (theta,phi)
count = 0
for m = 1:nphi
    for n = 1:ntheta
        count = count+1; 
        imagesc(reshape(recovered(count:ntheta*nphi:end),nx,ny)')
        %caxis([0 35])
        colormap gray
        axis image
        caxis([0 5])
        pause(1/9)
    end
    
end

%% Try refocusing. Because why not.
h8 = figure(8)
imup = 1
im_x = 1:nx;
im_y = 1:ny;
[IMX, IMY] = meshgrid(im_x,im_y);
for s = 1
    %s =  -1
    refocused = zeros(settings.N*imup,settings.M*imup);
    count = 0;
    for m = 1:nphi
        for n = 1:ntheta
            count = count+1;
            sub_im = transpose(imresize(reshape(recovered(count:ntheta*nphi:end),ny,nx),imup));
            
            %im_shifted = circshift(sub_im,[-round(s*(m-nphi/2-1/2)),-round(s*(n-ntheta/2-1/2))]);
            im_shifted = interp2(IMX,IMY,sub_im,IMX+s*(n-ntheta/2-1/2),IMY+s*(m-nphi/2-1/2));
            %round(s*(m-nphi/2-1/2))
            %round(s*(n-ntheta/2-1/2));
            refocused = im_shifted+refocused;
            %set(0,'CurrentFigure',h8)
            %imagesc(refocused/count)
            %colormap gray
            %caxis([0 3])
            %drawnow
            %pause(1/10)
        end
    end
    imagesc((refocused/prctile(refocused(:),99.9)).^(1/.7))
    colormap gray
    axis image
    caxis([0 1])
        pause(1/50)
    drawnow

end