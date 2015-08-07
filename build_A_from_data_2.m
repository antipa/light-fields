%%
%Find center of each screen pixel on camera from diffuser-free image. These
%centers will serve as the centers for cropping each diffuser psf.
precompute = 0;
compute_roi = 0;
vis_grid = 0;
npeaks = 4;
color = 2;   %Use 1 for red, 2 for green, 3 for blue.
demosaic_type = 'rggb';
cam_type = 'mono';
im_downsize = .1;   %Use 1 for no downsizing
roi_scale = 1;
bg_size = 0;

%Distortion numbers
% d = -0.03;
% %Define optical axis center manually
% OA_c = 1200;
% OA_r = 1200; 

%Global distortion correction
d = 0;
OA_r = 1000
OA_c = 1200;

%off axis distortion correction - This will only get used for correcting
%perspective and disortion that is pupil-coordinate dependent. It will
%apply a cubic distortion correction by the amount offaxis_distortion
%centered at OA_c+OA_dir*OA_c_offset and OA_r+OA_dir*OA_r_offset--OA_dir
%can be used to flip direction. WHereas global distortion correction will
%happen to all images from all angles, and must be applied to the
%multiplexed image prior to inversion, the off axis correction will be
%applied only to the cropping image and will NOT warp the image pixel
%space. Rather, it is intended to warp the grid-based cropping schemed TO
%the warped image (as opposed to warping the image to a grid). 

offaxis_distortion =  -.000000011;
offaxis_distortion_quadratic = .000021;
% offaxis_distortion =  0;
%  offaxis_distortion_quadratic = 0;

OA_c_offset = 900
OA_r_offset = 0;
OA_dir = 1;   %Use -1 for negative OA_offset

%Image used to determine grid
grid_image = '/Users/nick.antipa/Documents/Diffusers/1deg_stereo/calib_both_true.tif';
grid_image = [];



%fname = '/Users/nick.antipa/Google Drive/Diffuser data/5_by_5_4x4_shift_0p5_deg_closed.tif'
fname = '/Users/nick.antipa/Documents/Diffusers/half_deg_stereo_20150804/left_white_5px_6x6.tif';
info = imfinfo(fname);
im_start = 5; 
nshifts = 6;
im_list = im_start:2:(2*nshifts^2+im_start-1);
im_in = zeros(info(1).Height,info(1).Width,numel(im_list));
count = 0;
for m = 1:numel(im_list);
    count = count+1;
    k = im_list(m);
    im_in(:,:,count) = imread(fname, k, 'Info', info);
    imagesc(im_in(:,:,count))
    pause(1/100)
end


%Find upper left row and column, as well as width and height of region of
%interest
if compute_roi
    [UL_r, UL_c, s, mat] = findMaxSq(roi_o);
    ulr = UL_r(1);
    ulc = UL_c(1);
    roi_w = s;
    roi_h = s;
else
    ulr = 1;
    ulc = 1;
    roi_w = info(1).Width*roi_scale;
    roi_h = info(1).Height*roi_scale; 
end
im_w = info(1).Width;
im_h = info(1).Height;
%% Correct geometric distortion
if d~=0
    im_x = (1:im_w)-OA_c;
    im_y = (1:im_h)-OA_r;
    [im_X,im_Y] = meshgrid(im_x,im_y);
    im_rho = sqrt(im_X.^2+im_Y.^2);
    im_rho = im_rho;
    im_th = atan2(im_Y,im_X);
    im_rhop = im_rho + max(im_rho(:))*d*(im_rho/max(im_rho(:))).^3;
    X_p = im_rhop.*cos(im_th);
    Y_p = im_rhop.*sin(im_th);
    %To do: X and Y are currently normalized. FIx that!
    im_in_dist = zeros(size(im_in));
    for k = 1:size(im_in,3)
        im_in_dist(:,:,k) = interp2(im_X,im_Y,im_in(:,:,k),X_p,Y_p);
        imagesc(im_in_dist(:,:,k));
        axis image
        pause(1/100)
    end
else
    im_in_dist = im_in;
end
%%
clear im_in
clear im_rhop
clear im_rho


%% Figure out grid

x = 1:im_w;
y = 1:im_h;
pad = 5;
if isempty(grid_image)
    spect = fftshift(fft2(im_in_dist(:,:,1)-mean2(im_in_dist(:,:,1)),pad*im_h,pad*im_w));
else
    grid_im_in = imread(grid_image);
    spect = fftshift(fft2(grid_im_in-mean2(grid_im_in),pad*im_h,pad*im_w));
end
fx = linspace(-1/2/mean(diff(x)),1/2/mean(diff(x)),pad*numel(x));
fy = linspace(-1/2/mean(diff(y)),1/2/mean(diff(y)),pad*numel(y));

[FX, FY] = meshgrid(fx,fy);

%Find top 6 peaks within x roi
%%
peak_region = .002; %In frequency units
roi_width = .005;
peak_reg_idx = round(peak_region/mean(diff(fx)));
x_roi = abs(FX)<roi_width & abs(FY)>peak_region;
x_band = abs(spect).*x_roi;
f_roi = 1000;
pr = zeros(npeaks,1);
pc = pr;
for n = 1:npeaks
    mv = max(x_band(:));
    [pr(n), pc(n)] = find(x_band==mv,1,'first');
    x_band(pr(n)-peak_reg_idx:pr(n)+peak_reg_idx,pc(n)-peak_reg_idx:pc(n)+peak_reg_idx)=0;
    imagesc(x_band)
    caxis([0 1e9])
    axis([pad*middle(x)-f_roi pad*middle(x)+f_roi pad*middle(y)-f_roi pad*middle(y)+f_roi])
    pause(1/20)
end
xx = [0;mean(diff(fx))*(pc-max(x)*pad/2-1)];
xy = [0;mean(diff(fy))*(pr-max(y)*pad/2-1)];

%%
y_roi = abs(FY)<roi_width & abs(FX)>peak_region;
y_band = abs(spect).*y_roi;

pr = zeros(npeaks,1);
pc = pr;
for n = 1:npeaks
    mv = max(y_band(:));
    [pr(n), pc(n)] = find(y_band==mv,1,'first');
    y_band(pr(n)-peak_reg_idx:pr(n)+peak_reg_idx,pc(n)-peak_reg_idx:pc(n)+peak_reg_idx)=0;
    imagesc(y_band)
    caxis([0 1e9])
    axis([pad*middle(x)-f_roi pad*middle(x)+f_roi pad*middle(y)-f_roi pad*middle(y)+f_roi])
    pause(1/20)
end
yx = [0;mean(diff(fx))*(pc-max(x)*pad/2-1)];
yy = [0;mean(diff(fy))*(pr-max(y)*pad/2-1)];

%% figure out rotation
clf
imagesc(abs(spect),'XData',fx,'YData',fx)
set(gca,'YDir','normal')
hold on
caxis([0 4e9])
scatter(yx,yy,'rx')
scatter(xx,xy,'k+')

my = xx'*xy/(xx'*xx);
ly = line(1/my*[min(fx),max(fx)],[min(fx), max(fx)]);
ly.Color = 'k'

mx = yx'*yy/(yx'*yx);
lx = line([min(fx), max(fx)],mx*[min(fx),max(fx)]);
lx.Color = 'r'
hold off

clear spect
clear FX
clear FY
clear x_band
clear y_band
clear x_roi
clear y_roi
clear im_th
%% figure out period
ufx = [1;mx];
ufx = ufx/norm(ufx);

%Check for infinite slope
if ~isnan(my)  
    ufy = [1/my;1];
    ufy = ufy/norm(ufy);
else
    ufy = [0;
           1];
end

xproj = sort(ufy'*[xx';xy']);
yproj = sort(ufx'*[yx';yy']);

dfx = mean(diff(xproj));
dfy = mean(diff(yproj));

dpx = 1/dfx;
dpy = 1/dfy;



%% display grid
h3 = figure(3);
set(0,'CurrentFigure',h3)
clf

% xg = a_good(1:2:end);
% 
% yg = a_good(2:2:end);
im_calib = im_in_dist(:,:,1);
clf, imagesc(imresize(im_calib,im_downsize))
axis image
hold on

%Calculate grid line directions
if ~isnan(my)
    ux = [1,-1/my];
    ux = ux/norm(ux);   
else
    ux = [1 0];
end
uxp = ux';

uy = [-mx,1];
uy = uy/norm(uy);
uyp = uy';
%Find point in center (to minimize distortion issues with fitting grid)

% nlinesx = 14*2;
% nlinesy = 12*2;
xgm = [1100 1073]; %for 1/2 deg white 20150804 left
%xgm = [1011 1060];for 1/2 deg white 20150804 right
%xgm = [1084 1057];   %For 1 deg white 6x6 right 20150804
%xgm = [1868 1050];    %For 1 deg white 6x6 left 20150804
%xgm = [1269 1087];   1/2 deg %xy coord of grid center
%xgm = [1300 1013];  %xy coord of center for 1 deg diffuser
ul_saved = zeros(2,nshifts^2);
nx_saved = zeros(1,nshifts^2);
ny_saved = zeros(1,nshifts^2);
for ii = 1:nshifts
    for jj = 1:nshifts
        set(0,'CurrentFigure',h3)
        %clf
        %imagesc(im_in_dist(:,:,4*(ii-1)+jj))
        %hold on
        
        ind = nshifts*(ii-1)+jj;
        imagesc(imresize(im_in_dist(:,:,ind),im_downsize));
        axis image
        x3 = [];
        x5 = [];
        xgm_adj = xgm-(jj-1)*dpx*ux/nshifts-(ii-1)*dpy*uy/nshifts;
        
        %Use CVX to find coordinates of upper left corner by varying width,
        %height and deviation from known point, xgm_adj in integer steps of
        %uxp and uyp
        c = find_upper_left(xgm_adj,roi_w,roi_h,dpx,dpy,uxp,uyp,ulc,ulr);
        
        %Find optimal height and width of rectangle from fixed upper left
        %corner
        [nlinesy,nlinesx] = find_sizes_from_ul(c,ulc+roi_w,ulc+roi_h,dpx,dpy,uxp,uyp)
        
        x0 = c';
        %Generate coordinates of each point
        %xpp is a list of all xy points, each column is a point.
        xpp = grid_intersections(x0,ux,uy,dpx,dpy,nlinesx,nlinesy);
        
        if offaxis_distortion~=0
%             bg = double(opening(unif(im_calib,5),75));
%             b = watershed(2e16-gaussf(im_calib,21),2);
%             lab = bwlabel(~double(b));
%             msr = regionprops(lab,im_calib-bg,'WeightedCentroid');
            
            OA_total_offset_x = (OA_c+OA_dir*OA_c_offset);
            OA_total_offset_y = (OA_r+OA_dir*OA_r_offset);
            dist_x = xpp(1,:)-OA_total_offset_x;
            dist_y = xpp(2,:)-OA_total_offset_y;
            rho_pp = sqrt(dist_x.^2+dist_y.^2);
            rho_pp_d = rho_pp+offaxis_distortion_quadratic*rho_pp.^2+...
                offaxis_distortion*rho_pp.^3;
            im_th = atan2(dist_y,dist_x);
	
            %X_p = im_rhop.*cos(im_th);
            undist_x = rho_pp_d.*cos(im_th)+OA_total_offset_x;
            undist_y = rho_pp_d.*sin(im_th)+OA_total_offset_y;
            xpp = cat(1,undist_x,undist_y);
        end
            
        ul_saved(:,ind) = c;
        nx_saved(ind) = nlinesx;
        ny_saved(ind) = nlinesy;
        
        
        
        
        scatter(xpp(1,:)*im_downsize,xpp(2,:)*im_downsize,'w+')
        
        if vis_grid
            scatter(xgm_adj(1),xgm_adj(2),'k+')
            for m = 0:nlinesy
                x3(1:2,m+1) = c'+uy*m*dpy';%+nlinesx*ux*dpx;
            end

            x4 = x3+repmat((ux*nlinesx*dpx)',[1,size(x3,2)]);
            l1 =line([x3(1,:);x4(1,:)],[x3(2,:);x4(2,:)]);
            for mm = 1:numel(l1)
                l1(mm).Color = 'r';
            end


            for m = 0:nlinesx
                x5(1:2,m+1) = c'+ux*m*dpx';%+nlinesy*uy*dpy;
            end

            x6 = x5+repmat((uy*nlinesy*dpy)',[1,size(x5,2)]);
            l2 = line([x5(1,:);x6(1,:)],[x5(2,:);x6(2,:)]);
            for mm = 1:numel(l2)
                l2(mm).Color = 'r';
            end
            
        end
        axis(im_downsize*[ulc ulc+roi_w ulr ulr+roi_h])
        %axis([0 200 -20 70])
        %axis([30 500 -20 700])
        set(gca,'position',[0 0 1 .97],'units','normalized')
        axis off
        pause(1/5)
        
    end
end
%% Figure out index order for each image
xang = atan2(uxp(2),uxp(1));
rmat = [cosd(-xang) -sind(-xang); sind(-xang) cosd(-xang)];
ptsrot = rmat*ul_saved;
sub_est = round((ptsrot-repmat(rmat*xgm',[1,nshifts^2]))./...
    repmat([dpx/nshifts;dpy/nshifts],[1,nshifts^2]));
sub_est_z = sub_est+repmat([1-min(sub_est(1,:));1-min(sub_est(2,:))],[1,nshifts^2]);

%This vector holds the starting index for each sub image in row-column
%order (not x-y). i.e. if the entry 1 has index m corresponding to the mth
%image, then the mth image's upper-left-most point is the first entry in
%the list.
ind_est = sub2ind([nshifts,nshifts],sub_est_z(2,:),sub_est_z(1,:));

c0 = ul_saved(:,find(ind_est==1));
figure(5),clf
axis([c0(1) c0(1)+dpx c0(2) c0(2)+dpy])
hold on
[s si] = sort(ind_est);
for m = 1:nshifts^2
    scatter(ul_saved(1,si(m)),ul_saved(2,si(m)))
    pause(1/5)
end


%% Read in each diffuser image and crop
h4 = figure(4),clf
npts = sum((nx_saved+1).*(ny_saved+1));

[ny_total, nx_total] = find_sizes_from_ul(c0,ulc+roi_w,ulr+roi_h,dpx/nshifts,dpy/nshifts,uxp,uyp);
nrows = ny_total+1;
ncols = nx_total+1;
if npts~=(nrows*ncols)
    error('total rows and columns disagrees with total points')
end
%%
if ~precompute
    
    
    kernel = fspecial('Gaussian',65,15);
    kernel1 = fspecial('Gaussian',3511);
    
    SE = strel('disk',bg_size);    %background removal strel
    
    %Prepare cell arrays for sparse value storage.
  
    r_outc = cell(npts,1);
    c_outc = cell(npts,1);
    v_outc = cell(npts,1);
    
    for m = 1:length(im_list)
        %m = mm(n);
        set(0,'CurrentFigure',h4)
        k = im_list(m);
        im_in_diff = im_in_dist(:,:,m);
        
        %Filter image and subtract background
        if bg_size>1
            filtered1 = imfilter(im_in_diff,kernel1);

            bg = imopen(filtered1,SE);
            A_bgrm_full = double(im_in_diff-bg);
            A_bgrm_full(A_bgrm_full<0) = 0;
        else
            A_bgrm_full = double(im_in_diff);
            fprintf('using original\n')
        end
        A_bgrm = imresize(A_bgrm_full,im_downsize);
        A_bgrm = A_bgrm(ceil(im_downsize*ulr:im_downsize*(ulr+roi_h)-1),...
            ceil(im_downsize*ulc:im_downsize*(ulc+roi_w)-1));
        %imagesc(A_bgrm)
        
        ul = (ul_saved(:,m)-[ulc;ulr])*im_downsize;
        nlinesx = nx_saved(m);
        nlinesy = ny_saved(m);
        dpx_scl = im_downsize*dpx;
        dpy_scl = im_downsize*dpy;
        xpp = grid_intersections(ul',ux,uy,dpx_scl,dpy_scl,nlinesx,nlinesy);
        x = 1:size(A_bgrm,2);
        y = 1:size(A_bgrm,1);
%             im_x = (1:im_w)-OA_c;
%     im_y = (1:im_h)-OA_r;
%     [im_X,im_Y] = meshgrid(im_x,im_y);
%     im_rho = sqrt(im_X.^2+im_Y.^2);
%     im_rho = im_rho;
%     im_th = atan2(im_Y,im_X);
%     im_rhop = im_rho + max(im_rho(:))*d*(im_rho/max(im_rho(:))).^3;
%     X_p = im_rhop.*cos(im_th);
%     Y_p = im_rhop.*sin(im_th);
        if offaxis_distortion~=0
            im_x = x/im_downsize-OA_total_offset_x;
            im_y = y/im_downsize-OA_total_offset_y;
            [im_X,im_Y] = meshgrid(im_x,im_y);
            
            rho_pp = sqrt(im_X.^2+im_Y.^2);
            rho_pp_d = rho_pp-offaxis_distortion_quadratic*rho_pp.^2-...
                offaxis_distortion*rho_pp.^3;
            im_th = atan2(im_Y,im_X);
	
            %X_p = im_rhop.*cos(im_th);
            X = im_downsize*(rho_pp_d.*cos(im_th)+OA_total_offset_x);
            Y = im_downsize*(rho_pp_d.*sin(im_th)+OA_total_offset_y);
            
        else
            [X,Y] = meshgrid(x,y);
        end
        
        
        for n = 1:size(xpp,2)
            [r,c] = ind2sub([nlinesy+1,nlinesx+1],n);
            clf
            
            xmask = X-xpp(1,n);
            ymask = Y-xpp(2,n);
            
            pow = 10;
            supergauss = exp(-(xmask.^pow+ymask.^pow)/((dpx_scl/1.96).^pow));
            maskn = (abs(xmask)<=dpx_scl/1.5) & (abs(ymask)<=dpy_scl/1.5);
            maskn = maskn.*supergauss;
        
            masked = maskn.*A_bgrm;
%             imagesc(masked)
%             axis image
%             pause(1/100)
            maskn = zeros(size(maskn));
            colormap default
            %Row and col offset based on ordering
            r_off = sub_est_z(2,m);
            c_off = sub_est_z(1,m);
            r_adj = nshifts*(r-1)+r_off;
            c_adj = nshifts*(c-1)+c_off;
            idx = sub2ind([nrows,ncols],r_adj,c_adj);
            [r_outc{idx}, c_outc{idx}, v_outc{idx}] = ...
                build_A_matrix_sparse(masked,idx);
            fprintf('Done with number %i from image %i. Index: %i \n',n,m,idx)
        end
        clf
        imagesc(A_bgrm)
        hold on
        
        scatter(xpp(1,:),xpp(2,:),'w+')
        axis image
        grid on
        set(gca,'YDir','reverse')
        % axis([500 1000 500 1000])
        set(gca,'position',[0 0 1 1],'units','normalized')
        hold off
        pause(1/2)
    end
    
    r_out = vertcat(r_outc{:});
    c_out = vertcat(c_outc{:});
    v_out = vertcat(v_outc{:});
    A = sparse(r_out,c_out,v_out,numel(A_bgrm),npts);
else
    A_in = load('A_first_exp.mat');
    A = A_in.A;
end

%%
w = 101;
h = 101;
cman_c = 128;
cman_r = 80;

im4= imread('cameraman.tif');

im4 = im4(round(cman_r-h/2):round(cman_r+h/2)-1,round(cman_c-w/2):round(cman_c+w/2)-1);
im4 = imrotate(im4,-90);
figure(4),clf
imagesc(imrotate(im4,90))
set(gca,'position',[0 0 1 .97],'units','normalized')
axis image
colormap gray
title('Input image')

%Correct gamma
lf = double(im4(:)).^(2.2);
sensor = A*double(lf);
sensor_reshaped = reshape(sensor,floor((s+1)/sub),floor((s+1)/sub));

figure(5),clf
imagesc(imrotate(sensor_reshaped,90))
set(gca,'position',[0 0 1 .97],'units','normalized')
axis image
title('Predicted sensor image using forward model, gamma correction applied')
colormap default
caxis([0 prctile(sensor,99.8)])

%% Invert
invert_in = demosaic(imread('/Users/nick.antipa/Documents/Diffusers/20150601/cman_diff_small_aperture.tif'),'bggr');
%sub = 5;
kernel_sm = ones(sub);   %Smoothing kernel
coff = 0;
roff = -3;
to_invert_g = double(invert_in(ulr+roff:ulr+s+roff,ulc+coff:ulc+s+coff,2));
to_invert_sm = imfilter(to_invert_g,kernel_sm,'symmetric');
to_invert = to_invert_sm(1:sub:end-sub,1:sub:end-sub);

%to_invert = to_invert-prctile(to_invert(:),.05);
to_invert = to_invert-min(min(to_invert));

%to_invert = to_invert-mean2(to_invert);
%to_invert1 = sensor_reshaped(1:sub:end,1:sub:end);
% kernel_mat_sm = ones(sub,1);
A_sm = A;
%b = zeros(s+1,s+1);
%b(1:sub:end,1:sub:end) = 1;
%idx = find(b);
figure(2)
imagesc(imrotate(to_invert,90));
set(gca,'position',[0 0 1 .97],'units','normalized')
axis image

title('Measured sensor image')
caxis([prctile(to_invert(:),.1) prctile(to_invert(:),99.8)])
%caxis([0 2.4e4])
%axis([100 170 0 50])
%Subsample A
%A_sub = A_sm(idx,:);
A_sub = A_sm;


AtA = A_sub'*A_sub;



lambda = 1e5;

AtA_r = (AtA+lambda*speye(size(AtA)));

recovered = abs(AtA_r\(A_sub'*(to_invert(:))));

figure(3)
imagesc(imrotate(abs(reshape(recovered,size(xs,1),size(xs,2)).^(1/2.2)),90))
caxis([prctile(recovered.^(1/2.2),.1) prctile(recovered.^(1/2.2),100)])
set(gca,'position',[0 0 1 .97],'units','normalized')
title('recovered image')
axis image
caxis([prctile(abs(recovered).^(1/2.2),.1) prctile(abs(recovered).^(1/2.2),99.9)])

colormap(gray)




