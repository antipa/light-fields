%%
%Find center of each screen pixel on camera from diffuser-free image. These
%centers will serve as the centers for cropping each diffuser psf.
precompute = 0;
fname = '/Users/nick.antipa/Documents/Diffusers/20150601/c.tif';
info = imfinfo(fname);
num_images = numel(info);
max_proj = zeros(info(1).Height,info(1).Width);
im_start = 1;
im_list = 1:2:50;
for m = 1:numel(im_list);
    k = im_list(m);
    im_in(:,:,k) = imread(fname, k, 'Info', info);
    %imagesc(im_in(:,:,k))
    %pause(1/100)
end
[max_proj,im_id] = max(im_in,[],3);
clear A;
max_proj_dem = demosaic(max_proj,'rggb');
max_proj_g = max_proj_dem(:,:,2);
roi = max_proj_g>max(max(max_proj_g))*.05;
se = strel('disk',15);
roi_c = imclose(roi,se);
roi_o = imopen(roi_c,se);

[UL_r, UL_c, s, mat] = findMaxSq(roi_o);
ulr = UL_r(1);
ulc = UL_c(1);

%%
kernel = fspecial('Gaussian',25,2);
filtered = imfilter(max_proj_g,kernel);
B = watershed(max(max(filtered))-filtered,8);
B_c = B(ulr:ulr+s,ulc:ulc+s);
rgb = label2rgb(B_c,'lines',[.5 .5 .5]);
imshow(rgb);


%% clean up ROI
edges = zeros(size(B_c));
edges(1:end,1) = 1;
edges(1,:) = 1;
edges(end,:) = 1;
edges(:,end) = 1;
touch_edge = edges.*double(B_c);
to_remove = unique(touch_edge);
for m = 1:numel(to_remove)
    B_c(B_c==to_remove(m)) = 0;
end

%% background subtract
se1 = strel('disk',15);
bg = imopen(max_proj_g(ulr:ulr+s,ulc:ulc+s),se1);
max_bgrm = max_proj_g(ulr:ulr+s,ulc:ulc+s)-bg;
%%
p_stats = regionprops(B_c,max_bgrm,'WeightedCentroid');
imagesc(max_bgrm)
hold on
a = struct2array(p_stats);
a_good = a(~isnan(a));
scatter(a_good(1:2:end),a_good(2:2:end),'rx')
axis image
set(gca,'position',[0 0 1 1],'units','normalized')

%% Figure out grid
pad = 5;
spect = fftshift(fft2(max_bgrm-mean2(max_bgrm),5*size(max_bgrm,1),5*size(max_bgrm,2)));

x_roi = zeros(size(spect));
x = -floor(size(spect,2)/2):floor(size(spect,2))/2;
fx = linspace(-1/2/mean(diff(x)),1/2/mean(diff(x)),numel(x));


[FX, FY] = meshgrid(fx,fx);

%Find top 6 peaks within x roi
%%
peak_region = .02; %In frequency units
peak_reg_idx = round(peak_region/mean(diff(fx)));
x_roi = abs(FX)<.05 & abs(FY)>peak_region;
x_band = abs(spect).*x_roi;

pr = zeros(6,1);
pc = pr;
for n = 1:6
    mv = max(x_band(:));
    [pr(n), pc(n)] = find(x_band==mv,1,'first');
    x_band(pr(n)-peak_reg_idx:pr(n)+peak_reg_idx,pc(n)-peak_reg_idx:pc(n)+peak_reg_idx)=0;
    imagesc(x_band)
    caxis([0 1e7])
    pause(1/20)
end
xx = [0;mean(diff(fx))*(pc-max(x)-1)];
xy = [0;mean(diff(fx))*(pr-max(x)-1)];

%%
y_roi = abs(FY)<.05 & abs(FX)>peak_region;
y_band = abs(spect).*y_roi;

pr = zeros(6,1);
pc = pr;
for n = 1:6
    mv = max(y_band(:));
    [pr(n), pc(n)] = find(y_band==mv,1,'first');
    y_band(pr(n)-peak_reg_idx:pr(n)+peak_reg_idx,pc(n)-peak_reg_idx:pc(n)+peak_reg_idx)=0;
    imagesc(y_band)
    caxis([0 1e7])
    pause(1/10)
end
yx = [0;mean(diff(fx))*(pc-max(x)-1)];
yy = [0;mean(diff(fx))*(pr-max(x)-1)];

%% figure out rotation
clf
imagesc(abs(spect),'XData',fx,'YData',fx)
set(gca,'YDir','normal')
hold on
caxis([0 1e7])
scatter(yx,yy,'rx')
scatter(xx,xy,'k+')

my = xx'*xy/(xx'*xx);
ly = line(1/my*[min(fx),max(fx)],[min(fx), max(fx)]);
ly.Color = 'k'

mx = yx'*yy/(yx'*yx);
lx = line([min(fx), max(fx)],mx*[min(fx),max(fx)]);
lx.Color = 'r'
hold off


%% figure out period
ufx = [1;mx];
ufx = ufx/norm(ufx);

ufy = [1/my;1];
ufy = ufy/norm(ufy);

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

xg = a_good(1:2:end);

yg = a_good(2:2:end);

clf, imagesc(max_bgrm)
axis image
hold on
% l1 = line([xgm,xgm+1300],[ygm,ygm+1300*(-1/my)]);
% l1.Color = 'r'

ux = [1,-1/my];
ux = ux/norm(ux);

uy = [-mx,1]
uy = uy/norm(uy)

%Find point in center (to minimize distortion issues with fitting grid)
xc = xg-size(max_bgrm,2)/2;
yc = yg-size(max_bgrm,1)/2;
rc = sqrt(xc.^2+yc.^2);
[rm midx] = min(rc);
xgm = xg(midx);
ygm = yg(midx);

xgm_adj = [xgm,ygm]-ux*50*dpx
scatter(xgm_adj(1),xgm_adj(2),'k+')
for m = -50:50;
    x3(1:2,m+51) = (xgm_adj+uy*m*dpy)';
end

x4 = x3+repmat((ux*100*dpx)',[1,size(x3,2)]);
l1 =line([x3(1,:);x4(1,:)],[x3(2,:);x4(2,:)]);
for mm = 1:numel(l1)
    l1(mm).Color = 'r';
end

ygm_adj = [xgm,ygm]-uy*50*dpy;

for m = -50:50;
    x5(1:2,m+51) = (ygm_adj+ux*m*dpx)';
end

x6 = x5+repmat((uy*100*dpy)',[1,size(x5,2)]);
l2 = line([x5(1,:);x6(1,:)],[x5(2,:);x6(2,:)]);
for mm = 1:numel(l2)
    l2(mm).Color = 'r';
end

%% Make points in order at grid intersections
npx = floor(s/dpx)-1;
npy = floor(s/dpy)-1;
x0 = [xgm,ygm];
count = 0;
xs = [];
ys = [];
mmin = -floor(npx/2)+2;
nmin = -floor(npy/2);
for m = mmin:floor(npx/2)
    xp = x0+dpx*ux*m;
    for n = nmin:floor(npy/2)-2
        
        xpp = xp+dpy*uy*n;
        if (xpp(1)>0&&xpp(2)>0&&xpp(1)<=s&&xpp(2)<=s)
            
            
            %         scatter(xpp(1),xpp(2),'kx')
            %         pause(1/1000)
            xs(m-mmin+1,n-nmin+1) = xpp(1);
            ys(m-mmin+1,n-nmin+1) = xpp(2);
        end
    end
end

%We now have xs and ys as the centers of each pixel in raster order (row by
%column, starting from upper left);
%% Read in each diffuser image and crop

fname = '/Users/nick.antipa/Documents/Diffusers/20150601/c_diff.tif';
info = imfinfo(fname);
num_images = numel(info);
max_proj = zeros(info(1).Height,info(1).Width);
im_list = [1:2:13,14:2:50,51:2:89,90:2:124];
h4 = figure(4)

%Find m that corresponds to upper left corner

% for m = 1:length(im_list)
%     ystart = mod(7-m,8)+1;
%     xstart = mod(6+ceil((m)/8),8)+1;
%     if xstart == 1 && ystart == 1
%         mstart = m
%     end
% end
%
% mm = circshift(1:length(im_list),[0,-(mstart-1)])
sub = 10;
if ~precompute
    x = 1:floor(size(im_in,2)/sub);
    y = 1:floor(size(im_in,1)/sub);
    [X,Y] = meshgrid(x,y);
    maskn = zeros(floor(size(X)/sub));
    
    kernel = fspecial('Gaussian',65,15);
    kernel1 = fspecial('Gaussian',35,11);
    SE = strel('disk',75);    %background removal strel
    kernel_sm = ones(sub);   %Smoothing kernel
    
    %Prepare cell arrays for sparse value storage.
    r_outc = cell(numel(xs),1);
    c_outc = cell(numel(xs),1);
    v_outc = cell(numel(xs),1);
    
    for m = 1:length(im_list)
        %m = mm(n);
        set(0,'CurrentFigure',h4)
        k = im_list(m);
        im_in_diff = imread(fname, k, 'Info', info);
        A_dem = demosaic(im_in_diff,'bggr');
        A_g = A_dem(:,:,2);
        

        ystart = mod(7-m,8)+1;
        xstart = mod(6+ceil((m)/8),8)+1;
        xsub = xs(xstart:8:end,ystart:8:end);
        ysub = ys(xstart:8:end,ystart:8:end);
        xt = xsub'/sub;
        yt = ysub'/sub;
        filtered1 = imfilter(A_g,kernel1);
        bg = imopen(filtered1,SE);
        A_bgrm_full = double(A_g-bg);
        A_bgrm_full(A_bgrm_full<0) = 0;
        
        A_g_sm = imfilter(A_bgrm_full,kernel_sm);
        A_bgrm = A_g_sm(1:sub:end,1:sub:end);
        
        
        for n = 1:numel(xsub)
            clf
            xmask = (xt(n)+ulc/sub-X);
            ymask = (yt(n)+ulr/sub-Y);
            
            pow = 6;
            supergauss = exp(-(xmask.^pow+ymask.^pow)/((4*dpx/sub).^pow));
            maskn = (abs(xt(n)+ulc/sub-X)<=5*dpx/sub) & (abs(yt(n)+ulr/sub-Y)<=5*dpy/sub);
            maskn = maskn.*supergauss;
%             maskn(round(ulr+ysub(n)-4*dpy):round(ulr+ysub(n)+4*dpy),...
%                 round(ulc+xsub(n)-4*dpx):round(ulc+xsub(n)+4*dpx))=1;
            
%             imagesc(maskn.*A_bgrm);
%             
%             caxis([0 3000])
%             
%             hold on
%             scatter(xt(n)+ulc/sub,yt(n)+ulr/sub,'k+')
%             pause(1/10000)
            masked = maskn.*A_bgrm;
            masked = masked(ceil(ulr/sub):floor((ulr+s)/sub),ceil(ulc/sub):floor((ulc+s)/sub));
            
            maskn = zeros(size(maskn));
            row = ystart+8*mod(n-1,size(xsub,2));
            col = xstart+8*floor((n-1)/size(xsub,2));
            idx = sub2ind(size(xs),row,col);
            [r_outc{idx}, c_outc{idx}, v_outc{idx}] = ...
                build_A_matrix_sparse(masked,idx);
            fprintf('Done with number %i from image %i\n',n,m)
        end
        clf
        imagesc(im_in_diff)
        hold on
        
        scatter(ulc+xsub(:),ulr+ysub(:),'rx')
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
    A = sparse(r_out,c_out,v_out,(floor((s+1)/sub))^2,numel(xs));
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
%%
%% %Segment each image individually
% segmented = cell(1,numel(im_start:2:num_images));
% stats_gray = cell(size(segmented));
% kernel = fspecial('Gaussian',65,15);
% kernel1 = fspecial('Gaussian',35,11);
% SE = strel('disk',75);
% count = 0;
% for k = im_start:2:num_images
%     count = count+1;
%     A = imread(fname, k, 'Info', info);
%     A_dem = demosaic(A,'rggb');
%     A_g = double(A_dem(ulr:ulr+s,ulc:ulc+s,2));
%     filtered = imfilter(A_g,kernel);
%     filtered1 = imfilter(A_g,kernel1);
%     bg = imopen(filtered1,SE);
%     segmented{count} = watershed(max(max(filtered))-filtered+eps,8);
%     edges = zeros(size(A_g));
%
%     edges(1:end,1) = 1;
%     edges(1,:) = 1;
%     edges(end,:) = 1;
%     edges(:,end) = 1;
%     touch_edge = edges.*double(segmented{count});
%     to_remove = unique(touch_edge);
%     for m = 1:numel(to_remove)
%         segmented{count}(segmented{count}==to_remove(m)) = 0;
%     end
%     rgb = label2rgb(segmented{count},'colorcube',[.5,.5,.5]);
%     stats_gray{count} = regionprops(segmented{count},A_g-bg,'WeightedCentroid');
%     imagesc(A_g)
%     colormap(gray)
%     axis image
%     hold on
%     for n = 1:numel(stats_gray{count})
%         scatter(stats_gray{count}(n).WeightedCentroid(1),stats_gray{count}(n).WeightedCentroid(2))
%     end
%     pause(1/100)
% end




%%
