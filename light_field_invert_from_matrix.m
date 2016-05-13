% load('/Users/nick.antipa/Documents/MATLAB/Output/1deg_20151121_A_sub_128x128x5x5.mat')
% AtA = A_sub'*A_sub;
%im_crop(s_resh==0)=0;
%linkaxes([a1,a9],'xy')

%c = mean(sum(A_sub.^2,1))
A = A_sub/(1.5e5);
clear A_sub;
At = A';
% corred = corr2_fft(im_crop,s_resh);
% figure(99),
% contour(corred)
% axis image
% axis([1020 1028 1020 1028])
% grid on

%%
input_data = 'real'
downsamp = 1; % Binning integer. 2 uses every other, 3 uses every third etc.
monochrome = 1;
preload_matrices = 1;
save_file= 0;
%nsx = settings.npx;    %Number of sensor pixels in x
%nsy = settings.npy;  %Number of sensor pixels in y
%nx = settings.N;   %Number of light field x bins
%ny = settings.M;   %Number of light field y bins
%ntheta = settings.P; %Number of theta bins

%nphi = settings.Q; %Number of phi bins
switch lower(input_data)
    case('real')
        figure(1)
        a1 = gca;
        crop = 0;
        dist = 0;
        circ_aper = 0;
        shift = [0,1];
        if crop
            %crop_c = 1280;
            %crop_r = 1080;
            %crop_half_size = 2^10;
            crop_c = 1030;
            crop_r = 761;
            crop_half_size = 501;
            odd_size = 1;
        end
        
        
        im_in = imread('/Users/nick.antipa/Documents/Diffusers/1deg_yes4f_flea3_20160510/pics/4p5mmAper_z654/king_leaf1.tif');
        
        
        %in = load('../Output/1deg_20151129_A_sub_128_5_z260_gptiesensor_sim.mat');
        %im_in = in.s_resh;
        
        
        figure(9)
        
        if dist
            OA_c = crop_c-512+1;
            OA_r = crop_r-512+1;
            c = 1/.9955;
            r_off = -2.1;
            c_off = -(-1.58);
            im_x = (1:size(im_in,2))-OA_c;
            im_y = (1:size(im_in,1))-OA_r;
            [im_X,im_Y] = meshgrid(im_x/512,im_y/512);
            im_rho = sqrt(im_X.^2+im_Y.^2);
            im_th = atan2(im_Y,im_X);
            im_rhop = im_rho*c;
            X_p = im_rhop.*cos(im_th);
            Y_p = im_rhop.*sin(im_th);
            for k = 1
                im_in_dist = interp2(im_X*512-c_off,im_Y*512-r_off,double(im_in),X_p*512,Y_p*512);
                im_dist = im_in_dist;
            end
            im_dist(isnan(im_dist))=0;
        else
            im_dist = double(im_in);
        end
        
        if crop
            if odd_size ~=1 
                im_crop = double(im_dist(crop_r-crop_half_size:crop_r+crop_half_size-1,crop_c-crop_half_size:crop_c+crop_half_size-1));
            else
                im_crop = double(im_dist(crop_r-crop_half_size+1:crop_r+crop_half_size-1,crop_c-crop_half_size+1:crop_c+crop_half_size-1));
            end
        else
            im_crop = im_dist;
        end
 
        im_crop = circshift(im_crop,shift);

        %im_crop = im_crop-min(im_crop(:));
        %im_crop = im_crop-min(min(im_crop));
        
        if downsamp ~= 1
            if (downsamp-floor(downsamp))~=0
                error('You need to use an integer for downsampling')
            end
            im_conv = conv2(im_crop,ones(downsamp),'same');
            im_crop = im_conv(1:downsamp:end,1:downsamp:end);
            
        end
        colormap default
        imagesc_pctl(im_crop,5,99.9);
        axis image
        a9 = gca
        b = double(im_crop(:));
    case('simulation')
        if ~preload_matrices
            %in_inv = load('/Users/nick.antipa/Documents/MATLAB/Output/1deg_20151129_A_sub_128_5_z260_gptie.mat');
            in_inv = load('/Users/nick.antipa/Documents/MATLAB/Output/1deg_20151129_A_sub_512_5_z260_gptie.mat');
            %inv_str = '../Output/A_sub_38_15_600_for_square.mat';
            %inv_str = '../Output/A_sub_38_15_425_forCS_smallang.mat';
            in_inv = load(inv_str);
            A_sub = in_inv.A_sub;
            settings = in_inv.settings;
            clear in_inv
            in_for = load(inv_str);
            A_for = in_for.A_sub;
            clear in_for
            
            if size(A_for,1)~=size(A_sub,1)
                error('matrix sensor sizes do not match')
            end
            
            A = A_sub/settings.nrays;
            At = A';
            %clear A_for
            %clear A_sub;
        end
        lf_str = '../Output/38_15_4d_lf.mat';
        lf_in = load(lf_str);
        if monochrome
            b_raw = A*double(lf_in.lfg(:))*settings.nrays*12;
        end
        b = b_raw+randn(size(b_raw))*.01*prctile(b_raw,99);
        b(b<0) = 0;
        b = b*1e-9;
%         nsx = settings.npx;    %Number of sensor pixels in x
%         nsy = settings.npy;  %Number of sensor pixels in y
%         nx = settings.N;   %Number of light field x bins
%         ny = settings.M;   %Number of light field y bins
%         ntheta = settings.P; %Number of theta bins
%         
%         nphi = settings.Q; %Number of phi bins
        
end


b = b-min(b);
b(b<0)=0;
nphi = settings.Q; %Number of phi bins
clear im_dist;
clear im_in_dist;
clear im_rho;
clear im_rhop;
clear im_th;
clear im_X;
clear im_Y;
%%

nsx = settings.npx;    %Number of sensor pixels in x
nsy = settings.npy;  %Number of sensor pixels in y
nx = settings.N;   %Number of light field x bins
ny = settings.M;   %Number of light field y bins
ntheta = settings.P; %Number of theta bins
count= 0;
order_vec = [];
%prepares vector for reshaping into stack
for m = 1:nphi
    for n = 1:ntheta
        count = count+1;
        order_vec = cat(2,order_vec,count:ntheta*nphi:nx*ny*ntheta*nphi);
    end
end

count = 0;
deindex_vec = zeros(nx*ny*ntheta*nphi,1);
%Takes stack to vector
for q = 1:nx
    for p = 1:ny
        count = count+1;
        deindex_vec((count-1)*(ntheta*nphi)+1:(count)*(ntheta*nphi)) = count:nx*ny:nx*ny*ntheta*nphi;
    end
end

angle_deindex_vec = (1:nx*ny*ntheta*nphi)';   %Transform stack to vector
angle_index_vec = (1:nx*ny*ntheta*nphi)';  %vector to stack, switch nx,ntheta and ny,nphi
% Make indexing to go from angle stack to vector (for input to A)
%x0 = A_adj_tv(b,At,order_vec,nx,ny,ntheta,nphi);
%bp = A_tv(x0,A,deindex_vec);



lamb = 0;
write_gif = 0;
use_twist = 1;
if ~use_twist
    recovered = (AtA+lamb*speye (size(AtA)))\(A_sub'*double(im_crop(:)));
elseif use_twist
    %Phi_func = @(x)(1e8*x'*-(x<0));
    %A_func = @(x) A_tv(x,A,angle_deindex_vec);
    model_type = '3dtv';
    solver_type = 'TwIST';
    switch lower(model_type)
        case 'angle_tv'
            A_func = @(x) A_tv(x,A,angle_deindex_vec);
            A_adj_func = @(x) A_adj_tv(x,At,angle_index_vec,ntheta,nphi,nx,ny);
            Phi_func = @(x) TVnorm_lf(x);
            Psi_func = @(x,th)  tvdenoise_lf(x,2/th,tv_iters);
            debias = 0;
        case '3dtv'
            tv_iters = 3;
            A_func = @(x) A_tv(x,A,deindex_vec);
            A_adj_func = @(x) A_adj_tv(x,At,order_vec,nx,ny,ntheta,nphi);
            Phi_func = @(x) TVnorm3d(x);
            Psi_func = @(x,th) tvdenoise3d(x,2/th,tv_iters,1);
            debias = 0;
            tau = 10;
        case 'wavelet'
            dwtmode('per');
            wvlt = 'db9';
            dwvlt = 'db9';
            %wvlt = 'bior6.8';
            %dwvlt = 'rbio6.8';
            N = 4;
            [C,L] = wavedec2(rand(nx,ny),N,wvlt);
            wvltsize = length(C);
            %             for m = 1:ntheta*nphi
            %                 x = cat(1,x,wavedec2(recovered(:,:,m)',N,wvlt)');
            %             end
            A_func = @(x) A_wvlt2d(x,L,dwvlt,A,nx,ny,ntheta,nphi,deindex_vec);
            A_adj_func = @(x) A_adj_wvlt2d(x,At,wvlt,N,order_vec,ntheta,nphi,nx,ny,wvltsize);
            Phi_func = @(x) norm(x,1);
            Psi_func = @(x,tau) soft(x,tau);
            debias = 0;
            tau = 5;
        case 'hard_wavelet'
            dwtmode('per');
            %wvlt = 'db9';
            %dwvlt = 'db9';
            wvlt = 'bior6.8';
            dwvlt = 'rbio6.8';
            N = 4;
            [C,L] = wavedec2(rand(nx,ny),N,wvlt);
            wvltsize = length(C);
            %             for m = 1:ntheta*nphi
            %                 x = cat(1,x,wavedec2(recovered(:,:,m)',N,wvlt)');
            %             end
            A_func = @(x) A_wvlt2d(x,L,dwvlt,A,nx,ny,ntheta,nphi,deindex_vec);
            A_adj_func = @(x) A_adj_wvlt2d(x,At,wvlt,N,order_vec,ntheta,nphi,nx,ny,wvltsize);
            Phi_func = @(x) norm(x,1);
            Psi_func = @(x,tau) hard(x,tau);
            debias = 0;
            tau = 9000;
        case 'angle_wavelet'
            wvlt = 'haar';
            dwvlt = 'haar';
            N = 2;
            [C,L] = wavedec2(zeros(nphi,ntheta),N,wvlt);
            wvltsize = length(C);
            A_func = @(x) A_wvlt2d(x,L,dwvlt,A,ntheta,nphi,nx,ny,angle_deindex_vec);
            A_adj_func = @(x) A_adj_wvlt2d(x,At,wvlt,N,angle_index_vec,nx,ny,ntheta,nphi,wvltsize);
            Phi_func = @(x) norm(x,1);
            Psi_func = @(x,tau) soft(x,tau);
            debais = 0;
            tau = 30000;
        case 'lsmr'
            lam = 1e-10;
            
    end
end
%A_adj_func = @(x) A_adj_tv(x,At,angle_index_vec,ntheta,nphi,nx,ny);
%Phi_func = @(x) norm(x,1);
%Phi_func = @(x) TVnorm_lf(x);


% denoising function;

%Psi_func = @(x,th)  tvdenoise_lf(x,2/th,tv_iters);

%Psi_func = @(x,th) tvdenoise_lf(x,2/th,tv_iters);
%Phi_func = @(x) TVnorm_lf(x,ntheta,nphi,ny,nx);
%Phi_func = @(x) wave_norm(x,100,'db9');
%Psi_func = @(x,tau) soft_nn(x,tau);
%Psi_func = @(x,tau) nothing(x,tau);
%Phi_func = 6
switch lower(solver_type)
    case('twist')
        recovered =  TwIST(b,A_func,tau,...
            'AT', A_adj_func, ...
            'StopCriterion',1,...
            'Psi',Psi_func,...
            'Phi',Phi_func,...
            'Initialization',2,...
            'Debias',debias,...
            'MaxiterA',300,...
            'Monotone',1,...
            'ToleranceA',.0001);
    case('lsmr')
        recovered = lsmr(A,b,lam,[],[],[],[],[],1e6);
end
%20151205 These work really well with the 180x180x7x7x2048 4 degree
%matrix along with TVnorm3d and tvdenoise3d, tviters 5, A_tv and
%A_adj_tv with A=A_sub/sqrt(9e7) and At=A'
%     recovered =  TwIST(b,A_func,100,...
%         'AT', A_adj_func, ...
%         'StopCriterion',1,...
%         'Psi',Psi_func,...
%         'Phi',Phi_func,...
%         'Initialization',2,...
%         'Debias',0,...
%         'MaxiterA',1000,...
%         'Monotone',1,...
%         'ToleranceA',.0001);
rec_backup = recovered;
%%
recovered = rec_backup;
switch lower(model_type)
    case 'angle_tv'
        recovered_ang = recovered;
        recovered_vec = recovered_ang(:);
        recovered = reshape(recovered_vec(order_vec),[nx,ny,ntheta*nphi]);
    case 'angle_wavelet'
        recovered_ang = recovered;
        recovered_vec = recovered_ang(:);
        recovered = reshape(recovered_vec(order_vec),[nx,ny,ntheta*nphi]);
    case {'wavelet','hard_wavelet'}
        rec_w = zeros(nx,ny,ntheta*nphi);
        
        decsize = length(recovered)/ntheta/nphi;
        for m = 1:ntheta*nphi
            
            
            vind = (m-1)*decsize+1:m*decsize;
            rec_w(:,:,m) = waverec2(recovered(vind),L,wvlt);
        end
        recovered = rec_w;
    case 'lsmr'
        recovered = reshape(rec_backup(order_vec),ny,nx,ntheta*nphi);
end

%%
save_subap = 0;
if save_subap
    out_filename = input('Filename (type ''n'' to cancel) ','s')
    if strcmpi(out_filename,'n')
        error('cancelled')
    end
    dots = findstr(out_filename,'.');
    if isempty(dots)
        out_filename = [out_filename,'.gif'];
    end
    
end
for n =1:ntheta*nphi
    %     if ~ismember(n,good_idx)
    %         recovered(:,:,n)=0;
    %     end
    pct = prctile(recovered(:),99);
    subap_scl = recovered(:,:,n)'/pct*255;
    imagesc(subap_scl);
    colormap gray
    axis image
    caxis([prctile(subap_scl(:),1) 255])
    pause(1/10)
    drawnow
    if save_subap
        
        [imind,cm] = rgb2ind(repmat(uint8(subap_scl),[1,1,3]),256);
        if n == 1;
            imwrite(imind,cm,out_filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,out_filename,'gif','WriteMode','append');
        end
    end
end

%% Make subaperture layout
count = 0;
layout = zeros(nphi*ny,ntheta*nphi);
for n = 1:ntheta
    for m = 1:nphi
        count = count+1;
        layout((n-1)*ny+1:n*ny,(m-1)*nx+1:m*nx)=recovered(:,:,count)';
    end
end
imagesc(layout/prctile(layout(:),99.9))
axis image
caxis([0 1])

%% Make lenslet view
angle_index = zeros(size(A,2),1);  %Maps from subaperture stack to angle stack
count = 0;
for n = 1:ny
    for m = 1:nx
        count = count+1;
        angle_index((count-1)*ntheta*nphi+1:count*ntheta*nphi) = count:nx*ny:ntheta*nx*nphi*ny;
    end
end

%% display

h7 = figure(7);
set(0,'CurrentFigure',h7)

%Loop over theta and phi, display LF(x,y) at each (theta,phi)
count = 0
index_vector = [];
for m = 1:nphi
    for n = 1:ntheta
        count = count+1;
        index_vector = cat(2,index_vector,[count:ntheta*nphi:length(recovered)]);
        imagesc(reshape(recovered(count:ntheta*nphi:end),nx,ny)')
        %caxis([0 35])
        colormap gray
        axis image
        caxis([0 10])
        pause(1/9)
    end
    
end

%% make xtheta

xtheta = zeros(ntheta,nx);
m = 6
%if m == ceil(nphi/2)  %Make x theta slice
for nn = 1:ntheta
    ii = sub2ind([nphi,ntheta],nn,m)
    xtheta(nn,:) = wvlt_recovered.recovered(:,60,ii);
end
%end

imagesc(xtheta)

%% Try refocusing. Because why not.
h8 = figure(8)


imup = 1;
im_x = 1:nx;
im_y = 1:ny;
[IMX, IMY] = meshgrid(im_x,im_y);
save_gif = 0

if save_gif
    out_filename = input('Filename (type ''n'' to cancel) ','s')
    if strcmpi(out_filename,'n')
        error('cancelled')
    end
    dots = findstr(out_filename,'.');
    if isempty(dots)
        out_filename = [out_filename,'.gif'];
    end
    
end
fcount = 0;
xtheta = zeros(ntheta,nx);
scale_vec = -1
refocused_stack= zeros(settings.M,settings.N,length(scale_vec));
for s = scale_vec
    fcount = fcount+1;
    %s =
    refocused = zeros(nx*imup,ny*imup);
    count = 0;
    for m = 1:nphi
        for n = 1:ntheta
            count = count+1;
            
            %This deals with the fact that the angles of matrices are
            %oposite that of the POVRay simulations. It's a hack.
            switch lower(input_data)
                case('simulation')
                    xshift = -(n-ntheta/2-1/2);
                case('real')
                    xshift = (n-ntheta/2-1/2);
            end
            yshift = (m-nphi/2-1/2);
            r = sqrt(xshift^2+yshift^2);
            %amt = .25; use with r^2
            amt = 0;
            rho = r/ntheta*2;
            ang = atan2(yshift,xshift);
            %xc = (rho*cos(ang))^3*amt-amt*ntheta/2*xshift
            %yc = (rho*sin(ang))^3*amt+amt*nphi/2*yshift;
            xs = xshift/(ntheta-1)*2;
            ys = yshift/(nphi-1)*2;
            xc = ((xs)^3-xs)*amt;
            yc = ((ys)^3-ys)*amt;
            
            %sub_im = transpose(imresize(reshape(recovered(count:ntheta*nphi:end),ny,nx),imup));
            sub_im = recovered(:,:,count)';
            %im_shifted = circshift(sub_im,[-round(s*(m-nphi/2-1/2)),-round(s*(n-ntheta/2-1/2))]);
            im_shifted = interp2(IMX,IMY,sub_im,IMX+s*xshift+xc,IMY+s*yshift+yc,'nearest');
            %im_shifted(isnan(im_shifted)) = mean2(sub_im);
            im_shifted(isnan(im_shifted))=0;
            %round(s*(m-nphi/2-1/2))
            %round(s*(n-ntheta/2-1/2));
            refocused = im_shifted+refocused;
            refocused(refocused<0) = 0;
            %             if m == ceil(nphi/2)  %Make x theta slice
            %                     ii = sub2ind([nphi,ntheta],nn,m);
            %                     xtheta(n,:) = mean(im_shifted(80,:),1);
            %                     xtheta(n,:) = xtheta(n,:)/2+circshift(xtheta(n,:),[0,1])/2;
            %             end
            %             %             %set(0,'CurrentFigure',h8)
            %                         imagesc(im_shifted)
            %                         colormap gray
            %                         caxis([0 1e8])
            %                         drawnow
            %             pause(1/5)
        end
        
    end
    
    refocused_scl = real((refocused/prctile(refocused(:),99.8)).^(1))*255;
    refocused_stack(:,:,fcount) = refocused_scl;
    imagesc(refocused_scl);
    colormap gray
    
    axis image
    caxis([0 255])
    pause(1/20)
    drawnow
    if save_gif
        
        %         frame = getframe(h8);
        %         im = frame2im(frame);
        [imind,cm] = rgb2ind(repmat(uint8(refocused_scl),[1,1,3]),256);
        if fcount == 1;
            imwrite(imind,cm,out_filename,'gif', 'Loopcount',inf);
        else
            imwrite(imind,cm,out_filename,'gif','WriteMode','append');
        end
    end
end
% figure(15),clf
% imagesc(xtheta,'XData',(1:size(im_shifted,2))*px,'YData',linspace(-1/2,1/2,ntheta)*settings.th_range)
% ylabel('\theta degrees')
% caxis([5000 2e4])
% xlabel('x \mum')
% colormap gray