ntheta = 7;
nphi = 7;
ims_per_theta = 1;
ims_per_phi = 1;
ntheta_tot = ims_per_theta*ntheta;
nphi_tot = ims_per_phi*nphi;
nx = 141;
ny = 141;
upsamp = 1;
% n_y = 380;
% n_x = 380;
xstart = 1;
ystart = 1;
cols = floor(xstart*upsamp):floor(xstart*upsamp)+nx-1;
rows = floor(ystart*upsamp):floor(ystart*upsamp)+ny-1;
clf
save_gif=0;

figure(5),clf
count = 1;
a = zeros(nx,ny,nphi,ntheta);
r = a;
g = a;
b = a;
a_test = a;
subim_type = 'add';
monochrome = 0;  %1 for red, 2 green 3 blue
h8 = figure(8),clf
fprintf('blah')
for q = 1:nphi
    for p = 1:ntheta
        count = count+1;
        
        
        if monochrome
            im_in = imread(['/Users/nick.antipa/Documents/Light_field_data/Dragons/dragons-',...
                num2str(count,'%03d'),'.png']);
            a(:,:,q,p) = im_in(rows,cols,monochrome);
            imagesc(a(:,:,q,p))
            pause(1/24)
            %a_test(:,:,q,p) = reshape(count*nx*ny:count*nx*ny+nx*ny-1,[ny,nx]);
        else
            cint = mod(q-1,nphi)*ntheta+p;
            im_in = imread(['/Users/nick.antipa/Dropbox/PovRay/Scripts/SpheresAndBlock/01_spheres-',...
                num2str(cint,'%02d'),'.png']);
            %im_in = imread(['/Users/nick.antipa/Documents/Light_field_data/DragonBunnyFlatcam/dragons-flatcam-',...
                    %num2str(cint,'%03d'),'.png']);
            im_in = imresize(im_in,upsamp);
            r(:,:,q,p) = im_in(rows,cols,1);
            g(:,:,q,p) = im_in(rows,cols,2);
            b(:,:,q,p) = im_in(rows,cols,3);
            set(0,'CurrentFigure',h8)
            clf
            imagesc(uint16(round(cat(3,r(:,:,q,p),g(:,:,q,p),b(:,:,q,p)))))
            axis image
            drawnow
           
        end
    end
end




%str = input('Output mat file name: ','s')
str = '1_degree_lf'
filename = ['/Users/nick.antipa/Documents/Light_field_data/spheresAndBlock/1_degree/',str,'.mat']
if monochrome
    lf = permute(a,[4,3,2,1]);
else
    lfr = permute(r,[4,3,2,1]);
    lfg = permute(g,[4,3,2,1]);
    lfb = permute(b,[4,3,2,1]);
    save(filename,'lfr','lfg','lfb');
end

%%
figure()
for n = 1:nx
    for m = 1:ny
        if monochrome
            lf_im((m-1)*nphi+1:m*nphi,(n-1)*ntheta+1:n*ntheta) = (lf(:,:,n,m));
        else
            for p = 1:3
                if p == 1
                    lf_im((m-1)*nphi+1:m*nphi,(n-1)*ntheta+1:n*ntheta,p) = uint8(lfr(:,:,n,m));
                elseif p == 2
                    lf_im((m-1)*nphi+1:m*nphi,(n-1)*ntheta+1:n*ntheta,p) = uint8(lfg(:,:,n,m));
                else
                    lf_im((m-1)*nphi+1:m*nphi,(n-1)*ntheta+1:n*ntheta,p) = uint8(lfb(:,:,n,m));
                end
            end
        end
    end
end
imagesc(uint8(lf_im))
%axis([165 220 115 170])
%%
index = 0;
for mm = 1:ny
    for nn = 1:nx
        for pp = 1:nphi
            for qq = 1:ntheta
                index = index+1;
                lf_vec(index) = lf(pp,qq,nn,mm);
            end
        end
    end
end
%%
if save_gif
    str = input('Output file name: ','s')
    filename = ['../Output/',str,'.gif']
end

%%
count = 0;
set(0,'CurrentFigure',h8)
for n = 1:nphi
    for m = 1:ntheta
        count = count+1;
        im_int = (lf_im(n:5:end,m:5:end,:));
        im_n = double(im_int) + 255*.001*randn(size(im_int));
        imagesc(uint8(im_n))
        
        [imind,cm] = rgb2ind(uint8(lf_im(n:5:end,m:5:end,:)),256);
        if save_gif
            if count == 1;
                imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
            else
                imwrite(imind,cm,filename,'gif','WriteMode','append');
            end
        end
        axis image
        pause(1/15)
    end
end




