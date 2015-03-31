ntheta = 5;
nphi = 5;
nx = 50;
ny = 50;
xstart = 740;
ystart = 210;

cols = xstart:xstart+nx-1;
rows = ystart:ystart+ny-1;
clf

lf = zeros(ntheta,npix,3);
figure(5),clf
count = 0;
a = zeros(nx,ny,nphi,ntheta);
a_test = a;
for p = 1:nphi
    for q = 1:ntheta
        count = count+1;
    %     a = imread(['/Users/nick.antipa/Documents/Light field/Data/dice/dice-',...
    %         num2str(n,'%02d'),'.png']);
        %a = imread(['/Users/nick.antipa/Documents/Light field/Data/xyzrgb_dragon/xyzrgb_dragon-',...
            %num2str(n,'%02d'),'.png']);
         im_in = imread(['/Users/nick.antipa/Documents/Light field/Data/DragonAndBunnies/DragonsAndBunnies_5x5_ap6.6/dragons-',...
             num2str(count,'%02d'),'.png']);
         a(:,:,q,p) = im_in(rows,cols,3);
        %lf(n,:,:) = a(row,cols,:);
        imagesc(a(:,:,q,p))
        %a_test(:,:,q,p) = reshape(count*nx*ny:count*nx*ny+nx*ny-1,[ny,nx]);        
    end
end

%%
lf = permute(a,[4,3,2,1]);
figure()
for n = 1:nx
    for m = 1:ny
        lf_im((m-1)*nphi+1:m*nphi,(n-1)*ntheta+1:n*ntheta) = (lf(:,:,n,m));
    end
end
imagesc(lf_im)
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
count = 0;
filename = 'dragon_bunny_animation_original.gif';
for n = 1:nphi
    for m = 1:ntheta
        count = count+1;
        imagesc(uint8(lf_im_rgb(n:5:end,m:5:end,:)))
        [imind,cm] = rgb2ind(uint8(lf_im_rgb(n:5:end,m:5:end,:)),256);
          if count == 1;
              imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
          else
              imwrite(imind,cm,filename,'gif','WriteMode','append');
          end
        axis image
        pause(1/7)
    end
end




