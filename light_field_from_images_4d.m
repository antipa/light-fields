ntheta = 5;
nphi = 5;
nx = 50;
ny = 50;
xstart = 740;
ystart = 210;
cols = xstart:xstart+nx-1;
rows = ystart:ystart+ny-1;
clf
save_gif=0;

figure(5),clf
count = 0;
a = zeros(nx,ny,nphi,ntheta);
r = a;
g = a;
b = a;
a_test = a;

monochrome = 0;  %1 for red, 2 green 3 blue
h8 = figure(8),clf
for p = 1:nphi
    for q = 1:ntheta
        count = count+1;
        if monochrome
        
        %     a = imread(['/Users/nick.antipa/Documents/Light field/Data/dice/dice-',...
        %         num2str(n,'%02d'),'.png']);
            %a = imread(['/Users/nick.antipa/Documents/Light field/Data/xyzrgb_dragon/xyzrgb_dragon-',...
                %num2str(n,'%02d'),'.png']);
             %im_in = imread(['/Users/nick.antipa/Documents/Light field/Data/DragonAndBunnies/DragonsAndBunnies_5x5_ap6.6/dragons-',...
                 %num2str(count,'%02d'),'.png']);
             a(:,:,q,p) = im_in(rows,cols,monochrome);
            %lf(n,:,:) = a(row,cols,:);
            imagesc(a(:,:,q,p))
            pause(1/24)
            %a_test(:,:,q,p) = reshape(count*nx*ny:count*nx*ny+nx*ny-1,[ny,nx]);
        else
            
             im_in = imread(['/Users/nick.antipa/Documents/Light field/Data/DragonAndBunnies/DragonsAndBunnies_5x5_ap6.6/dragons-',...
                num2str(count,'%02d'),'.png']);
            r(:,:,q,p) = im_in(rows,cols,1);
            g(:,:,q,p) = im_in(rows,cols,2);
            b(:,:,q,p) = im_in(rows,cols,3);                
            imagesc(uint8(cat(3,r(:,:,q,p),g(:,:,q,p),b(:,:,q,p))))
            pause(1/24)
        end
    end
end

%%

str = input('Output mat file name: ','s')
filename = ['./Output/',str,'.mat'] 
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
    filename = ['./Output/',str,'.gif'] 
end

%%
count = 0;
set(0,'CurrentFigure',h8)
for n = 1:nphi
    for m = 1:ntheta
        count = count+1;
        imagesc(uint8(lf_im(n:5:end,m:5:end,:)))
        
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




