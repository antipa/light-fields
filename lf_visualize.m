
%%
monochrome = 1;
if monochrome
    %recin = load('/Users/nick.antipa/Documents/Light_field_data/spheresAndBlock/15_deg_recon_green_medreg_no_noise_10x_1000iter.mat');
    recin = load('/Users/nick.antipa/Documents/MATLAB/Output/ruler_leaf_proxsterov_nonnegonly.mat');
else
    recin_red = load('/Volumes/nick.antipa/MATLAB/Output/7_deg_recon_red_medreg_no_noise_10x.mat');
    recin_green = load('/Volumes/nick.antipa/MATLAB/Output/7_deg_recon_green_medreg_no_noise_10x.mat');
    recin_blue = load('/Volumes/nick.antipa/MATLAB/Output/7_deg_recon_blue_medreg_no_noise_10x.mat');
    recin.recovered = cat(5,recin_red.recovered,recin_green.recovered,recin_blue.recovered);
end

recin.settings = settings;
recin.recovered = recovered_2dtv;
%recin.recovered = flip(lf_in.lfg,3);
nangles = recin.settings.P;
d = linspace(-1,1,nangles);
[dX,dY] = meshgrid(d,d);
R = dX.^2+dY.^2;
aper = R<2;
recap = zeros(size(recin.recovered));

lfsettings.D = '4d';
lfsettings.nx = recin.settings.N;
lfsettings.ny = recin.settings.M;
lfsettings.ntheta = recin.settings.P;
lfsettings.nphi = recin.settings.Q;
lfsettings.snake = 1;

switch lower(lfsettings.D)
    case('3d')
        for mm = 1:size(recin.recovered,3)
            [rb,cb] = ind2sub([nangles,nangles],mm);


            imagesc(transpose(recin.recovered(:,:,mm)));

            caxis([0 prctile(recin.recovered(:),99)])
            colormap gray
                axis image
            drawnow

            pause(1/100)
        end
    case('4d')
        if monochrome
            for m = 1:size(recin.recovered,1)
                for n = 1:size(recin.recovered,2)
                    imagesc(transpose(squeeze(recin.recovered(n,m,:,:))))
                    caxis([0 prctile(recin.recovered(:),99)])
                    colormap gray
                    axis image
                    drawnow
                    pause(1/24)
                end
            end
        else
            for m = 1:size(recin.recovered,1)
                for n = 1:size(recin.recovered,2)
                    subap = cat(3,transpose(squeeze(recin.recovered(n,m,:,:,1))),...
                        transpose(squeeze(recin.recovered(n,m,:,:,2))),...
                        transpose(squeeze(recin.recovered(n,m,:,:,3))));
                    imagesc(subap/2^16);
                    caxis([0 prctile(recin.recovered(:),99)])
                    colormap gray
                    axis image
                    drawnow
                    pause(1/24)
                end
            end
        end
end
%%

for s = 1
   switch lower(lfsettings.D)
       case('3d')
            refocused = shift_and_add(permute(recin.recovered,[2,1,3]),lfsettings,s);
       case('4d')
           if monochrome
               refocused = shift_and_add(permute(recin.recovered,[1,2,4,3]),lfsettings,s);
           else
               refocused_r = shift_and_add(permute(recin.recovered(:,:,:,:,1),[1,2,4,3]),lfsettings,s);
               refocused_g = shift_and_add(permute(recin.recovered(:,:,:,:,2),[1,2,4,3]),lfsettings,s);
               refocused_b = shift_and_add(permute(recin.recovered(:,:,:,:,3),[1,2,4,3]),lfsettings,s);
               refocused = cat(3,refocused_r,refocused_g,refocused_b);
           end
   end

   imagesc(refocused/max(refocused(:)))
      %caxis([0 prctile(refocused(:),99.5)])
      caxis([0 1])
   axis image
   colormap gray
   drawnow
end

