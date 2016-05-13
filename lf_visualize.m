%%
recin = load('/Users/nick.antipa/Documents/Diffusers/1deg_yes4f_flea3_20160510/king_queen_3dtv_9x155_5p2_ap.mat');
nangles = recin.settings.P;
d = linspace(-1,1,nangles);
[dX,dY] = meshgrid(d,d);
R = dX.^2+dY.^2;
aper = R<3;
recap = zeros(size(recin.recovered));
for mm = 1:size(recin.recovered,3)
    [rb,cb] = ind2sub([nangles,nangles],mm);
    
    recap(:,:,mm) = (transpose(recin.recovered(:,:,mm)));
    imagesc(recap(:,:,mm));

    caxis([0 prctile(recin.recovered(:),99)])
    colormap gray
        axis image
    drawnow

    pause(1/100)
end
%%
lfsettings.D = '3d';
lfsettings.nx = recin.settings.N;
lfsettings.ny = recin.settings.M;
lfsettings.ntheta = nangles;
lfsettings.nphi = nangles;

for s = -2.2:.2:2.4
   refocused = shift_and_add(permute(recap,[1,2,3]),lfsettings,s);
   imagesc(refocused)
   axis image
   drawnow
end

