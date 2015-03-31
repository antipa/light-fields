ntheta = 5;
nphi = 5;
nx = 50;
ny = 50;
xstart = 600;
ystart = 205;

cols = xstart:xstart+nx-1;
rows = ystart:ystart+ny-1;
clf
lf = zeros(ntheta,npix,3);
figure(5),clf
for n = 1:ntheta
%     a = imread(['/Users/nick.antipa/Documents/Light field/Data/dice/dice-',...
%         num2str(n,'%02d'),'.png']);
    %a = imread(['/Users/nick.antipa/Documents/Light field/Data/xyzrgb_dragon/xyzrgb_dragon-',...
        %num2str(n,'%02d'),'.png']);
     a = imread(['/Users/nick.antipa/Documents/Light field/Data/DragonAndBunnies/DragonsAndBunnies_5x5_ap6.6/dragons-',...
         num2str(n,'%02d'),'.png']);
    lf(n,:,:) = a(row,cols,:);
    imagesc(a)
    pause(1)
end

