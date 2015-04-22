A_in = load('./Output/A_sub_50_50_5_5_5e4_1_deg_CORRECT.mat');
save_file = 0;
nsx = 500;
nsy = 500;
nx = 50;
ny = 50;
ntheta = 5;
nphi = 5;
monochrome = 1;
A_sub1 = A_in.A_sub;
clear A_in
A_sub1 = A_sub1./mean(sum(A_sub1,1));   %Normalize A matrix
%%
lambda = .05;    %regularization
noise = .0005;  %Sensor noise2

h5 = figure(5);
clf
h6 = figure(6);
clf
h7 = figure(7);
clf

lf_final = zeros(ny*nphi,nx*ntheta,1);
in = load('./Output/dragon_bunny_50x50.mat');
if monochrome == 1
    lf = in.lfr;
elseif monochrome == 2
    lf = in.lfg;
else 
    lf = in.lgb;   
end
clear in
sensor = A_sub1 * lf(:);
set(0,'CurrentFigure',h5)
sensor_reshaped = reshape(sensor,[nsx,nsy]);
imagesc(sensor_reshaped);
intens_noisy = sensor + abs(noise*max(sensor)*randn(size(sensor)));

%%
tic
cvx_begin
    variable lf_cvx(numel(lf))
    minimize(norm(A_sub1*lf_cvx-intens_noisy,2)+lambda*norm(lf_cvx,2))
%     subject to
%     lf_cvx >= 0
cvx_end
toc

%%
lf_cvx_reshaped = reshape(lf_cvx,[nphi,ntheta,nx,ny]);
for n = 1:nx
    for m = 1:ny
        if monochrome
            lf_final((n-1)*ntheta+1:n*ntheta,(m-1)*nphi+1:m*nphi) = (lf_cvx_reshaped(:,:,m,n));
        else
            for ncol = 1:3
                lf_final((n-1)*ntheta+1:n*ntheta,(m-1)*nphi+1:m*nphi,ncol) = uint8(recovered_reshaped{ncol}(:,:,m,n));            
            end
        end
    end
end
set(0,'CurrentFigure',h6)
imagesc(uint8(lf_final))
axis image

title('full recovered light field')

count = 0;

set(0,'CurrentFigure',h7)
if save_file
    str = input('Output file name','s')
    filename = ['./Output/',str,'.gif'] 
end
for n = 1:nphi
    for m = 1:ntheta
        count = count+1;
        imagesc(uint8(lf_final(n:5:end,m:5:end,:)))
        if save_file
            if ~monochrome
                [imind,cm] = rgb2ind(uint8(lf_final(n:5:end,m:5:end,:)),256);
                if count == 1;
                  imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
                else
                  imwrite(imind,cm,filename,'gif','WriteMode','append');
                end
            else
                if count == 1;
                  imwrite(uint8(lf_final(n:5:end,m:5:end)),filename,'gif', 'Loopcount',inf);
                else
                  imwrite(uint8(lf_final(n:5:end,m:5:end)),filename,'gif','WriteMode','append');
                end
            end
        end
        axis image
        pause(1/7)
    end
end