nsx = 1000;
nsy = 1000;
nx = 50;
ny = 50;
ntheta = 5;
nphi = 5;
lfr = load('./recovered_50_50_5_5_lambda_00005_noise_005.mat');
lfg = load('./recovered_50_50_5_5_lambda_00005_noise_005_green.mat');
lfb = load('./recovered_50_50_5_5_lambda_00005_noise_005_blue.mat');
lfr_reshaped = reshape(lfr.recovered,[nphi,ntheta,nx,ny]);
lfg_reshaped = reshape(lfg.recovered,[nphi,ntheta,nx,ny]);
lfb_reshaped = reshape(lfb.recovered,[nphi,ntheta,nx,ny]);
sensor_r = load('./dragon_bunny_sensor_50_50_5_5_red.mat');
sensor_g = load('./dragon_bunny_sensor_50_50_5_5_green.mat');
sensor_b = load('./dragon_bunny_sensor_50_50_5_5_blue.mat');
sensor_rgb = (zeros(size(sensor_r.sensor_reshaped,1),size(sensor_r.sensor_reshaped,2),3));
sensor_rgb(:,:,1) = (sensor_r.sensor_reshaped);
sensor_rgb(:,:,2) = (sensor_g.sensor_reshaped);
sensor_rgb(:,:,3) = (sensor_b.sensor_reshaped);
sensor_rgb = uint8(sensor_rgb*255/max(max(max(sensor_rgb))));
lf_final = uint8(zeros(nx*ntheta,ny*nphi,3));
for n = 1:nx
    for m = 1:ny
        for l = 1:3
            if l == 1
                lf_final((n-1)*ntheta+1:n*ntheta,(m-1)*nphi+1:m*nphi,l) = uint8(lfr_reshaped(:,:,m,n));
            elseif l == 2
                lf_final((n-1)*ntheta+1:n*ntheta,(m-1)*nphi+1:m*nphi,l) = uint8(lfg_reshaped(:,:,m,n));
            elseif l ==3
                lf_final((n-1)*ntheta+1:n*ntheta,(m-1)*nphi+1:m*nphi,l) = uint8(lfb_reshaped(:,:,m,n));
            end
        end
    end
end
imagesc(lf_final)

%%
for n = 1:nphi
    for m = 1:ntheta
        imagesc(lf_final(n:5:end,m:5:end,1:3))        
        
    end
end


