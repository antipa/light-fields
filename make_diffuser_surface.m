
%function surface_out = make_diffuser_surface(X,Y,dtheta_max,freq)
%Generate diffuser surface with characteristic angle spread and spatial
%frequency over 2d grid defined by vectors x and y
write_file = 0
save_m = 1
[X, Y] = meshgrid(-2540:1:2540,-2540:1:2540);   %in microns
period = 20;                %microns
nrows = size(X,1);
ncols = size(X,2);
fname = '../Output/half_deg.mat';
if size(X)~=size(Y)
    error('X and Y must be same size')
end
delta_n = .5;
dtheta_max = 1*pi/180;
dx = mean(diff(X(1,:)));    
period_p = round(period/dx);  %width of gaussian
%amplitude of gaussian dtheta = delta_n*dz/dx
%dz = delta_n*period_p/dtheta_max
dz = dtheta_max*period_p/delta_n;
seeds = randn(nrows,ncols)*dz*.3;
fsize = round(min(period_p*100,size(seeds,1)));
K = fspecial('gaussian',fsize,period_p/2);
K = K*sum(sum(K));

%filtered = imfilter(seeds,K);
filtered = ifft2(fft2(K,size(seeds,1),size(seeds,2)).*(fft2(seeds)));
%[Fx Fy] = gradient(filtered)*;
Fx = zeros(size(filtered));
Fy = Fx;
Fxy = Fx;
if write_file
    %[Fxx Fxy] =gradient(Fx);
    fid = fopen('diffuser_high.DAT','w')
    fprintf(fid,'%i\t%i\t%.3f\t%.3f\t%i\t%.4f\t%.4f\n',ncols,nrows,dx/1000,dx/1000,0,0,0)

    for m = 1:numel(Fx)
       fprintf(fid,'%.8f\t%.8f\t%.8f\t%.8f\n',filtered(m),Fx(m),Fy(m),Fxy(m));
    end
    fclose(fid)
end
if save_m
    x = X(1,:);
    save(fname,'x','filtered')
end
%%
%Trace ray
% for n = 1:numel(Y(:,1))
%     for m = 1:numel(X(1,:))
%         
%     x_sensor = -1000:5:1000;
%     y_sensor = x_sensor;
%     [X_sensor, Y_sensor] = meshgrid(x_sensor,y_sensor);
%     dz = 200;  %microns
%     %m = 200;
%     %n = 200;
%     x = X_sensor(m,n)
%     y = Y_sensor(m,n)
%     NA = .25;
%     delta_theta_max = max(max(sqrt(Fx.^2+Fy.^2)*dx));
%     th_max = delta_theta_max + asin(NA);
%     theta = linspace(-th_max,th_max,35);
%     phi = theta;
%     [TH PH] = meshgrid(theta,phi);
%     for mm = 1:numel(theta);
%         for nn = 1:numel(phi);
%             delta_y = dz*phi(nn);
%             delta_x = dz*theta(mm);
%             x_diff(mm,nn) = f_of(X(1,:),X(1,:),x+delta_x);
%             y_diff(mm,nn) = f_of(Y(:,1),Y(:,1),y+delta_y);
%             [r c] = find(X==x_diff(mm,nn)&Y==y_diff(mm,nn));
%             dthx = Fx(r,c)*dx*10;
%             dphy = Fy(r,c)*dx*10;
%             thx(mm,nn) = theta(mm)+dthx;
%             phy(mm,nn) = phi(nn)+dphy;
%         end
%     end
%     figure(1),hold on
%     scatter(x_diff(:),thx(:))
%     end
% end
