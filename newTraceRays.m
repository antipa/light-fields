%also added index variable
%gets data from tracerays_settings.m
use_defaults = 1;
if use_defaults
    try
        tracerays_settings;
    catch exception
       display('Error loading settings. Read instructions in tracerays_settings.m.template.'); 
    end
else
    %You need to load settings externally!
end

%get diffuser surface

%Enable code profiler for timing etc. Doesn't work with parfor
if settings.code_profiler
    profile on
end

%If diffuser has not been loaded externally, load this file
if ~settings.preload_diff
    in = load('/Users/nick.antipa/Documents/MATLAB/Output/1deg_flea_4f_gptie_dz381um_surface_1552.mat');
    diffuser_in = in.filtered * settings.strength;   %Scale height
    x = in.x;   %get physical units grid from file
end

%Upsample diffuser if needed
if settings.diff_upsample
    diffuser = imresize(diffuser_in,settings.diff_upsample,'bicubic');
    x = linspace(min(x),max(x),numel(x)*settings.diff_upsample);
else
    diffuser = diffuser_in;
end


y = x;
px = mean(diff(x)); %diffuser "pixel" size in um/pixel
clear in
range_idx = floor(settings.x_range/px);     %Compute the range over which the light field illuminates the diffuser, in index units
dx_idx = floor(settings.x_range/settings.N/px);    %Lateral extent of each bundle in index units

%Don't allow index less than 1
if dx_idx<1
    dx_idx = 1;
end

%Repeat same calculations for y direction
y_idx = (y-min(y))/px;   %x vector as index
range_idy = floor(settings.y_range/px);
dy_idx = floor(settings.y_range/settings.M/px);

%If y window size is less than 1, assume we are doing a 1d problem
if dy_idx<1
    dy_idx = 0;
end

% Calculate the sensor size in physical units
ssize = [max(settings.y_range,1), max(settings.x_range,1)];   %Sensor size in microns
dpx = ssize(2)/settings.npx;    %pixel size on sensor in micron/pixel
dpy = ssize(1)/settings.npy;
xs = min(x):dpx:max(x);   %Create x and y vectors in physical units for sensor
ys = xs;

%Use y_range=0 to tell the code to run a 1d problem
if settings.y_range == 0
    settings.npy = 1;   %Setup a single pixel in 1 for the 1d problem
    settings.ph_range = 0;   %Zero out the y-direction angular range to all rays stay within 1d sensor
end

%Repeat same for x if 1d in x direction
if settings.x_range == 0
    settings.npx = 1;
    settings.th_range = 0;
end

%Angular ranges -- ray angles are distributed uniformly across this range
dph = settings.ph_range/settings.P;
dth = settings.th_range/settings.Q;

%If not using parfor, preallocate the A matrix
% Po
%if ~settings.useParfor
    %creates sparse form of matrix
%    A_sub = sparse(settings.npx*settings.npy,settings.N*settings.M*settings.P*settings.Q);
%end

%preallocate cell to hold values and indices for sparse matrix creation
%A_vals = cell(settings.M*settings.N*settings.P*settings.Q,1);    %DELETE THIS?
r_outc = cell(settings.M*settings.N*settings.P*settings.Q,1);    %Row index goes here
c_outc = cell(settings.M*settings.N*settings.P*settings.Q,1);    %Columns
v_outc = cell(settings.M*settings.N*settings.P*settings.Q,1);    %values

%Setup gradient of surface for refraction calculations
if dy_idx == 0
    Fx = gradient(diffuser);
    Fy = zeros(size(Fx));  %No Fy if 1d 
elseif dx_idx == 0
    Fy = gradient(diffuser);
    Fx = zeros(size(Fy));
else 
    [Fx, Fy] = gradient(diffuser);   %2d gradient
end

%Convert gradient to correct units
Fx = Fx/px;
Fy = Fy/px;

%Setup figure windows and handles for visualization
if settings.vis
    h1 = figure(1);
    clf
    h2 = figure(2);
    clf
end

if settings.vis_sensor
    h4 = figure(4);
    clf
end

if settings.vis_prop
    h3 = figure(3);
    clf
    if settings.vis_hist
        h6 = figure(6);
        clf
    end
    if settings.vis_xcorr
        h7 = figure(7);
        clf
    end
end

%Setup a waitbar for non parfor
if ~settings.useParfor
    if settings.waitbar
        h5 = waitbar(0,'beginning');
    end
end

%Loop over grid positions on surface
[Xg, Yg] = meshgrid(-dx_idx*px/2:px:dx_idx*px/2,-dx_idx*px/2:px:dx_idx*px/2);
tstart1 = tic;
if settings.useParfor
    %Matlab makes you use one index in parfor loops, so linearly index each
    %variable
    parfor LF_index = 1:settings.M*settings.N*settings.P*settings.Q 
        ytstart = tic;
        qq = 1+mod(LF_index-1,settings.Q);    %Calculate the ind
        pp = 1+mod(floor((LF_index-1)/settings.Q),settings.P);
        nn = 1+mod(floor((LF_index-1)/settings.P/settings.Q),settings.N);
        mm = 1+mod(floor((LF_index-1)/settings.P/settings.Q/settings.N),settings.M);
        %Get indices of region in y direction
        midx = (((mm-1)*dy_idx+1):(mm*dy_idx+1))';
        midx_1 = midx(1);
        %Generate y vector in physical units in local coordinates
        %(i.e. origin shifts wich each new region
        yg = px*(midx-midx(1));
        
        nidx = (((nn-1)*dx_idx+1):(nn*dx_idx+1))';
        nidx_1 = nidx(1);
        xg = (nidx-nidx(1))*px;
        
        Fx_crop = Fx(midx,nidx);
        Fy_crop = Fy(midx,nidx);
        
        ph = dph*(rand(settings.nrays,1)-pp+settings.P/2);
        th = dth*(rand(settings.nrays,1)-qq+settings.Q/2); %random theta points
        
        
        %Generate random (x,y) positions within dx_idx*px wide
        %square -- these are the physical locations where each ray will
        %strike the diffuser
        yr = (dy_idx)*px*rand(settings.nrays,1);
        xr = (dx_idx)*px*rand(settings.nrays,1);
           
        %Estimate slope of diffuser at each (xr,yr) position
        if dy_idx==0
            Fxr = interp1(xg,Fx_crop,xr);  %Interpolate x gradient
            Fyr = zeros(size(Fxr));
        elseif dx_idx==0
            Fyr = interp1(yg,Fy_crop,yr);  %Interpolate x gradient
            Fxr = zeros(size(Fyr));
        else
            %Choose nearest neightbor interpolation for speed, or default
            %for accuracy
            if strcmpi(settings.interp_method,'nearest')
                Fxr = interp2(xg,yg,Fx_crop,xr,yr,settings.interp_method);  %Interpolate x gradient
                Fyr = interp2(xg,yg,Fy_crop,xr,yr,settings.interp_method);  %Interpolate y gradiet
            else
                Fxr = interp2(xg,yg,Fx_crop,xr,yr);  %Interpolate x gradient
                Fyr = interp2(xg,yg,Fy_crop,xr,yr);  %Interpolate y gradiet
            end
        end
        
        %Choose paraxial approximation or real vector based tracing
        if ~settings.paraxial
            %Refraction starts here ---------------------
            [uxp, uyp, uzp] = refraction(Fxr, Fyr, th, ph, settings.indexDiff, settings.indexEnv,'angles');
            %End refraction-------------

            %propagate to output plane by a distance z
            [yo, xo] = propagation(uyp, uzp, settings.z0, yr, uxp, xr);
            %                 yo = uyp./uzp*settings.z0+yr;
            %                 xo = uxp./uzp*settings.z0+xr;
        elseif settings.paraxial
            
            xo = xr+settings.z0*(settings.dn*Fxr+th*pi/180*settings.indexDiff);
            yo = yr+settings.z0*(settings.dn*Fyr+ph*pi/180*settings.indexDiff);
        end
        
        gatherer = gather_rays_nohist(xo,yo,settings.npx,settings.npy,dpx,dpy,nidx_1,midx_1,px);
        %Only gather rays if we've got more than 1 ray. It sounds
        %dumb, but the hist code doesn't work when there's 1 ray,
        %and I don't care about a single ray anyway.
        if nnz(gatherer)>1
            %Gather rays on sensor. This is essentially a 2D histogram, but
            %outputs appropriate row, column and value cells for passing
            %into sparse matrix construction at end of loop
            [r_outc{LF_index}, c_outc{LF_index}, v_outc{LF_index}] = ...
                build_A_matrix_sparse(gatherer,LF_index);
        end
    end
else
    for mm = 1:max(1,settings.M)
        ytstart = tic;
        %Get indices of region in y direction
        midx = (((mm-1)*dy_idx+1):(mm*dy_idx+1))';
        %Generate y vector in physical units in local coordinates
        %(i.e. origin shifts wich each new region
        yg = px*(midx-midx(1));
        %yr = dx_idx*px*(rand(settings.nrays,1));
        for nn = 1:max(1,settings.N)
            nidx = (((nn-1)*dx_idx+1):(nn*dx_idx+1))';   %x direction (column) index for pulling diffuser patch out of master diffuser surface array
            xg = (nidx-nidx(1))*px;
            diff_crop = diffuser(midx,nidx);
            
            [NIDX, MIDX] = meshgrid(nidx,midx);  %The meshgridded version of nidx and midx are only used for plotting
            
            Fx_crop = Fx(midx,nidx);
            Fy_crop = Fy(midx,nidx);
            for pp = 1:settings.P
                ph = dph*(rand(settings.nrays,1)-pp+settings.P/2);    %Random phi points
                for qq = 1:settings.Q
                    tstart = tic;
                    th = dth*(rand(settings.nrays,1)-qq+settings.Q/2); %random theta points
                    
                    %Generate random (x,y) positions within dx_idx*px wide
                    %square
                    yr = (dy_idx)*px*rand(settings.nrays,1);
                    xr = (dx_idx)*px*rand(settings.nrays,1);
                    
                    
                    %Calculate surface norm at random positions by
                    %interpolation within gradient
                    
                    %Get z and surface normal at each random (xr(i),yr(i)) pair
                    if dy_idx==0
                        Fxr = interp1(xg,Fx_crop,xr);  %Interpolate x gradient
                        Fyr = zeros(size(Fxr));
                        zr = interp1(xg,diff_crop,xr);   %Interpolate surface
                    elseif dx_idx==0
                        Fyr = interp1(yg,Fy_crop,yr);  %Interpolate x gradient
                        Fxr = zeros(size(Fyr));
                        %zr = interp1(yg,diff_crop,yr);   %Interpolate surface
                    else
                        %zr = interp2(xg,yg,diff_crop,xr,yr);   %Interpolate surface
                        Fxr = interp2(xg,yg,Fx_crop,xr,yr);  %Interpolate x gradient
                        Fyr = interp2(xg,yg,Fy_crop,xr,yr);  %Interpolate y gradiet
                    end
                    
                    %Choose between paraxial and real snell's law
                    if ~settings.paraxial
                        %Refraction starts here ---------------------
                        [uxp, uyp, uzp] = refraction(Fxr, Fyr, th, ph, settings.indexDiff, settings.indexEnv,'angles');
                        %End refraction-------------

                        %propagate to output plane by a distance z
                        [yo, xo] = propagation(uyp, uzp, settings.z0, yr, uxp, xr);
                        %                 yo = uyp./uzp*settings.z0+yr;
                        %    
                        
                        xo = uxp./uzp*settings.z0+xr;
                    elseif settings.paraxial
                        %yoxo = [yr,xr]+settings.z0*[settings.dn*Fyr+ph*pi/180*settings.indexDiff,settings.dn*Fxr+th*pi/180*settings.indexDiff];
                        xo = xr+settings.z0*(settings.dn*Fxr+th*pi/180*settings.indexDiff);
                        yo = yr+settings.z0*(settings.dn*Fyr+ph*pi/180*settings.indexDiff);
                        %x_retina = x_pupil+z*(dW(x_pupil)/dx+theta_pupil)
                        %n1*sin(phi1)=n2*sin(phi2) ~=
                        %n1*phi1=n2*phi2
                        %W = (n_lens - n_air)*Z(x,y)
               
                        
                    end
                    
                    gatherer = gather_rays_nohist(xo,yo,settings.npx,settings.npy,dpx,dpy,nidx(1),midx(1),px);
                    %Gather rays on sensor
                    if settings.npy == 1
                        gatherer = gather_rays(xo,yo,settings.npx,settings.npy,dpx,dpy,nidx(1),midx(1),px);
                    end
                    
                    
                    %Only gather rays if we've got more than 1 ray. It sounds
                    %dumb, but the hist code doesn't work when there's 1  ray,
                    %and I don't care about a single ray anyway.
                    if nnz(gatherer)>1
                        
                        %make light field index. This is so it matches with
                        %parfor version
                        LF_index = settings.P*settings.Q*settings.N*(mm-1)+settings.P*settings.Q*(nn-1)+settings.Q*(pp-1)+qq;

                        [r_outc{LF_index}, c_outc{LF_index}, v_outc{LF_index}] = ...
                            build_A_matrix_sparse(gatherer,LF_index);
                        

                        if settings.vis_prop %animation to visualize propagation after refraction
                        	
                            zvec = 0:settings.zstep:settings.zmax;
                            v_mat = zeros(length(zvec),settings.npx);
                            for z = zvec;
                                set(0,'CurrentFigure',h3)
                                 yo = uyp*z./uzp+yr;
                                xo = uxp*z./uzp+xr;
                                if settings.npy~=1
                                    gatherer1 = hist4(xo,yo,settings.npx,settings.npy,dpx,dpy,nidx(1),midx(1),px);
                                elseif settings.npy == 1
                                    gatherer1 = gather_rays(xo,yo,settings.npx,settings.npy,dpx,dpy,nidx(1),midx(1),px);
                                end
                                if nnz(gatherer1)>1
                                    x_sensor = [0:settings.npx-1]*dpx;
                                    y_sensor = [0:settings.npy-1]*dpy;
                                    thcorr = tand(settings.th_range)*settings.zmax;
                                    apcorrx = max(max(abs(Fx_crop)))*settings.zmax;
                                    phcorr = tand(settings.ph_range)*settings.zmax;
                                    apcorry = max(max(abs(Fy_crop)))*settings.zmax;
                                    xmin = max(0,xg(1)+nidx(1)*px-thcorr-apcorrx);
                                    xmax = min(x_sensor(end),xg(end)+nidx(1)*px+thcorr+apcorrx);
                                    ymin = max(0,yg(1)+midx(1)*px-phcorr-apcorry);
                                    ymax = min(y_sensor(end),yg(end)+midx(1)*px+phcorr+apcorry);
                                    
                                    if dy_idx>1 && dx_idx>1
                                        
                                        imagesc(gatherer1,'XData',x_sensor,'YData',y_sensor)
                                        hold on
                                        axis image
                                        axis([xmin xmax ymin ymax])
                                        if settings.vis_hist
                                            set(0,'CurrentFigure',h6)
                                            clf
                                            hist(full(gatherer(gatherer1>0)))
                                            title(['histogram of intensity z=',num2str(z)])
                                            drawnow
                                        end
                                    elseif dy_idx==0
                                        stairs(x_sensor,gatherer1)
                                        xlim([xmin xmax])
                                        xlabel('\mum')
                                        v_mat(zvec==z,:) = gatherer1;
                                    elseif dx_idx==0
                                        stairs(y_sensor,gatherer1)
                                        xlim([ymin ymax])
                                        xlabel('\mum')
                                    end
                                    set(0,'CurrentFigure',h3)
                                    title(['intensity at z=',num2str(z),' for ',num2str(settings.nrays),...
                                        ' rays, \theta=',num2str(dth*(-qq+settings.Q/2+.5)),...
                                        ' \phi=',num2str(dph*(-pp+settings.P/2+.5))])
                                    hold off
                                    
                                    if settings.vis_xcorr
                                        [rxcorr,cxcorr] = find(gatherer1);
                                        rv = min(rxcorr):max(rxcorr);
                                        cv = min(cxcorr):max(cxcorr);
                                        %xcorred = xcorr2(full(gatherer1(rv,cv)));
                                        xcorred = (ifftshift(ifft2(abs(fft2(full(gatherer1)-mean2(gatherer1))).^2)));
                                        set(0,'CurrentFigure',h7)
                                        if find(z==zvec(1))
                                            zleg = {num2str(z)};
                                            clf
                                        else
                                            zleg{end+1} = num2str(z);
                                        end
                                        subplot(2,1,1)
                                        imagesc(xcorred)
                                        title(['autocorrelation at z=',num2str(z)])
                                        subplot(2,1,2)
                                        %hold on
                                        lout = xcorred(middle(xcorred),:);
                                        plot(lout(middle(lout):middle(lout)+round(dx_idx)))
                                        title('Lineout')
                                        %legend(zleg);
                                        drawnow
                                    end
                                    
                                    if ~isempty(settings.save_prop)
                                        filename = settings.save_prop;
                                        frame = getframe(h3);
                                        im = frame2im(frame);
                                        [imind,cm] = rgb2ind(im,256);
                                        if z==zvec(1);
                                            imwrite(imind,cm,filename,'gif', 'Loopcount',inf);
                                        else
                                            imwrite(imind,cm,filename,'gif','WriteMode','append');
                                        end
                                    end
                                    
                                    pause(1/100)
                                end
                            end
                        end  %End propagation visualization
                        
                        if settings.vis_sensor
                            set(0,'CurrentFigure',h4)
                            x_sensor = [0:settings.npx-1]*dpx;
                            y_sensor = [0:settings.npy-1]*dpy;
                            thcorr = tand(settings.th_range)*settings.z0;
                            apcorrx = max(max(abs(Fx_crop)))*settings.z0;
                            phcorr = tand(settings.ph_range)*settings.z0;
                            apcorry = max(max(abs(Fy_crop)))*settings.z0;
                            xmin = max(0,xg(1)+nidx(1)*px-thcorr-apcorrx);
                            xmax = min(x_sensor(end),xg(end)+nidx(1)*px+thcorr+apcorrx);
                            ymin = max(0,yg(1)+midx(1)*px-phcorr-apcorry);
                            ymax = min(y_sensor(end),yg(end)+midx(1)*px+phcorr+apcorry);
                            if dy_idx>1 && dx_idx>1
                                imagesc(gatherer,'XData',x_sensor,'YData',y_sensor)
                                hold on
                                xlabel('\mum')
                                ylabel('\mum')
                                axis image
                                axis([xmin xmax ymin ymax])
                            elseif dy_idx==0
                                stairs(x_sensor,gatherer)
                                xlim([xmin xmax])
                                xlabel('\mum')
                            elseif dx_idx==0
                                stairs(y_sensor,gatherer)
                                xlim([ymin ymax])
                                xlabel('\mum')
                            end
                            title(['intensity at z=',num2str(settings.z0),' for ',num2str(settings.nrays),...
                                ' rays, \theta=',num2str(dth*(-qq+settings.Q/2+.5)),...
                                ' \phi=',num2str(dph*(-pp+settings.P/2+.5))])
                            hold off
                            pause(1/24)
                        end
                    end
                    
                    %               %% Visualization
                    if settings.vis
                        %Calculate global x and y index locations for rays
                        scl = 1;
                        nnr = xr/px+nidx(1);
                        mmr = yr/px+midx(1);
                        set(0,'CurrentFigure',h1)
                        if dy_idx == 0
                            
                            clf
                            plot(xg+px*nidx(1),diff_crop)
                            hold on
                            quiver(xr+px*nidx(1),zr,-Fxr,ones(size(Fxr)),5,'ShowArrowHead','off','LineStyle','-','Color',[1 .95 .8])
                            quiver(xr+px*nidx(1),zr,Fxr,-ones(size(Fxr)),5,'ShowArrowHead','off','LineStyle','-','Color',[1 .95 .8])
                            quiver(xr+px*nidx(1),zr,uxp,uzp,10,'Color','b')
                            quiver(xr+px*nidx(1),zr,uxn,-uzn,10,'ShowArrowHead','off','Color','b')
                            axis equal
                            
                        elseif dx_idx == 0
                        else
                            clf
                            surf(NIDX,MIDX,diff_crop/scl,'linestyle','none')
                            hold on
                            quiver3(nnr,mmr,zr/scl,-Fxr/scl,-Fyr/scl,ones(size(Fxr)))
                            quiver3(nnr,mmr,zr/scl,-uxn,-uyn,-uzn,100)
                            quiver3(nnr,mmr,zr/scl,uxp,uyp,uzp,100)
                            scatter3(nnr,mmr,zr/scl)
                            axis equal
                            view([0,0])
                            xlim([-2000 2000])
                            ylim([-2000 2000])
                            zlim([0 5000])
                            hold off
                            
                            set(0,'CurrentFigure',h2)
                            surf(NIDX,MIDX,diff_crop/scl,'linestyle','none')
                            hold on
                            quiver3(nnr,mmr,zr,uxn,uyn,uzn,1/scl)
                            axis equal
                            grid on
                            view([11,58])
                        end
                    end
                end
            end
        end
        ytend = toc(ytstart);
        tleft = (settings.M*settings.N*settings.P*settings.Q-LF_index)/(settings.N*settings.P*settings.Q)*ytend;
        tmin = floor(tleft/60);
        tsec = mod(tleft/60,tmin)*60;
        if settings.waitbar
            waitbar(LF_index/settings.M/settings.N/settings.P/settings.Q,h5,...
                [num2str(100*(LF_index)/(settings.M*settings.N*settings.P*settings.Q)),'% done. ',num2str(settings.M*settings.N*settings.P*settings.Q),...
                ' total, ',num2str(ytend),' seconds per pass. ~',...
                num2str(tmin),':',num2str(tsec,'%.0f'),' to go.'])
        end
    end
end

r_out = vertcat(r_outc{:});
c_out = vertcat(c_outc{:});
v_out = vertcat(v_outc{:});
A_sub = sparse(r_out,c_out,v_out,settings.npx*settings.npy,settings.N*settings.M*settings.P*settings.Q);
total_runtime = toc(tstart1);
fprintf([num2str(total_runtime),' seconds for ',num2str(settings.P*settings.Q*settings.N*settings.M*settings.nrays),' rays. ',...
    num2str(settings.P*settings.Q*settings.N*settings.M*settings.nrays/total_runtime/1e6),' million rays/sec\n'])

if settings.code_profiler
    profile viewer
end