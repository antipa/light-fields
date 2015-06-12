%also added index variable
%gets data from tracerays_settings.m
try
    tracerays_settings;
catch exception
   display('Error loading settings. Read instructions in tracerays_settings.m.template.'); 
end

%get diffuser surface

%Notes:
%Make sensor size dynamic

profile on
in = load('../Output/diffuser.mat');
diffuser_in = in.filtered * strength;
x = in.x;
if diff_upsample
    diffuser = imresize(diffuser_in,diff_upsample,'bicubic');
    x = linspace(min(x),max(x),numel(x)*diff_upsample);
else
    diffuser = diffuser_in;
end
y = x;
px = mean(diff(x)); %diffuser "pixel" size in um/pixel

range_idx = floor(x_range/px);
dx_idx = floor(x_range/N/px);
if dx_idx<1
    dx_idx = 1;
end

y_idx = (y-min(y))/px;   %x vector as index
range_idy = floor(y_range/px);
dy_idx = floor(y_range/M/px);
if dy_idx<1
    dy_idx = 0;
end

ssize = [max(y_range,1), max(x_range,1)];   %Sensor size in microns
dpx = ssize(2)/npx;
dpy = ssize(1)/npy;
xs = min(x):dpx:max(x);
ys = xs;

if y_range == 0
    npy = 1;
    ph_range = 0;
end

if x_range == 0
    npx = 1;
    th_range = 0;
end
dph = ph_range/P;
dth = th_range/Q;
if ~useParfor
    %creates sparse form of matrix
    A_sub = sparse(npx*npy,N*M*P*Q);
end
%preallocate
A_row_index = cell(M*N*P*Q,1);  %Row index fo
A_vals = A_row_index;
r_outc = cell(M*N*P*Q,1);
c_outc = cell(M*N*P*Q,1);
v_outc = cell(M*N*P*Q,1);

%Setup gradients
if dy_idx == 0
    Fx = gradient(diffuser);
    Fy = zeros(size(Fx));
elseif dx_idx == 0
    Fy = gradient(diffuser);
    Fx = zeros(size(Fy));
else
    [Fx, Fy] = gradient(diffuser);
end

%Setup figure windows and handles for visualization
if vis
    h1 = figure(1);
    clf
    h2 = figure(2);
    clf
end
if vis_prop
    h3 = figure(3);
    clf
end
if vis_sensor
    h4 = figure(4);
    clf
end
if ~useParfor
    h5 = waitbar(0,'beginning');
end
%Loop over grid positions on surface
[Xg, Yg] = meshgrid(-dx_idx*px/2:px:dx_idx*px/2,-dx_idx*px/2:px:dx_idx*px/2);
tstart1 = tic;
if useParfor
    parfor LF_index = 1:M*N*P*Q
        ytstart = tic;
        qq = 1+mod(LF_index-1,Q);
        pp = 1+mod(floor(LF_index/Q),P);
        nn = 1+mod(floor(LF_index/P/Q),N);
        mm = 1+mod(floor(LF_index/P/Q/N),M);
        %Get indices of region in y direction
        midx = (((mm-1)*dy_idx+1):(mm*dy_idx+1))';
        midx_1 = midx(1);
        %Generate y vector in physical units in local coordinates
        %(i.e. origin shifts wich each new region
        yg = px*(midx-midx(1));
        %yr = dx_idx*px*(rand(nrays,1));
        
        nidx = (((nn-1)*dx_idx+1):(nn*dx_idx+1))';
        nidx_1 = nidx(1);
        xg = (nidx-nidx(1))*px;
        diff_crop = diffuser(midx,nidx);
        
        %[NIDX, MIDX] = meshgrid(nidx,midx);
        
        Fx_crop = Fx(nidx,midx);
        Fy_crop = Fy(nidx,midx);
        
        %LF_index_start = P*Q*N*(mm-1)+P*Q*(nn-1);
        %LF_index_start = P*Q*N*(mm-1)+P*Q*(nn-1)+Q*(pp-1);
        
        
        
        ph = dph*(rand(nrays,1)-pp+P/2);
        th = dth*(rand(nrays,1)-qq+Q/2); %random theta points
        
        
        %Generate random (x,y) positions within dx_idx*px wide
        %square
        yr = (dy_idx)*px*rand(nrays,1);
        xr = (dx_idx)*px*rand(nrays,1);
        
        
        %Calculate surface norm at random positions by
        %interpolation within gradient
        
        %Get z and surface normal at each random (xr(i),yr(i)) pair
        if dy_idx==0
            Fxr = interp1(xg,Fx_crop,xr);  %Interpolate x gradient
            Fyr = zeros(size(Fxr));
            %zr = interp1(xg,diff_crop,xr);   %Interpolate surface
        elseif dx_idx==0
            Fyr = interp1(yg,Fy_crop,yr);  %Interpolate x gradient
            Fxr = zeros(size(Fyr));
            %zr = interp1(yg,diff_crop,yr);   %Interpolate surface
        else
            %zr = interp2(xg,yg,diff_crop,xr,yr);   %Interpolate surface
            Fyr = interp2(xg,yg,Fx_crop',xr,yr);  %Interpolate x gradient
            Fxr = interp2(xg,yg,Fy_crop',xr,yr);  %Interpolate y gradiet
        end
        
        %Refraction starts here ---------------------
        [uxp, uyp, uzp] = refraction(Fxr, Fyr, th, ph, index);
        %End refraction-------------
        
        %propagate to output plane by a distance z
        [yo, xo] = propagation(uyp, uzp, z0, yr, uxp, xr);
        %                 yo = uyp./uzp*z0+yr;
        %                 xo = uxp./uzp*z0+xr;
        
        %Gather rays on sensor
        
        gatherer = gather_rays_nohist(xo,yo,npx,npy,dpx,dpy,nidx_1,midx_1,px);
        
        
        %Only gather rays if we've got more than 1 ray. It sounds
        %dumb, but the hist code doesn't work when there's 1 ray,
        %and I don't care about a single ray anyway.
        if nnz(gatherer)>1
            
            [r_outc{LF_index}, c_outc{LF_index}, v_outc{LF_index}] = ...
                build_A_matrix_sparse(gatherer,LF_index);
            if vis_prop %animation to visualize propagation after refraction
                
                for z = 0:zstep:zmax
                    set(0,'CurrentFigure',h3)
                    yo = uyp*z+yr;
                    xo = uxp*z+xr;
                    gatherer1 = gather_rays_nohist(xo,yo,npx,npy,dpx,dpy,nidx(1),midx(1),px);
                    
                    if nnz(gatherer1)>1
                        x_sensor = [0:npx-1]*dpx;
                        y_sensor = [0:npy-1]*dpy;
                        xmin = xg(1)+nidx(1)*px-tand(th_range)*zmax;
                        xmax = xg(end)+nidx(1)*px+tand(th_range)*zmax;
                        ymin = yg(1)+midx(1)*px-tand(ph_range)*zmax;
                        ymax = yg(end)+midx(1)*px+tand(ph_range)*zmax;
                        if dy_idx>1 && dx_idx>1
                            
                            imagesc(gatherer1,'XData',x_sensor,'YData',y_sensor)
                            hold on
                            axis image
                            axis([xmin xmax ymin ymax])
                        elseif dy_idx==0
                            stairs(x_sensor,gatherer1)
                            xlim([xmin xmax])
                            xlabel('\mum')
                        elseif dx_idx==0
                            stairs(y_sensor,gatherer1)
                            xlim([ymin ymax])
                            xlabel('\mum')
                        end
                        title(['indensity at z=',num2str(z),' for ',num2str(nrays),...
                            ' rays, \theta=',num2str(dth*(-qq+Q/2+.5)),...
                            ' \phi=',num2str(dph*(-pp+P/2+.5))])
                        hold off
                        pause(1/100)
                    end
                end
            end  %End propagation visualization
            
            if vis_sensor
                set(0,'CurrentFigure',h4)
                x_sensor = [0:npx-1]*dpx;
                y_sensor = [0:npy-1]*dpy;
                xmin = xg(1)+nidx(1)*px-tand(th_range)*z0;
                xmax = xg(end)+nidx(1)*px+tand(th_range)*z0;
                ymin = yg(1)+midx(1)*px-tand(ph_range)*z0;
                ymax = yg(end)+midx(1)*px+tand(ph_range)*z0;
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
                title(['indensity at z=',num2str(z0),' for ',num2str(nrays),...
                    ' rays, \theta=',num2str(dth*(-qq+Q/2+.5)),...
                    ' \phi=',num2str(dph*(-pp+P/2+.5))])
                hold off
                pause(1/24)
            end
        end
        
        %               %% Visualization
        if vis
            %Calculate global x and y index locations for rays
            scl = 1;
            nnr = xr/px+nidx(1);
            mmr = yr/px+midx(1);
            set(0,'CurrentFigure',h1)
            if dy_idx == 0
                %%
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
            pause(1/24)
        end
    end
else
    for mm = 1:max(1,M)
        ytstart = tic;
        %Get indices of region in y direction
        midx = (((mm-1)*dy_idx+1):(mm*dy_idx+1))';
        %Generate y vector in physical units in local coordinates
        %(i.e. origin shifts wich each new region
        yg = px*(midx-midx(1));
        %yr = dx_idx*px*(rand(nrays,1));
        for nn = 1:max(1,N)
            nidx = (((nn-1)*dx_idx+1):(nn*dx_idx+1))';
            xg = (nidx-nidx(1))*px;
            diff_crop = diffuser(midx,nidx);
            
            [NIDX, MIDX] = meshgrid(nidx,midx);
            
            Fx_crop = Fx(nidx,midx);
            Fy_crop = Fy(nidx,midx);
            for pp = 1:P
                ph = dph*(rand(nrays,1)-pp+P/2);    %Random phi points
                for qq = 1:Q
                    tstart = tic;
                    th = dth*(rand(nrays,1)-qq+Q/2); %random theta points
                    
                    %Generate random (x,y) positions within dx_idx*px wide
                    %square
                    yr = (dy_idx)*px*rand(nrays,1);
                    xr = (dx_idx)*px*rand(nrays,1);
                    
                    
                    %Calculate surface norm at random positions by
                    %interpolation within gradient
                    
                    %Get z and surface normal at each random (xr(i),yr(i)) pair
                    if dy_idx==0
                        Fxr = interp1(xg,Fx_crop,xr);  %Interpolate x gradient
                        Fyr = zeros(size(Fxr));
                        %zr = interp1(xg,diff_crop,xr);   %Interpolate surface
                    elseif dx_idx==0
                        Fyr = interp1(yg,Fy_crop,yr);  %Interpolate x gradient
                        Fxr = zeros(size(Fyr));
                        %zr = interp1(yg,diff_crop,yr);   %Interpolate surface
                    else
                        %zr = interp2(xg,yg,diff_crop,xr,yr);   %Interpolate surface
                        Fyr = interp2(xg,yg,Fx_crop',xr,yr);  %Interpolate x gradient
                        Fxr = interp2(xg,yg,Fy_crop',xr,yr);  %Interpolate y gradiet
                    end
                    
                    
                    %Refraction starts here ---------------------
                    [uxp, uyp, uzp] = refraction(Fxr, Fyr, th, ph, index);
                    %End refraction-------------
                    
                    %propagate to output plane by a distance z
                    [yo, xo] = propagation(uyp, uzp, z0, yr, uxp, xr);
                    %                 yo = uyp./uzp*z0+yr;
                    %                 xo = uxp./uzp*z0+xr;
                    
                    %Gather rays on sensor
                    
                    gatherer = gather_rays_nohist(xo,yo,npx,npy,dpx,dpy,nidx(1),midx(1),px);
                    
                    
                    %Only gather rays if we've got more than 1 ray. It sounds
                    %dumb, but the hist code doesn't work when there's 1 ray,
                    %and I don't care about a single ray anyway.
                    if nnz(gatherer)>1
                        
                        %make light field index
                        LF_index = P*Q*N*(mm-1)+P*Q*(nn-1)+Q*(pp-1)+qq;
                        %[A_row_index{LF_index}, A_vals{LF_index}] = find(gatherer(:));
                        %A_sub(:,LF_index) = gatherer(:);
                        %A_sub = build_A_matrix(A_sub,gatherer,LF_index);
                        %[r,c,v] = build_A_matrix_sparse(gatherer,LF_index);
                        %r_out = cat(1,r_out,r);
                        %c_out = cat(1,c_out,c);
                        %v_out = cat(1,v_out,v);
                        [r_outc{LF_index}, c_outc{LF_index}, v_outc{LF_index}] = ...
                            build_A_matrix_sparse(gatherer,LF_index);
                        if vis_prop %animation to visualize propagation after refraction
                            for z = 0:zstep:zmax
                                set(0,'CurrentFigure',h3)
                                yo = uyp*z./uzp+yr;
                                xo = uxp*z./uzp+xr;
                                gatherer1 = gather_rays_nohist(xo,yo,npx,npy,dpx,dpy,nidx(1),midx(1),px);
                                
                                if nnz(gatherer1)>1
                                    x_sensor = [0:npx-1]*dpx;
                                    y_sensor = [0:npy-1]*dpy;
                                    xmin = xg(1)+nidx(1)*px-tand(th_range)*zmax;
                                    xmax = xg(end)+nidx(1)*px+tand(th_range)*zmax;
                                    ymin = yg(1)+midx(1)*px-tand(ph_range)*zmax;
                                    ymax = yg(end)+midx(1)*px+tand(ph_range)*zmax;
                                    if dy_idx>1 && dx_idx>1
                                        
                                        imagesc(gatherer1,'XData',x_sensor,'YData',y_sensor)
                                        hold on
                                        axis image
                                        axis([xmin xmax ymin ymax])
                                    elseif dy_idx==0
                                        stairs(x_sensor,gatherer1)
                                        xlim([xmin xmax])
                                        xlabel('\mum')
                                    elseif dx_idx==0
                                        stairs(y_sensor,gatherer1)
                                        xlim([ymin ymax])
                                        xlabel('\mum')
                                    end
                                    title(['indensity at z=',num2str(z),' for ',num2str(nrays),...
                                        ' rays, \theta=',num2str(dth*(-qq+Q/2+.5)),...
                                        ' \phi=',num2str(dph*(-pp+P/2+.5))])
                                    hold off
                                    pause(1/100)
                                end
                            end
                        end  %End propagation visualization
                        
                        if vis_sensor
                            set(0,'CurrentFigure',h4)
                            x_sensor = [0:npx-1]*dpx;
                            y_sensor = [0:npy-1]*dpy;
                            xmin = xg(1)+nidx(1)*px-tand(th_range)*z0;
                            xmax = xg(end)+nidx(1)*px+tand(th_range)*z0;
                            ymin = yg(1)+midx(1)*px-tand(ph_range)*z0;
                            ymax = yg(end)+midx(1)*px+tand(ph_range)*z0;
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
                            title(['indensity at z=',num2str(z0),' for ',num2str(nrays),...
                                ' rays, \theta=',num2str(dth*(-qq+Q/2+.5)),...
                                ' \phi=',num2str(dph*(-pp+P/2+.5))])
                            hold off
                            pause(1/24)
                        end
                    end
                    
                    %               %% Visualization
                    if vis
                        %Calculate global x and y index locations for rays
                        scl = 1;
                        nnr = xr/px+nidx(1);
                        mmr = yr/px+midx(1);
                        set(0,'CurrentFigure',h1)
                        if dy_idx == 0
                            %%
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
                        pause(1/24)
                    end
                end
            end
        end
        ytend = toc(ytstart);
        tleft = (M*N*P*Q-LF_index)/(N*P*Q)*ytend;
        tmin = floor(tleft/60);
        tsec = mod(tleft/60,tmin)*60;
        waitbar(LF_index/M/N/P/Q,h5,...
            [num2str(100*(LF_index)/(M*N*P*Q)),'% done. ',num2str(M*N*P*Q),...
            ' total, ',num2str(ytend),' seconds per pass. ~',...
            num2str(tmin),':',num2str(tsec,'%.0f'),' to go.'])
    end
end

r_out = vertcat(r_outc{:});
c_out = vertcat(c_outc{:});
v_out = vertcat(v_outc{:});
A_sub = sparse(r_out,c_out,v_out,npx*npy,N*M*P*Q);
total_runtime = toc(tstart1);
fprintf([num2str(total_runtime),' seconds for ',num2str(P*Q*N*M*nrays),' rays\n'])

profile viewer


%fit polynomial locally

%calculate exact refraction angle

%propagate to output plane