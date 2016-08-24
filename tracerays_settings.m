%%
%%   This file contains the settings for the project
%%   You should make a copy of tracerays_settings.m.template as tracerays_settings.m and make your changes there
%%
%'forward'
%'backward'
direction = 'forward';   %'forward' or 'backward'. Forward traces rays from light field bundles
                         % to sensor. Backward traces from pixels back to
                         % light field space. This is potentially useful
                         % if we want to ditch the discretization in LF
                         % space.
sync = false;
%TODO:adjust syncing for 4d
switch direction
    case 'forward'
        %must make vis_prop = 0
        settings.paraxial = false;   %Use paraxial approximation to speed up refraction.
        settings.useParfor = false;    %use parallel for loop. Faster if set to true, but no plots can be made
        settings.code_profiler = false;    %Leave set to false by default. Useful for speeding up new modules
        settings.interp_method = 'linear';  %Interpolation method for gradient. Nearest is quick!
        if settings.useParfor  %Override code profiler if using parfor
            settings.code_profiler = false;
        end
        %settings.gather_interp = true;   %No longer used
        settings.strength = 1;   %z-direction (amplitude) scaling applied to raw diffuser surface data
        settings.diff_upsample = 0;  %Upsample diffuser. If set to 0, don't upsample, if ~=0, this is the resampling factor
        settings.vis = 0;   %Visualize surface, normals, and refracted rays. Good for slides.
        settings.vis_prop = 0;   %Visualize gatherer as a function of z. 
        settings.save_prop = [];   %Filename for outputing propagation video. Put in a string here, and you'll get a gif of propagation saved to that filename.
        settings.vis_sensor = 1;   %Visualize image created in each pass -- this is a 2d version of each column in A
        %Range and step size for propagation movie
        settings.zmax = 1500;   %Max prop distance for vis_prop
        settings.zstep = 100;    %step size
        settings.preload_diff = 0;   %Set to 0 to load diffuser file specified in newTraceRays, set to 1 if you loaded one yourself and it's in memory already
        settings.vis_hist = 0;   %Histogram of sensor data
        settings.vis_xcorr = 0;   %Show autocorrelation of caustic pattern--closer to delta is good for inversion 
        settings.waitbar = 0;   %Waitbar that estimates time left. Doesn't work with parfor
        %Setup tracing grid
        %number of Y points use 1 for 1d case. This is the number of
        %spatial samples you'll get in your resulting light field. 
        settings.M = 64;
        % number of X points
        settings.N = 64;
        %number of phi points (angle in y (M) direction)   use 1 for 1d.
        %This is the number of subaperture images.
        settings.P = 12;
        %number of theta points (angle in x direction)
        settings.Q = 12;
        settings.nrays = 1e5;   %Number of rays to cast over each unique (x,y,theta,phi). 
                                %Total number of rays is nrays*P*Q*N*M
        %how far along x to go in same units as pixels (micron). Not index!
        settings.x_range = 2.5*1024;
        %This will be divided into N steps
        %in M steps. Use 0 for 1d.
        settings.y_range = 2.5*1024;
        %Setup sensor parameters for a sensor that is the same size as the diffuser
        settings.npx = 800;
        %for 1d use 1
        settings.npy = 800;
        %z0 = 228.6
        %settings.z0 = 350; Works with strength-.446
        settings.z0 = 500;  %Propagation distance to caustic plane
        %how far in angle to go in degrees
        %Divided into P steps.
        settings.ph_range = 5;   %Angle spraed in degrees. Phi is in y direction
        settings.th_range = 5;   %x direction angle (theta) spread
        settings.indexEnv = 1;   %Environmental index of refraction (between surface and sensor)
        settings.indexDiff = 1.68;    %Diffuser index of refraction
        settings.dn = settings.indexDiff-settings.indexEnv;  %Store difference since this is what's used in equations
        if sync  %Setup corresponding parameters for backward ray tracing. This traces from sensor back to light field space. 
            % Currently, forward is the way to go. 
            settings.pixelSize = x_range ./ npx; %pixel size in physical units
            settings.sensorSizeX = npx; %given in pixels
            settings.sensorSizeY = npy;
            %2.5 is for strength 50
            settings.thetaSpread = th_range + 2.5; %in degrees
            settings.phiSpread = ph_range + 2.5;
            settings.rays = nrays * (M*N*P*Q);
            
            %distance to diffuser
            settings.z = z0;
            
            %diffuser strength
            settings.strengthB = -strength;
            settings.gridX = N;
            settings.gridT = Q;
            settings.xRange = [0 x_range];
            settings.tRange = [(th_range ./ 2) (-th_range ./ 2)];
        end
    case 'backward'
        pixelSize = 2; %pixel size in physical units
        sensorSizeX = 50; %given in pixels
        sensorSizeY = 50;
        thetaSpread = 2; %in degrees
        phiSpread = 2;
        rays = 1000000;
        
        %distance to diffuser
        z = 100;
        %diffuser strength
        strengthB = -50;
        %second index of refraction
        indexEnv = 1;
        
        %first index of refraction
        indexDiff = 1.5;
        
        gridX = 10;
        gridT = 3;
        xRange = [0 pixelSize*sensorSizeX];
        tRange = [.5 -.5];
        gridP = 3;
        pRange = [.5 -.5];
        gridY = 10;
        yRange = [0 pixelSize*sensorSizeY];
        if sync
            %must make vis_prop = 0
            useParfor = false;
            strength = -strengthB;
            diff_upsample = 0;
            vis = 0;
            vis_prop = 0;
            save_prop = 0;
            vis_sensor = 0;
            %Range and step size for propagation movie
            zmax = 1000;
            zstep = 10;
            %Setup tracing grid
            %number of Y points use 1 for 1d case
            M = gridY;
            % number of X points
            N = gridX;
            %number of phi points (angle in y (M) direction)   use 1 for 1d
            P = gridP;
            %number of theta points (angle in x direction)
            Q = gridT;
            nrays = round(rays ./ (N*M*P*Q));
            %how far along x to go in same units as pixels (micron)
            x_range = pixelSize * sensorSizeX;
            %This will be divided into N steps
            %in M steps. Use 0 for 1d.
            y_range = pixelSize * sensorSizeY;
            %Setup sensor parameters for a sensor that is the same size as the diffuser
            npx = sensorSizeX;
            %for 1d use 1
            npy = sensorSizeY;
            %z0 = 228.6
            z0 = z;
            %how far in angle to go in degrees
            %Divided into P steps.
            ph_range = max(pRange) - min(pRange);
            th_range = max(tRange) - min(tRange);
        end
end
