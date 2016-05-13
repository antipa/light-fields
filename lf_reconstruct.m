function recovered = lf_reconstruct(b,A,At,model_type,solver_type,settings,solve_settings)

nx = settings.N;   %Number of light field x bins
ny = settings.M;   %Number of light field y bins
ntheta = settings.P; %Number of theta bins
nphi = settings.Q;  %Number of phi bins

% Make ordering vectors for regularizer reshaping
count= 0;
order_vec = [];

%prepares vector for reshaping into stack
for m = 1:nphi
    for n = 1:ntheta
        count = count+1;
        order_vec = cat(2,order_vec,count:ntheta*nphi:nx*ny*ntheta*nphi);
    end
end

count = 0;
deindex_vec = zeros(nx*ny*ntheta*nphi,1);

%Takes stack to vector
for q = 1:nx
    for p = 1:ny
        count = count+1;
        deindex_vec((count-1)*(ntheta*nphi)+1:(count)*(ntheta*nphi)) = count:nx*ny:nx*ny*ntheta*nphi;
    end
end

angle_deindex_vec = (1:nx*ny*ntheta*nphi)';   %Transform stack to vector
angle_index_vec = (1:nx*ny*ntheta*nphi)';  %vector to stack, switch nx,ntheta and ny,nphi



    switch lower(model_type)
        case 'angle_tv'
            A_func = @(x) A_tv(x,A,angle_deindex_vec);
            A_adj_func = @(x) A_adj_tv(x,At,angle_index_vec,ntheta,nphi,nx,ny);
            Phi_func = @(x) TVnorm_lf(x);
            Psi_func = @(x,th)  tvdenoise_lf(x,2/th,tv_iters);
            debias = 0;
        case '3dtv'
            if ~isfield(solve_settings,'tv_iters')
                tv_iters = 3;
            else
                tv_iters = solve_settings.tv_iters;
            end
            A_func = @(x) A_tv(x,A,deindex_vec);
            A_adj_func = @(x) A_adj_tv(x,At,order_vec,nx,ny,ntheta,nphi);
            Phi_func = @(x) TVnorm3d(x);
            Psi_func = @(x,th) tvdenoise3d(x,2/th,tv_iters,solve_settings.nneg);
            debias = 0;
            if ~isfield(solve_settings,'tau')
                error('solve_settings.tau must exist!')
            end
            
        case 'wavelet'
            dwtmode('per');
            if ~isfield(solve_settings,'wvlt')
                fprintf('using db9 default wavelet\n')
                wvlt = 'db9';
            else
                wvlt = solve_settings.wvlt;
            end
            if ~isfield(solve_settings,'dwvlt')
                fprintf('using db9 adjoint \n')
                dwvlt = 'db9';
            else
                dwvlt = solve_settings.dwvlt;
            end
            %wvlt = 'bior6.8';
            %dwvlt = 'rbio6.8';
            if ~isfield(solve_settings,'N')
                N = 4;
            else
                N = solve_settings.N;
            end
            [C,L] = wavedec2(rand(nx,ny),N,wvlt);
            wvltsize = length(C);
            %             for m = 1:ntheta*nphi
            %                 x = cat(1,x,wavedec2(recovered(:,:,m)',N,wvlt)');
            %             end
            A_func = @(x) A_wvlt2d(x,L,dwvlt,A,nx,ny,ntheta,nphi,deindex_vec);
            A_adj_func = @(x) A_adj_wvlt2d(x,At,wvlt,N,order_vec,ntheta,nphi,nx,ny,wvltsize);
            Phi_func = @(x) norm(x,1);
            Psi_func = @(x,tau) soft(x,tau);
            debias = 0;
            tau = 5;
        case 'hard_wavelet'
            dwtmode('per');
            %wvlt = 'db9';
            %dwvlt = 'db9';
            wvlt = 'bior6.8';
            dwvlt = 'rbio6.8';
            N = 4;
            [C,L] = wavedec2(rand(nx,ny),N,wvlt);
            wvltsize = length(C);
            %             for m = 1:ntheta*nphi
            %                 x = cat(1,x,wavedec2(recovered(:,:,m)',N,wvlt)');
            %             end
            A_func = @(x) A_wvlt2d(x,L,dwvlt,A,nx,ny,ntheta,nphi,deindex_vec);
            A_adj_func = @(x) A_adj_wvlt2d(x,At,wvlt,N,order_vec,ntheta,nphi,nx,ny,wvltsize);
            Phi_func = @(x) norm(x,1);
            Psi_func = @(x,tau) hard(x,tau);
            debias = 0;
            tau = 9000;
        case 'angle_wavelet'
            wvlt = 'haar';
            dwvlt = 'haar';
            N = 2;
            [C,L] = wavedec2(zeros(nphi,ntheta),N,wvlt);
            wvltsize = length(C);
            A_func = @(x) A_wvlt2d(x,L,dwvlt,A,ntheta,nphi,nx,ny,angle_deindex_vec);
            A_adj_func = @(x) A_adj_wvlt2d(x,At,wvlt,N,angle_index_vec,nx,ny,ntheta,nphi,wvltsize);
            Phi_func = @(x) norm(x,1);
            Psi_func = @(x,tau) soft(x,tau);
            debais = 0;
            tau = 30000;
        case 'lsmr'
            if ~isfield(solve_settings,'lam')
                lam = 1e-5;
            else
                lam = solve_settings.lam;
            end
            
    end

switch lower(solver_type)
    case('twist')
        if ~isfield(solve_settings,'ToleranceA')
            tolA = .0001;
        else
            tolA = solve_settings.ToleranceA;
        end
        
        if ~isfield(solve_settings,'MaxiterA')
            miter = 300;
        else
            miter = solve_settings.MaxiterA;
        end
        tau = solve_settings.tau;
        
        recovered =  TwIST(b,A_func,tau,...
            'AT', A_adj_func, ...
            'StopCriterion',1,...
            'Psi',Psi_func,...
            'Phi',Phi_func,...
            'Initialization',2,...
            'Debias',debias,...
            'MaxiterA',miter,...
            'Monotone',1,...
            'ToleranceA',tolA);
    case('lsmr')
        recovered = lsmr(A,b,lam,[],[],[],[],[],1e6);
end


% Reshape vector to 3D light field
switch lower(model_type)
    case 'angle_tv'
        recovered_ang = recovered;
        recovered_vec = recovered_ang(:);
        recovered = reshape(recovered_vec(order_vec),[nx,ny,ntheta*nphi]);
    case 'angle_wavelet'
        recovered_ang = recovered;
        recovered_vec = recovered_ang(:);
        recovered = reshape(recovered_vec(order_vec),[nx,ny,ntheta*nphi]);
    case {'wavelet','hard_wavelet'}
        rec_w = zeros(nx,ny,ntheta*nphi);
        
        decsize = length(recovered)/ntheta/nphi;
        for m = 1:ntheta*nphi
            
            
            vind = (m-1)*decsize+1:m*decsize;
            rec_w(:,:,m) = waverec2(recovered(vind),L,wvlt);
        end
        recovered = rec_w;
    case 'lsmr'
        recovered = reshape(rec_backup(order_vec),ny,nx,ntheta*nphi);
end
