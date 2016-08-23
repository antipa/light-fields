data_type = 'simulation';

%TwIST
solve_settings.tau = 30;

%lsmr stuff
solve_settings.lam = 5e-1;
solve_settings.lsmr_tol = 1e-12;
solve_settings.dc = 1;   %Include DC bias offset in solution
solve_settings.precondition = 1; %Preconditioning: divide each column by its
%2-norm first, then multiply by preconditioning
%vector to recover solution

%TV stuff
solve_settings.tv_iters = 8;

solve_settings.nneg = 1;
solve_settings.N = 7;
solve_settings.MaxiterA = 500;
solve_settings.ToleranceA = .000005
solve_settings.wvlt = 'db9';
solve_settings.dwvlt = 'db9';
solve_settings.angweight = 1;
solve_settings.clipping = 1;

switch lower(data_type)
    case 'real'
        bin = load('../Output/king_queen.mat');
        b = bin.b-000;
        b(b<0) = 0;
    case 'simulation'
        fprintf('simulated sensor data. Loading light field')
        lf_in = load('/Users/nick.antipa/Documents/Light_field_data/spheresAndBlock/mat_files/15_degree_lf.mat');
        lf_orig = lf_in.lfr;
        lf_vec = lf_orig(:);
        b_raw = A_sub*lf_vec(1:size(A_sub,2));
        noise_var = .01;
        noise_add = randn(size(b))*max(b)*noise_var;
        b = b_raw+noise_add;
        b(b<0) = 0;
end
        
fprintf('data loaded\n')
recovered = lf_reconstruct(b,A_sub,At,'4dtv','twist',settings,solve_settings);
switch lower(data_type)
    case 'simulation'
        psnr_val = psnr(recovered,lf_orig);
        fprintf('psnr: %.2f',psnr_val)
end
