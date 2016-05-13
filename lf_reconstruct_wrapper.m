solve_settings.tau = 30;
solve_settings.tv_iters = 8;
solve_settings.nneg = 1;
solve_settings.N = 7;
solve_settings.ToleranceA = .0001
bin = load('./Data/king_leaf_b.mat');
b = bin.b;

fprintf('data loaded\n')
recovered = lf_reconstruct(b,A,At,'3dtv','twist',settings,solve_settings);
