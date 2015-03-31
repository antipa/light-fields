function b = LF_project(x)
global A_sub_sparse
global intensity_point
global beta
b = norm((A_sub_sparse*x-intensity_point(:)),2)%;+beta*norm(diff(x),1);
end