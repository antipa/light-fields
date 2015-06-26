theta = th(1:10);
x = xo(1:10);

crop = diffuser_in(1,1:ceil(max(x)));
grad = gradient(diffuser_in(1,:));
gradx = interp1(0:length(grad)-1,grad,x);
gamma = atand(-gradx);
alpha = gamma - theta;
beta = asind(sind(alpha) ./ 1.5);
out = gamma - beta;

in2 = -out;
crop = diffuser_in(1,1:ceil(max(x)));
grad = -gradient(diffuser_in(1,:));
gradx = interp1(0:length(grad)-1,grad,x);
gamma = atand(-gradx);
alpha = gamma - in2;
beta = asind(sind(alpha) .* 1.5);
in2 = -(gamma - beta);
%in and th should be equal
