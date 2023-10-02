function c = curvature_approx(x, y, lambda)
% compute curvature of the curve (x(t),y(t))
%  using the formula: (x(t), y(t))
%                      x' y'' - y'x''
%            kappa = ------------------- 
%                     (x'^2+y'^2)^(3/2)
n = numel(lambda);
x_p  = zeros(n, 1);
x_pp = zeros(n, 1);
y_p  = zeros(n, 1);
y_pp = zeros(n, 1);
c    = zeros(n, 1);
for i = 2:n-1
    x_p(i)  = (x(i)- x(i-1))/(lambda(i)- lambda(i-1));
    x_pp(i) = (x(i-1) + x(i+1) - 2 * x(i))/ (lambda(i+1)- lambda(i-1));
    y_p(i)  = (y(i)- y(i-1))/(lambda(i)- lambda(i-1));
    y_pp(i) = (y(i-1) + y(i+1) - 2 * y(i))/ (lambda(i+1)- lambda(i-1));
    c(i)    = 2*abs(x_p(i) * y_pp(i) - x_pp(i) * y_p(i))/((x_p(i)^2 +y_p(i)^2)^(3/2));
end
end