function basis = spline_get_basis(lb, rb, knot_num, deg,  dx)
% This function generates the spline basis function


% L = 10;
% knot_num = 10;
% deg = 2;
% dx = 0.01;
if nargin == 4
    dx = 0.001;
end

u = linspace(lb,rb,knot_num+1);
u = [lb*ones(1,deg), u, rb*ones(1,deg)];
xgrid = (lb:dx:rb)';
N = cell(deg+1, 1);

for p = 1:deg + 1
    for i = 1 : length(u) - p
        if p == 1
            if i ~= length(u) - deg -1
                N{p}(:,i) = (xgrid>=u(i)).*(xgrid<u(i+1));
            else
                N{p}(:,i) = (xgrid>=u(i)).*(xgrid<=u(i+1));
            end
        else
            if u(i+p-1) == u(i);        temp1 = 0;
            else;                       temp1 = (xgrid-u(i))/(u(i+p-1)-u(i));       end
            
            if u(i+p) == u(i+1);        temp2 = 0;
            else;                       temp2 = (u(i+p)-xgrid)/(u(i+p)-u(i+1));     end
            
            N{p}(:, i) = temp1.*N{p-1}(:, i) + temp2.*N{p-1}(:, i+1);
        end
    end
end

temp = N{end};

dim = knot_num + deg;
basis = cell(dim, 1);

for i = 1:dim
    f = fit(xgrid, temp(:, i), 'linearinterp');
    basis{i} = @(x) reshape(f(x), size(x));
end


end
