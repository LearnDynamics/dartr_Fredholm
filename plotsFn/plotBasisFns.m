function plotBasisFns(xgrid,basis_funs, title_str)
% plot the basis functions; 

nBasis  = length(basis_funs);
%xmin  = min(xgrid); xmax = max(xgrid); 
%xgrid = xmin:(xmax-xmin)/1000:xmax; 

fval = zeros(length(xgrid), nBasis);

for i = 1: nBasis   
    fval(:, i) = arrayfun(basis_funs{i}, xgrid);    
end
%{
if strcmp(inferInfo.basisType,'Bspline')
    figure; plot(xgrid, sum(fval,2), 'r'); title('Test partition unity: sum = 1');
end
%}

plot(xgrid, fval,'linewidth', 2);  title(title_str); 

end