function [lambda] = Opt_lambda_Quanjun(E, R, all_lambda, plotON, titl)
% find the maximal curvature of the curve (log(E), log(R))
% 1. we fit this curve using FIT method of MATLAB
% 2. Use uniform mesh for the E axis and compute the curvature using
% discrete mesh
% 3. We can get the coordinate for E for the maximal curvature point. But
% then we have to find the corresponding lambda for this perticular E. The
% method is to fit another curve for lambda and E.
% We can do this is because (log(E), log(R)) is continuous and
% (lambda, E) is monotone (!)
%% 0. Make  sure data are all positive
% Due to numerical issue, small lambda will cause res < eps. Hence we
% truncate the range for the lambda

% N = length(all_lambda);
% 
% I = 1:N;
% temp = I(E < 10*eps);
% if length(temp) == 0
%     ind = 0;
% else
%     ind = temp(end);
% end
% 
% temp = I(R < 10*eps);
% if length(temp) == 0
%     tail = N+1;
% else
%     tail = temp(1);    
% end
% 
% E = E(ind+1:tail-1);
% R = R(ind+1:tail-1);
% all_lambda = all_lambda(ind+1:tail-1);

%% 1. Fit the curve (log(E), log(R))
E(E<1e-50) = 1e-50;  E(E>1e20) = 1e20;   % avoid to NAN or inf  
R(R<1e-50) = 1e-50;  R(R>1e20) = 1e20;   % avoid to NAN or inf  
x = log10(E);
y = log10(R);

f = fit(x, y, 'smoothingspline');

% equally spacing over x axis
xgrid = linspace(min(x), max(x), 200);
dx = xgrid(2) - xgrid(1);
zeta = xgrid';
ita = f(xgrid);

% compute the curvature using curvature formula of a function
% k = y''/(1 + y'^2)^(3/2)

g = gradient(ita, dx);
gg = gradient(g, dx);
k = gg ./ ((1 + g.^2).^(3/2));
[~, ind] = max(k);




%% choose the best lambda from zeta_opt (old)
% h = fit(x(2:end), log10(all_lambda(2:end))', 'smoothingspline');
% 
% if ind == 1
%     lambda = 0;
% else
%     lambda = 10^h(zeta(ind));
% end

%% choose the best lambda from zeta_opt (old)
h = fit(x, log10(all_lambda)', 'smoothingspline');


lambda = 10^h(zeta(ind));
if ~ isfinite(lambda)
    lambda = 0;
end

% a version using provided input sequence of lambda
% not recommended.
% We can just use the eigenvalues of A as input lambda sequence
% and fit a best lambda from the curves.


% I = sum(x < zeta(ind));
% lambda_opt = all_lambda(I);

%% plot the curvatures
if plotON
    %%
    figure;
    subplot(4,2,[1,3])
    plot(zeta, ita, '-*');hold on;plot(zeta(ind), ita(ind), 'or', 'MarkerSize',10, 'LineWidth', 2);title('THE L-CURVE')
    ylabel log_{10}(||x||_{B});  xlabel log_{10}(E);
   % subplot(4,2,2)
   % plot(zeta, g, '-*');hold on;plot(zeta(ind), g(ind), 'or', 'MarkerSize',10, 'LineWidth', 2);title('df/dx')
   % xlabel log_{10}(E);
   % subplot(4,2,3)
   % plot(zeta, gg, '-*');hold on;plot(zeta(ind), gg(ind), 'or', 'MarkerSize',10, 'LineWidth', 2);title('d^2f/dx^2')
   % xlabel log_{10}(E);
    subplot(4,2,[2,4])
    plot(xgrid, k, '-*');hold on;plot(zeta(ind), k(ind), 'or', 'MarkerSize',10, 'LineWidth', 2);title('curvature')
    xlabel log_{10}(E);
    legend('data',['log_{10}(E_0) = ', num2str(zeta(ind))],'Location','best')
    %     sgtitle('Use fitted curve and evaulate curvature on equally spacing xgrids')
    subplot(4,2,5);scatter(log10(all_lambda), log10(R));xlabel log_{10}(\lambda);ylabel log_{10}(R);
	subplot(4,2,6);scatter(log10(all_lambda), log10(E));xlabel log_{10}(\lambda);ylabel log_{10}(E);

    subplot(4,2,[7,8])
    plot(x, h(x), '-*');hold on;plot(zeta(ind), log10(lambda), 'or', 'MarkerSize',10, 'LineWidth', 2)
    xlabel('log_{10}(E)')
    ylabel('log_{10}(\lambda)')
    legend('data',['\lambda_0 = ', num2str(lambda)],'Location','best')
    max_temp = compose('%1.1e', max(all_lambda));
    min_temp = compose('%1.1e', min(all_lambda));
    title(['\lambda \in (', min_temp{1}, ',', max_temp{1}, ')'])
    sgtitle(titl)
    %%
end
end