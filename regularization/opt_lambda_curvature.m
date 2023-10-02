
function lambda_opt = opt_lambda_curvature(res_seq,reg_seq, lambda_seq,plotON)
% get optimal lambda using curvature of (x(t), y(t))

reg_seq   = log10(reg_seq)/2;
res_seq   = log10(res_seq)/2;
kappa     = curvature_approx( res_seq,reg_seq, lambda_seq);
len           = numel(lambda_seq);
indkeep       = 1:len; 
[max_kappa, max_ind] = max(kappa(indkeep));
lambda_opt   = lambda_seq(max_ind); % attains max curvature of L curve

if plotON
    figure;
    subplot(121);plot(res_seq, reg_seq, 'o-');hold on; scatter(res_seq(max_ind), reg_seq(max_ind), 'd','filled');
    title('L curve');  ylabel log_{10}(||x||_{B});  xlabel log_{10}(||Ax-b||);hold on;
    subplot(122);semilogx(lambda_seq, kappa,'o-'); hold on; scatter(log10(lambda_opt),kappa(max_ind),'d','filled');xlabel('lambda');
    title('Curvature of L curve');hold off;
end
end