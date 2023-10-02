function [f_reg,lambda_opt,loss_reg] = lse_by_opt_con(L_operator,y,B,dx, normtype)
% inversion L_operator* f = y by constraint optimization with regularization
%{
          loss(x,lambda)   = \|L_operator*f - y \|^2 + lambda* f'*C*f 
                           = f'*A*f - 2*f'*L_operator'*y + y'*y + lambda* f'*C*f
with constraint: f is a density, i.e., sum f*dx  = 1, and f nonnegative  
%}
[m,n_x] = size(L_operator); 

 A = L_operator'*L_operator; 
       
% define loss function 
switch normtype
    case 'l2'     %% 1. Tikhonov regulariztion, norm_type = 'l2';
        Cmat      = eye(n_x);
        lossFn = @(f_lambda) norm(L_operator*f_lambda(1:end-1)-y)^2 + f_lambda(end)* f_lambda(1:end-1)'*Cmat *f_lambda(1:end-1);
    case 'L2'     %% 2. Tikhonov regulariztion, norm_type = 'L2';
        Cmat  = B; 
        lossFn = @(f_lambda) norm(L_operator*f_lambda(1:end-1)-y)^2 + f_lambda(end)* f_lambda(1:end-1)'*Cmat *f_lambda(1:end-1);
    case 'RKHS'   %% 3. Lcurve-rkhs regulariztion, norm_type = 'rkhs'       
        [V, eigL]  = eig(A,B);  % generalized eigenvalue   A V = B*V*eigL; V'*B*V =I;  V'*A*V = eigL; >>>  B_rkhs = inv(V*diag(eigL)*V')
        eigL       = real(diag(eigL));  [~,ind] = sort(eigL,'descend'); % V'*A*V = eigL; V'*B*V =I;
        eigL       = eigL(ind); V = V(:,ind);
        tol        = 1e-20;        eigL(eigL<tol)  = tol;   % if eigL < 1e-20, set it to be 1e-20;
        sqrtBinv = V*diag(real(eigL).^(1/2)); % sqrtBinv = sqrtm(V*diag(real(eigL))*V');    % Note: V is not unitary, thus, we cannot use V*sqrt(diag)*V'. Also, sqrtm has error when eigL is singular
        sqrtBinv = real(sqrtBinv);
        lossFn = @(f_lambda) norm(L_operator*f_lambda(1:end-1)-y)^2 + f_lambda(end)* norm(pinv(sqrtBinv)*f_lambda(1:end-1))^2;
end

Aeq = [ones(1, n_x),0];
beq = 1/dx;
Acs =[];
bcs =[];
f0 = ones(n_x+1, 1)/n_x/dx;
lb = [0*ones(n_x,1); 1e-14] ;
ub = [1e6*ones(n_x,1); max( eig(A) )];
nonlcon = [];
options = optimoptions('fmincon', 'MaxFunctionEvaluations',1000000,'StepTolerance',1e-6);
[f_lambda,loss_reg,exitflag_reg,output_reg,~,grad_reg,hessian_reg] = fmincon(lossFn, f0, Acs, bcs, Aeq, beq, lb, ub, nonlcon, options);
f_reg = f_lambda(1:end-1); lambda_opt = f_lambda(end);
end