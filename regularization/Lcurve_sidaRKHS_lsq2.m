function [x_reg,lambda_opt,E_regu] = Lcurve_sidaRKHS_lsq2(A,b,Bmat,titl,plotON)
% L-curve regularization with RKHS using lsqminnorm minimizing rkhs-norm,
% but without computing Brkhs, which involves inversion of singular/ill-conditioned matrices (Lcurve_with_Norm_lsq2.m)

%{
[V, eigL]  = eig(Abar,Bmat);  % generalized eigenvalue   A = B*V*eigL; V'*B*V =I;  V'*A*V = eigL; >>>  B_rkhs = inv(V*diag(eigL)*V')
eigL       = real(diag(eigL));  [~,ind] = sort(eigL,'descend');
eigL       = eigL(ind); V = V(:,ind);
B_rkhs     = pinv(V*diag(real(eigL))*V');   % pinv(Abar);  B_rkhs= inv(Abar) when Bmat=Id;
% [U,S,V] = svd(B_rkhs); sqrtBinv = U*pinv(sqrtm(real(S)))*U'; % sqrtBinv = sqrtm(pinv(B_rkhs)); 
% [U,S] = eig(B_rkhs); sqrtBinv = U*pinv(sqrtm(real(S)))*U'; 
>>> sqrtBinv = V*diag(sqrt(real(eigL)))*V'
%}

[V, eigL]  = eig(A,Bmat);  % generalized eigenvalue   A V = B*V*eigL; V'*B*V =I;  V'*A*V = eigL; >>>  B_rkhs = inv(V*diag(eigL)*V')
eigL       = real(diag(eigL));  [~,ind] = sort(eigL,'descend'); % V'*A*V = eigL; V'*B*V =I;
eigL       = eigL(ind); V = V(:,ind);
tol  = 1e-20; 
eigL(eigL<tol)  = tol;   % if eigL < 1e-20, set it to be 1e-20; 

% B_rkhs = pinv(V*diag(eigL)*V');             % round off error is too large when eigL is almost singular
% sqrtBinv = sqrtm(V*diag(real(eigL))*V');    % Note: V is not unitary, thus, we cannot use V*sqrt(diag)*V'. Also, sqrtm has error when eigL is singular
sqrtBinv = V*diag(real(eigL).^(1/2));  
sqrtBinv = real(sqrtBinv); 

%% === the rest are the same as 
b1       = sqrtBinv'*b;
A1       = sqrtBinv'*A * sqrtBinv; 
n = length(b1);
%% Treatment for lambda range and LossFn: so that the Lcurve to work the best
%   1. lambda  in a range depending on the eigenvalues
%   2. compute Loss using norm(sqrt(A)x-sqrt(A)\b), to avoid estimating const 
% Some observations:
% 1. c'Ac - 2b'c + b'A^-1b = (Ac - b)'A^-1(Ac-b)
% 2. (A + lambda B)c = b  ====> Ac-b = -lambda * B * c
N  = 1000;
% lambda_seq = 10.^(linspace(-20, log10(max(eigL)), N)); 
lambda_seq = 10.^(linspace(log10(min(eigL(eigL>1e-12)))-2, log10(max(eigL)), N));
% lambda_seq = linspace(min(eigL)*1e-3, max(eigL(4:end)), N);
len = length(lambda_seq);
E   = zeros(len,1);
R   = zeros(len,1);

% % % % change back to least sqaures x'*A*x -2*x'*b +C = | A2 x - b2|^2, such that A=  A2'*A2, b = A2'*b2
% [U,eigA] = eig(A);   
% eigA     = real(diag(eigA));  [~,ind] = sort(eigA,'descend'); eigA = eigA(ind); U = U(:,ind);
% A2       = U*diag(real(eigA).^(1/2)); 
% b2       =pinv(A2)*b; % lsqminnorm(A2,b);  % A2\b;      % or lsqminnorm(A2,b); 

 A2 = sqrtm(A);   b2 = pinv(A2)*b;    % so that E = norm(A2*x_l - b2);   %%% this gives better results for L2 in invLap

for ll = 1:len
    lambda1 = lambda_seq(ll);
    % for each lambda compute c_lambda = (A+lambda)\b by lsqminnorm
     A_l      = A1+lambda1*eye(n);
     x_rkhs   = lsqminnorm(A_l,b1);
     x_l      = sqrtBinv*x_rkhs;
    %  % A_l_root = sqrtm(A_l);   x_l     = sqrtBinv*lsqminnorm(A_l_root, pinv(A_l_root)*b1);  
    E(ll) = norm(A2*x_l - b2); %This is the best method
    R(ll) = x_rkhs'*x_rkhs;
end
R  = sqrt(R); 

%% L-curve
if ~exist('curvatureType','var');  curvatureType = '3ptCircle'; end
switch curvatureType
    case 'fit_curve'
        lambda_opt = Opt_lambda_Quanjun(E, R, lambda_seq, plotON, titl);
    case '3ptCircle'      % use 3point circle radius
        lambda_opt = Opt_lambda_curvature3pt(E,R, lambda_seq, plotON,titl);
    otherwise
        lambda_opt = opt_lambda_curvature(E,R, lambda_seq, plotON);
end

%% select from pinv and lsqminnorm
% if rank(A1) < length(A1(1,:)) && lambda_opt< min(eigL)
%    x_pinv = sqrtBinv*pinv(A1+lambda_opt*eye(n)) * b1; % lsqminnorm(A+lambda_opt*B, b);    % % >>> > lsqminnorm and pinv are similar here bc. A+lambda_opt*B is well-conditioned.
%    x_reg  = x_pinv; 
% elses
    x_lsq  = sqrtBinv*lsqminnorm(A1+lambda_opt*eye(n), b1);
    x_reg  = x_lsq;
% end
% if  norm(A2*x_lsq - b2)< norm(A2*x_pinv - b2)
%     x_reg = x_lsq; 
% else;      x_reg = x_pinv; 
% end  


%% deal the case when nsr = 0: compare without regularization
x_no_regu = lsqminnorm(A,b); 
E_no_regu = norm(A2*x_no_regu - b2);
E_regu    = norm(A2*x_reg - b2); 
if E_regu >  1e1*E_no_regu  && lambda_opt <1e-12 
    x_reg = x_no_regu; E_regu = E_no_regu; fprintf(['No regularization: RKHS\n']);
end

end
