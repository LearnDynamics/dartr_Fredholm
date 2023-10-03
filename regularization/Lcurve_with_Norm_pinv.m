function [x_reg,lambda_opt,E_regu] = Lcurve_with_Norm_pinv(Abar, bbar, B, titl, plotON,normtype)
% L-curve regularization with a given norm and specified range for parameter
%                (A + lambda B)c = b
%% evaluate error/loss and regularization term
          % lambda in range adapative to eigenvalues
          % loss evaluated by norm(sqrt(A)*c- sqrt(A)\b):  no need to estimate the constant 
          % x_lambda solved by pinv or lsqmininorm
% [~, eigA]  = eig(A);  % generalized eigenvalue   A V = B*V*eigL; V'*B*V =I;  V'*A*V = eigL; >>>  B_rkhs = inv(V*diag(eigL)*V')
% eigA       = real(diag(eigA)); 

[V, eigL]  = eig(Abar,B);  % generalized eigenvalue   A V = B*V*eigL; V'*B*V =I;  V'*A*V = eigL; >>>  B_rkhs = inv(V*diag(eigL)*V')
eigL       = real(diag(eigL));  [~,ind] = sort(eigL,'descend'); % V'*A*V = eigL; V'*B*V =I;
eigL       = eigL(ind); V = V(:,ind);
tol  = 1e-20; 
eigL(eigL<tol)  = tol;   % if eigL < 1e-20, set it to be 1e-20; 


% Some observations:
% 1. c'Ac - 2b'c + b'A^-1b = (Ac - b)'A^-1(Ac-b)
% 2. (A + lambda B)c = b  ====> Ac-b = -lambda * B * c
N = 1000;
if plotON==1; N=200; end
lambda_seq = 10.^(linspace(log10(min(eigL(eigL>1e-12)))-2, log10(min(1,max(eigL))), N));
%lambda_seq = 10.^(linspace(-16, 10, N));
len = length(lambda_seq);
E = zeros(len,1);    % error/loss 
R = zeros(len,1);    % regularization norm^2


% % % % change back to least sqaures x'*A*x -2*x'*b +C = | A2 x - b2|^2, such that A=  A2'*A2, b = A2'*b2
% [U,eigA] = eig(A);   
% eigA     = real(diag(eigA));  [~,ind] = sort(eigA,'descend'); eigA = eigA(ind); U = U(:,ind);
% A2       = U*diag(real(eigA).^(1/2)); 
% b2       =pinv(A2)*b; % lsqminnorm(A2,b);  % A2\b;      % or lsqminnorm(A2,b); 

 A2 = sqrtm(Abar);   b2 = pinv(A2)*bbar;    % so that E = norm(A2*x_l - b2);   %%% this gives better results for L2 in invLap

% compare_computE(lambda_seq,A,B,b); 

for ll = 1:len
    lambda1 = lambda_seq(ll);
    x_l     =  lsqminnorm(Abar+lambda1*B, bbar);% pinv(Abar+lambda1*B)*b; %
%     E(ll) = (Abar*x_l-b)'*lsqminnorm(Abar, Abar*x_l-b);
%     E(ll) = x_l'*Abar*x_l -2*b'*x_l + b'*(pinv(Abar) * b);
    E(ll) = norm(A2*x_l - b2); % This is the best method
    R(ll) = x_l'*B*x_l;
end
R = sqrt(R); 


if 0   % plot the scattered L-curve points
    figure;
    subplot(2,1,1);scatter(log10(lambda_seq), log10(R));xlabel log_{10}(\lambda);ylabel log_{10}(R);
	subplot(2,1,2);scatter(log10(lambda_seq), log10(E));xlabel log_{10}(\lambda);ylabel log_{10}(E);
end


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
x_reg =  real(pinv(Abar+lambda_opt*B)*bbar); % (A+lambda_opt*B)\b;


%% when the nsr = 0: no regularization will be the best 
x_no_regu = lsqminnorm(Abar,bbar); 
E_no_regu = norm(A2*x_no_regu - b2);
E_regu    = norm(A2*x_reg - b2); 
if E_regu >  1e1*E_no_regu  && lambda_opt <1e-12 
    x_reg = x_no_regu; E_regu = E_no_regu; fprintf(['No regularization,',normtype,'\n']);
end


end


function compare_computE(lambda_seq,A,B,b)
% compare two ways computing E: (1) by c'*A*c - 2*c'*b + const; (2) norm(sqrt(A)*c- sqrt(A)\b)
% conclusion: when the noise presents or when lambda is small.  
len = length(lambda_seq);
E2 = zeros(len,1);    % error/loss 
R2 = zeros(len,1);    % regularization norm^2
E  = zeros(len,1);
% Approach 2: loss evaluated by norm(sqrt(A)x-b): does not need to estimate the constant 
% Approach 1: direct evaluation of loss: with const estimated 
A2    = sqrtm(A); 
b2    = pinv(A2)*b;
const = b'*(pinv(A) * b);  
for ll = 1:len
    lambda1 = lambda_seq(ll);
    x_l     = pinv(A+lambda1*B)*b; % lsqminnorm(A+lambda1*B, b);
%     E(ll) = (A*x_l-b)'*lsqminnorm(A, A*x_l-b);
%     E(ll) = x_l'*A*x_l -2*b'*x_l + b'*(pinv(A) * b);
    E2(ll) = norm(A2*x_l - b2)^2; % This is the best method
    E(ll)  = x_l'*A*x_l -2*b'*x_l + const;
    R2(ll) = x_l'*B*x_l;
end
figure; 
loglog(lambda_seq,E,lambda_seq,E2,'-.','linewidth',1);xlabel('lambda'); ylabel('loss values');
legend('E with estimated const','E= norm(sqrt(A)-b)'); 
end