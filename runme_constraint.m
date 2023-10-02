% This is the main function of this project
clc; close all; clear all;
add_mypaths_discrete;    % get SAVE_DIR = local dir for saving data
rng(1)
%% Load system settings
sysInfo = system_settings();
% sysInfo = update_system_settings(sysInfo);

%% Generate Data
lb = sysInfo.lb;
rb = sysInfo.rb;
T = sysInfo.T;
tn = sysInfo.tn;
xn = sysInfo.xn;
L_operator = sysInfo.L_operator;
tgrid = sysInfo.tgrid;
xgrid = sysInfo.xgrid;
dt = sysInfo.dt;
dx = sysInfo.dx;
phi = sysInfo.phi;

%% Get regression matrix A, 
% to get vector b later for different f
A = L_operator'*L_operator;  

%% Get rho and L2(rho) basis matrix B 
rho = sum(L_operator);  rho = rho/(sum(rho)*dx);  % normalize, does not seem necessary for this example, maybe other examples
figure; % plot the exploration measure 
plot(xgrid, rho,'linewidth',1); xlabel('u');ylabel('rho');
B   = diag(rho);

%% analysis function space of identifability 
[V_A,eigA,V_AB, eigAB,r]= EigenAB_fsoi(A,B); 

%% Get regression vector for each f. 
nsr   = 0.01;    % noise to signal ratio    
           % -- issue when nsr=0: the optimal lambda=0, but numerical error in inversion prevents us from get to it.  
           % Solution: add lambda =0 estimator, and select between lambda_opt and 0 by min-loss ( a factor (1e2) to be robust) >> estimator

% f_true from functions 
f_true_func = @(x) (sin(x-6)).^2+1;   % True f
%  f_true_func = @(x) (x-6).^2+1;   % True f
% f_true_func = @(x) 15*((sin(x-6)).^2-3);   % True f
% f_true_func = @(x) 15*((sin(2*x-6)).^2-3);

norm = integral(f_true_func, lb, rb);
f_true_func = @(x) f_true_func(x) / norm;


f_true = f_true_func(xgrid);

% f_true insider and outside of FSOI: 
% f_true = V_AB(:,4);   f_true = V_AB(:,5);    % ind>5 ( i.e, eigAB<1e-9): rkhs not good, l2+L2 can slightly tolerate more, bc. not using rkhs inversion).



f_true_norm = sum(f_true.^2.*rho); 

y_true = L_operator*f_true; 
y_true_norm = sqrt(sum(y_true.^2)*sysInfo.dt/sysInfo.T);
y      = y_true + y_true_norm*nsr*randn(tn, 1);
b      = L_operator'*y;

%% estimate by constraint optimization: fmincon
lambda_reg = 0.0001;
% fun = @(c) c'*A*c - 2*b'*c + b'*pinv(A)*b + lambda * c'*C*c;
C = eye(xn);
fun_reg = @(c) c'*A*c - 2*b'*c + b'*pinv(A)*b + lambda_reg * c'*C*c;
Aeq = ones(1, xn);
beq = 1/dx;
Acs = -eye(xn);
bcs = zeros(xn, 1);
f0 = ones(xn, 1)/xn/dx;
lb = [];
ub = [];
nonlcon = [];
options = optimoptions('fmincon', 'MaxFunctionEvaluations',1000000);
[f_reg,fval_reg,exitflag_reg,output_reg,~,grad_reg,hessian_reg] = fmincon(fun_reg, f0, Acs, bcs, Aeq, beq, lb, ub, nonlcon, options);

fun = @(c) c'*A*c - 2*b'*c;
[f,fval,exitflag,output,~,grad,hessian] = fmincon(fun, f0, Acs, bcs, Aeq, beq, lb, ub, nonlcon, options);
% plot the result
lgnd      = {'f no reg','f reg','True'}; 
est_array = [f, f_reg, f_true];
plot_estimators(est_array,lgnd,xgrid,rho); 
plot_estimators_image(est_array,tgrid,L_operator,y,y_true,lgnd); 



%%
% figure;plot(f_reg);hold on;plot(f_true)
% figure;plot(L_operator*f_reg);hold on;plot(y);plot(y_true)


