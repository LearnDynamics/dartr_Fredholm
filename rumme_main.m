% Adaptive RKHS regularization for discrete Fredhold integral equation
%{
 Dicrete-Fredholm integral equation:
      \int_lb^rb K(t,x) f(x)dx + noise = y(t)
 Goal: Given y(t) at discrete-times, to estimate f
%}
%{
In general: weighted Deconvolution in the form of inversion
               L f = y,   size(L) =  n_y x n_u       size(f) = n_u 
Solution: least square with RKHS-regularization 
Key: space of identifiability, exploration measure, RKHS regularization 
@Copyright: Fei Lu, feilu@math.jhu.edu. 2023/10/2
%}



clc; close all; clear all;
add_mypaths_discrete;                    % get SAVE_DIR = local dir for saving data
rng(1)
%% Load system settings
sysInfo    = system_settings();
L_operator = sysInfo.L_operator;
xn    = sysInfo.xn;
tn    = sysInfo.tn;
tgrid = sysInfo.tgrid;
xgrid = sysInfo.xgrid;
dx    = sysInfo.dx;

exp_poly = sysInfo.kernel_type;          % decay of the spectrum of the operator of inversion (L_G)

%% Get regression matrix A, 
% to get vector b later for different f
A = L_operator'*L_operator;  

%% Get rho and L2(rho) basis matrix B 
rho = sum(L_operator);  rho = rho/(sum(rho)*dx);  % normalize, does not seem necessary for this example, maybe other examples
figure; % plot the exploration measure 
plot(xgrid, rho,'linewidth',1); xlabel('u');ylabel('rho');
if strcmp(sysInfo.kernel_type,'exp')
    rho_exact     = @(u) u.^(-3).*(1-exp(-u*sysInfo.T ));
    rho_exact_val = rho_exact(xgrid); rho_exact_val = rho_exact_val/(sum(rho_exact_val)*dx);
    hold on;
    plot(xgrid, rho_exact_val,'-.','linewidth',1);
end
B          = diag(rho);

%% analysis function space of identifability 
method = 'svdA'; % 'svdA' 'svdAB': should use svdA, which uses eig(A,B), because otherwise, the G-eig does not satisify AV= BVS, V'BV=I.
[V_A,eigA,V_AB, eigAB,r]= EigenAB_fsoi(A,B,1,method,exp_poly); 
% % % compuareGSVD_Geig(A,B);   % compare gsvd and eig(A,B) for G-eig. Should use eig(A,B), not svd(A/B)

%% unconstained LSE with regularuzations: l2, L2, RKHS
% includes a single test for demonstration and tuning and multiple tests for robustness
reguLSE_unconstrained; 

return; 

%% y-mesh tests: it is a different topic: the model changes with mesh. 
y_mesh_tests; 

%% additional explorations: weighted loss; constraied inversion. NOT included in paper. 

%% test weighted loss: not good.
%{
%}
script_weighted_loss; 

%% constrained inversion: 
%{
===  tested 5 methods below:
  - LSE+normalization;    >> best, why???
  - augment to transform (strong/weak constraints)
  - augment to A,b (strong/weak constraints)
=== tested fmincon(f,lambda): not good for RKHS because of the singular C_rkhs
%}
script_constrained_inversion; 

