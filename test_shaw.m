% Using the example shaw for testing. 
% Adaptive RKHS regularization for Discrete Fredhold integral equation
%{
 Dicrete-Fredholm integral equation:
      \int_lb^rb K(t,x) f(x)dx + noise = y(t)
 Goal: Given y(t) at discrete-times, to estimate f
%}
%{
In general: weighted Deconvolution in the form of inversion
               L f = y,   size(L) =  n_y x n_u       size(f) = n_u 
Solution: least square with RKHS-regularization 
Key words: space of identifiability, exploration measure, RKHS regularization 
@Copyright: Fei Lu, feilu@math.jhu.edu. 2023/10/2-2024/12/20

-- TO Further update: use SVD as input for regularizers to avoid repetitive computation.  
%}



clc; close all; clear all;
add_mypaths_discrete;                    % get SAVE_DIR = local dir for saving data
rng(1)
%% Load system settings

exp_poly = 'shaw';
n=200; 
[L_operator,b,f_true] = shaw(n); 
dx = pi/n;
xgrid = -pi/2 + (.5:n-.5)*dx;  tgrid = xgrid; 

xn    = n; tn    = n; 
sysInfo.tn    = n; 
sysInfo.tgrid = xgrid; 
sysInfo.xgrid = xgrid; 
sysInfo.dx    = dx;
sysInfo.dt    = dx; sysInfo.T = pi; 
sysInfo.L_operator = L_operator; 

%% Get regression matrix A, 
% to get vector b later for different f
A = L_operator'*L_operator;  

%% Get rho and L2(rho) basis matrix B 
rho = sum(L_operator);  rho = rho/(sum(rho)*dx);  % normalize: not helpful if rho~unif; like pre-conditioning in ill-posed settings. 
figure; % plot the exploration measure 
plot(xgrid, rho,'linewidth',1); xlabel('u');ylabel('rho');
B          = diag(rho);

%% analysis function space of identifability 
method = 'svdA'; % 'svdA' 'svdAB': should use svdA, which uses eig(A,B), because otherwise, the G-eig does not satisify AV= BVS, V'BV=I.
[V_A,eigA,V_AB, eigAB,r]= EigenAB_fsoi(A,B,1,method,exp_poly); 


%% unconstained LSE with regularuzations: l2, L2, RKHS
% includes a single test for demonstration and tuning and multiple tests for robustness
nsr_seq  = [0.125,0.25,0.5,1,2];     % noise to signal ratio    
           % -- issue when nsr=0: the optimal lambda=0, but numerical error in inversion prevents us from get to it.  
           % Solution: add lambda =0 estimator, and select between lambda_opt and 0 by min-loss ( a factor (1e2) to be robust) >> estimator
normType  = {'l2','L2','RKHS'}; 

%% 1. f_true outside the FSOI: 



 % f_true = 0.1*V_AB(:,2) +2*V_AB(:,30) ; 
% figure; plot(xgrid,f_true/(dx*sum(f_true))); 
 
file_str  =['outsideFSOI_',method,exp_poly];   % 'outsideFSOI'; % %'outsideFSOI_Gaussian_mix'; 
data_name = [SAVE_DIR,'/data_',file_str,'.mat']; 
fig_dir  = [SAVE_DIR,'figures/']; if ~exist(fig_dir,'dir'), mkdir(fig_dir); end  

% % single simulultion demo and tuning. 
 single_simul_demo; 


% multiple simulations 
fprintf('Multiple tests: f true outside the FSOI\n ');
[err_L2rho_projAB,err_cells,err_l2_projA,loss_array] = multi_simuls(f_true,A,B,V_AB,V_A,r,xgrid,rho,eigAB,sysInfo,nsr_seq,normType,data_name);

newfigure =1; xlabelstr = 'nsr'; 
% L2 norm and L2(rho) norm

figure; 
subplot(121); 
label_y = 'L^2(\rho) error (Log10)'; 
plot_mean_std(log10(nsr_seq),log10(err_L2rho_projAB),0,normType,label_y,[],xlabelstr);  
subplot(122)
label_y = 'Loss value (Log10)';
plot_mean_std(log10(nsr_seq),log10(loss_array),0,normType,label_y,[],xlabelstr); 

figname = [fig_dir,file_str,'_L2rho_Loss2']; 
set_positionFontsAll;   print([figname,'.pdf'],'-dpdf', '-bestfit'); 



%% 2. f_true inside the FSOI: 
f_true = V_AB(:,1:5)*(1:5)';           case_num= '2';  % numerically: not in RKHS, not in FSOI since eig5 is large
 % f_true = V_AB(:,5);    % ind>5 ( i.e, eigAB<1e-9): rkhs not good, l2+L2 can slightly tolerate more, bc. not using rkhs inversion).
 % f_true = f_true/sqrt(dx*sum(f_true)); 
f_true = V_AB(:,2);                    case_num= '';   % numerically: in RKHS, in FSOI; no decay coefs  >> used in paper
% f_true = V_AB(:,1:5)*sqrt(eigAB(1:5)); case_num= '3';  % numerically: in RKHS, in FSOI; with decaying coefs >>> sharp rate


file_str  = ['insideFSOI_',case_num,'_',method,exp_poly]; 
data_name = [SAVE_DIR,'/data_',file_str,'.mat']; 
fig_dir  = [SAVE_DIR,'figures/']; if ~exist(fig_dir,'dir'), mkdir(fig_dir); end  

 % single simulultion demo 
single_simul_demo; 

% multiple simulations 
fprintf('Multiple tests: f true inside the FSOI\n ');
[err_L2rho_projAB2,err_cells,err_l2_projA2,loss_array2] = multi_simuls(f_true,A,B,V_AB,V_A,r,xgrid,rho,eigAB,sysInfo,nsr_seq,normType,data_name); 

% L2 norm and L2(rho) norm 
figure; 
subplot(121); 
label_y = 'L^2(\rho) error (Log10)'; 
plot_mean_std(log10(nsr_seq),log10(err_L2rho_projAB2),0,normType,label_y,[],xlabelstr);  
subplot(122)
label_y = 'Loss value (Log10)';
plot_mean_std(log10(nsr_seq),log10(loss_array2),0,normType,label_y,[],xlabelstr); 

figname = [fig_dir,file_str,'_L2rho_Loss2']; 
set_positionFontsAll;   print([figname,'.pdf'],'-dpdf', '-bestfit'); 

%% not used parts

function boxplot_err(err_L2rho_outside,nsr_seq)  % not good for view
 % boxplot of the errors  : not good 
figure; 
err = squeeze(err_L2rho_outside(:,1,:));
subplot(131); boxplot(err,nsr_seq); hold on;  xlabel('nsr'); ylabel('L2(rho) error'); title('L2rho errors-l2');
err = squeeze(err_L2rho_outside(:,2,:));
subplot(132); boxplot(err,nsr_seq); hold on;  xlabel('nsr'); ylabel('L2(rho) error'); title('L2rho errors: L2(rho)');
err = squeeze(err_L2rho_outside(:,3,:));
subplot(133); boxplot(err,nsr_seq); hold on;  xlabel('nsr'); ylabel('L2(rho) error'); title('L2rho errors: RKHS');
set_positionFontsAll;
end
