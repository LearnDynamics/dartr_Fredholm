% clc; close all; clear all;
add_mypaths_discrete;    % get SAVE_DIR = local dir for saving data
rng(1)

%% different y-mesh grids

gap_seq  =10*[1,2,4,8,16];     % different y-mesh to demonstrate convergence   

normType  = {'l2','L2','RKHS'}; 

% % in discrete: sysInfo.tn = 5000; for [0,5]
sysInfo     = system_settings();
sysInfo.tn  = 10000;   % 5000
sysInfo     = update_system_settings(sysInfo); 
L_operator_fine = sysInfo.L_operator; 
xn    = sysInfo.xn;
tn    = sysInfo.tn;
tgrid = sysInfo.tgrid;
xgrid = sysInfo.xgrid;
dx    = sysInfo.dx;
ymesh_seq  = sysInfo.dt*gap_seq; 

exp_poly = sysInfo.kernel_type; 

%% Get: A_conti rho and L2(rho) basis matrix B 
A_continuous = L_operator_fine'*L_operator_fine;  

% rho = sum(L_operator_fine);  rho = rho/(sum(rho)*dx);  % normalize, does not seem necessary for this example, maybe other examples
 figure; % plot the exploration measure 
% plot(xgrid, rho,'linewidth',1); xlabel('u');ylabel('rho'); hold on;  
rho_exact     = @(u) u.^(-3).*(1-exp(-u*sysInfo.T ));  
rho_exact_val = rho_exact(xgrid); rho_exact_val = rho_exact_val'/(sum(rho_exact_val)*dx);
plot(xgrid, rho_exact_val,'-.','linewidth',1); 
% max(rho-rho_exact_val);
rho = rho_exact_val; 
B   = diag(rho);


%% analysis function space of identifability 
method = 'svdA'; % 'svdA' 'svdAB': should use svdA, which uses eig(A,B), because otherwise, the G-eig does not satisify AV= BVS, V'BV=I.
[V_A,eigA,V_AB, eigAB,r]= EigenAB_fsoi(A_continuous,B,1,method,exp_poly); 


%% 1. f_true outside the FSOI: 
  % f_true from functions 
% f_true_func = @(x) (sin(x-6)).^2+1;   % True f
% f_true_func = @(x) 0.7*exp(-(x-2).^2/0.25)*sqrt(1/(2*.5*pi)) +.3*exp(-(x-4).^2/0.09)*sqrt(1/(2*0.3*pi));   % True f
% f_true_func = @(x) 15*((sin(x-6)).^2-3);   % True f
% f_true_func = @(x) 15*((sin(2*x-6)).^2-3);
  f_true_func = @(x) x.^2;   % True f
f_true = f_true_func(xgrid);
% f_true = 0.1*V_AB(:,2) +2*V_AB(:,30) ; 
% figure; plot(xgrid,f_true/(dx*sum(f_true))); 
 
file_str  = ['ymesh_outsideFSOI_',method,exp_poly]; % 'outsideFSOI_Gaussian_mix'; 
data_name = [SAVE_DIR,'/data_',file_str,'.mat']; 
fig_dir   = [SAVE_DIR,'/figures/']; if ~exist(fig_dir,'dir'), mkdir(fig_dir); end  

%  single simulultion demo and tuning. 
single_simul_demo_ymesh; 

nsr = 1; 
% multiple simulations 
fprintf('Multiple tests: f true outside the FSOI\n ');
[err_L2rho_projAB,err_cells,err_l2_projA,loss_array,disc_est] = multi_simuls_ymesh(f_true,A,B,V_AB,V_A,r,xgrid,rho,eigAB,sysInfo,gap_seq,nsr,normType,data_name);

 newfigure =0; xlabelstr = '\Delta t'; 
% label_y = 'L^2(\rho) error'; 
% fig_name = [fig_dir,file_str,'_L2rho_error']; 
% plot_mean_std(log10(ymesh_seq),err_L2rho_projAB,newfigure,normType,label_y,fig_name,xlabelstr);  
% 
% label_y = 'l2 error';
% fig_name = [fig_dir,file_str,'_l2error']; 
% plot_mean_std(log10(ymesh_seq),err_l2_projA,newfigure,normType,label_y,fig_name,xlabelstr);  
% 
% label_y = 'Loss value';
% fig_name = [fig_dir,file_str,'_loss']; 
% plot_mean_std(log10(ymesh_seq),log10(loss_array),newfigure,normType,label_y,fig_name,xlabelstr);  

% % fig_name = [fig_dir,file_str,'_full_estimator_err_'];
% % full_estimator_err(err_cells,normType,fig_name,ymesh_seq,xlabelstr);

% % discA estimator
% label_y = 'L^2(\rho) error (Log10)'; 
% fig_name = [fig_dir,file_str,'_L2rho_error_discA']; 
% plot_mean_std(log10(ymesh_seq),log10(disc_est.err_L2rho_projAB),newfigure,normType,label_y,fig_name,xlabelstr);  
% 
% label_y = 'l2 error (Log10) ';
% fig_name = [fig_dir,file_str,'_l2error_discA']; 
% plot_mean_std(log10(ymesh_seq),log10(disc_est.err_l2_projA),newfigure+1,normType,label_y,fig_name,xlabelstr);  
% 
% label_y = 'Loss value (Log10)';
% fig_name = [fig_dir,file_str,'_loss_discA']; 
% plot_mean_std(log10(ymesh_seq),log10(disc_est.loss_array),newfigure,normType,label_y,fig_name,xlabelstr); 


figure; 
subplot(121); 
label_y = 'L^2(\rho) error (Log10)'; 
plot_mean_std(log10(ymesh_seq),log10(disc_est.err_L2rho_projAB),0,normType,label_y,[],xlabelstr);  
subplot(122)
label_y = 'Loss value (Log10)';
plot_mean_std(log10(ymesh_seq),log10(disc_est.loss_array),0,normType,label_y,[],xlabelstr); 

figname = [fig_dir,file_str,'_L2rho_Loss_discA2']; 
set_positionFontsAll;   print([figname,'.pdf'],'-dpdf', '-bestfit'); 



%% 2. f_true inside the FSOI: 
f_true = V_AB(:,1:5)*(1:5)';           case_num= '2';  % numerically: not in RKHS, not in FSOI since eig5 is large
 % f_true = V_AB(:,5);    % ind>5 ( i.e, eigAB<1e-9): rkhs not good, l2+L2 can slightly tolerate more, bc. not using rkhs inversion).
 % f_true = f_true/sqrt(dx*sum(f_true)); 
f_true = V_AB(:,2);                    case_num= '';   % numerically: in RKHS, in FSOI; no decay coefs as theorem 3
f_true = V_AB(:,1:5)*sqrt(eigAB(1:5)); case_num= '3';  % numerically: in RKHS, in FSOI; with decaying coefs >>> sharp rate


file_str  = ['ymesh_insideFSOI_',case_num,'_',method,exp_poly]; 
data_name = [SAVE_DIR,'/data_',file_str,'.mat']; 
fig_dir  = [SAVE_DIR,'/figures/']; if ~exist(fig_dir,'dir'), mkdir(fig_dir); end  

% single simulultion demo 
single_simul_demo_ymesh;

% multiple simulations 
 fprintf('Multiple tests: f true inside the FSOI\n ');
[err_L2rho_projAB,err_cells,err_l2_projA,loss_array,disc_est] = multi_simuls_ymesh(f_true,A,B,V_AB,V_A,r,xgrid,rho,eigAB,sysInfo,gap_seq,nsr,normType,data_name);

newfigure =1;
xlabelstr = '\Delta t'; 
% label_y = 'L^2(\rho) error'; 
% fig_name = [fig_dir,file_str,'_L2rho_error']; 
% plot_mean_std(log10(ymesh_seq),log10(err_L2rho_projAB),newfigure,normType,label_y,fig_name,xlabelstr);  
% 
% label_y = 'l2 error';
% fig_name = [fig_dir,file_str,'_l2error']; 
% plot_mean_std(log10(ymesh_seq),log10(err_l2_projA),newfigure,normType,label_y,fig_name,xlabelstr);  
% 
% label_y = 'Loss value';
% fig_name = [fig_dir,file_str,'_loss']; 
% plot_mean_std(log10(ymesh_seq),log10(loss_array),newfigure,normType,label_y,fig_name,xlabelstr);  


% % discA estimator
% label_y = 'L^2(\rho) error (Log10)'; 
% fig_name = [fig_dir,file_str,'_L2rho_error_discA']; 
% plot_mean_std(log10(ymesh_seq),log10(disc_est.err_L2rho_projAB),newfigure,normType,label_y,fig_name,xlabelstr);  
% 
% label_y = 'l2 error (Log10)';
% fig_name = [fig_dir,file_str,'_l2error_discA']; 
% plot_mean_std(log10(ymesh_seq),log10(disc_est.err_l2_projA),newfigure,normType,label_y,fig_name,xlabelstr);  
% 
% label_y = 'Loss value (Log10)';
% fig_name = [fig_dir,file_str,'_loss_discA']; 
% plot_mean_std(log10(ymesh_seq),log10(disc_est.loss_array),newfigure,normType,label_y,fig_name,xlabelstr); 


figure; 
subplot(121); 
label_y = 'L^2(\rho) error (Log10)'; 
plot_mean_std(log10(ymesh_seq),log10(disc_est.err_L2rho_projAB),0,normType,label_y,[],xlabelstr);  
subplot(122)
label_y = 'Loss value (Log10)';
plot_mean_std(log10(ymesh_seq),log10(disc_est.loss_array),0,normType,label_y,[],xlabelstr); 

figname = [fig_dir,file_str,'_L2rho_Loss_discA2']; 
set_positionFontsAll;   print([figname,'.pdf'],'-dpdf', '-bestfit'); 

