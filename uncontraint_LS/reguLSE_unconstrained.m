%% 

nsr_seq  = [0.125,0.25,0.5,1,2];     % noise to signal ratio    
           % -- issue when nsr=0: the optimal lambda=0, but numerical error in inversion prevents us from get to it.  
           % Solution: add lambda =0 estimator, and select between lambda_opt and 0 by min-loss ( a factor (1e2) to be robust) >> estimator
normType  = {'l2','L2','RKHS'}; 

%% 1. f_true outside the FSOI: 
  % f_true from functions 
% f_true_func = @(x) (sin(x-6)).^2+1;   % True f
% f_true_func = @(x) 0.7*exp(-(x-2).^2/0.25)*sqrt(1/(2*.5*pi)) +.3*exp(-(x-4).^2/0.09)*sqrt(1/(2*0.3*pi));   % True f
 f_true_func = @(x) 15*((sin(x-6)).^2-3);   % True f
% f_true_func = @(x) 15*((sin(2*x-6)).^2-3);
 f_true_func = @(x) x.^2;   %    outside FSOI: components nonzero, but has decaying coefs inside FSOI
 f_true = f_true_func(xgrid);

 % f_true = 0.1*V_AB(:,2) +2*V_AB(:,30) ; 
%  figure; plot(xgrid,f_true/(dx*sum(f_true))); 
 
file_str  =['outsideFSOI_',method];   % 'outsideFSOI'; % %'outsideFSOI_Gaussian_mix'; 
data_name = [SAVE_DIR,'/data_',file_str,'.mat']; 
fig_dir  = [SAVE_DIR,'figures/']; if ~exist(fig_dir,'dir'), mkdir(fig_dir); end  

% % single simulultion demo and tuning. 
 single_simul_demo; 


% multiple simulations 
fprintf('Multiple tests: f true outside the FSOI\n ');
[err_L2rho_projAB,err_cells,err_l2_projA,loss_array] = multi_simuls(f_true,A,B,V_AB,V_A,r,xgrid,rho,eigAB,sysInfo,nsr_seq,normType,data_name);

newfigure =1; xlabelstr = 'nsr'; 
%% L2 norm and L2(rho) norm
% label_y = ' L^2(\rho) error  (Log10)'; 
% fig_name = [fig_dir,file_str,'_L2rho_error']; 
% plot_mean_std(log10(nsr_seq),log10(err_L2rho_projAB),newfigure,normType,label_y,fig_name,xlabelstr);  
% 
% 
% label_y = 'l2 error  (Log10)';
% fig_name = [fig_dir,file_str,'_l2error']; 
% plot_mean_std(log10(nsr_seq),log10(err_l2_projA),newfigure,normType,label_y,fig_name,xlabelstr);  
% 
% label_y = 'Loss value  (Log10)';
% fig_name = [fig_dir,file_str,'_loss']; 
% plot_mean_std(log10(nsr_seq),log10(loss_array),newfigure,normType,label_y,fig_name,xlabelstr);  

% fig_name = [fig_dir,file_str,'_full_estimator_err_'];
% full_estimator_err(err_cells,normType,fig_name,nsr_seq,xlabelstr);

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


file_str  = ['insideFSOI_',case_num,'_',method]; 
data_name = [SAVE_DIR,'/data_',file_str,'.mat']; 
fig_dir  = [SAVE_DIR,'figures/']; if ~exist(fig_dir,'dir'), mkdir(fig_dir); end  

 % single simulultion demo 
single_simul_demo; 

% multiple simulations 
fprintf('Multiple tests: f true inside the FSOI\n ');
[err_L2rho_projAB2,err_cells,err_l2_projA2,loss_array2] = multi_simuls(f_true,A,B,V_AB,V_A,r,xgrid,rho,eigAB,sysInfo,nsr_seq,normType,data_name); 

%% L2 norm and L2(rho) norm 
% newfigure =1; label_y = 'L^2(\rho) error (Log10)'; 
% fig_name = [fig_dir,file_str,'_L2rho_error']; 
% plot_mean_std(log10(nsr_seq),log10(err_L2rho_projAB2),newfigure,normType,label_y,fig_name,xlabelstr);  
% 
% label_y = 'l2 error  (Log10)';
% fig_name = [fig_dir,file_str,'_l2error']; 
% plot_mean_std(log10(nsr_seq),log10(err_l2_projA2),newfigure,normType,label_y,fig_name,xlabelstr);  
% 
% label_y = 'Loss value  (Log10)';
% fig_name = [fig_dir,file_str,'_loss']; 
% plot_mean_std(log10(nsr_seq),log10(loss_array2),newfigure,normType,label_y,fig_name,xlabelstr);  


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
