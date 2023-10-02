%% weighted loss
weight_t  = 1+exp(-(sysInfo.rb+sysInfo.lb)/2*tgrid);  % no weight seems to be the best; exp(a*t) is not good: it amplifies the noise; 
weight_t = weight_t /sum(weight_t ); 
wt_mat    = diag(weight_t);

A_wt   = L_operator'*wt_mat*L_operator;  
rho_wt = sum(wt_mat*L_operator); rho_wt = rho_wt/(sum(rho_wt)*dx);  % normalize, does not seem necessary for this example, maybe other examples
B_wt   = diag(rho_wt);
figure; % plot the exploration measure 
plot(xgrid, rho_wt,'linewidth',1); xlabel('u');ylabel('rho');

[V_A_wt,eigA_wt,V_AB_wt, eigAB_wt,r_wt]= EigenAB_fsoi(A_wt,B_wt); 


b_wt      = L_operator'*wt_mat*y;
[x_reg,lambda_opt,cLcurveAll] = reg_Lcurve_3in1(A_wt,b_wt,B_wt,plotON,normType); 
lgnd      = {'Est-l2','Est-L2','Est-RKHS','True'}; 
est_array = [cLcurveAll.creg_l2,cLcurveAll.creg_L2,cLcurveAll.creg_RKHS,f_true];
plot_estimators(est_array,lgnd,xgrid,rho_wt); 
plot_estimators_image(est_array,tgrid,L_operator,y,y_true,lgnd); % comment: we can use it for image denoising. 
