%% test constrained inversion: 
%{
===  tested 5 methods below:
  - LSE+normalization;    >> best, why???
  - augment to transform (strong/weak constraints)
  - augment to A,b (strong/weak constraints)
=== tested fmincon(f,lambda): not good for RKHS because of the singular C_rkhs
%}
f_true2 = f_true/(dx*sum(f_true)); 
y_true2 = L_operator*f_true2; 
y_true_norm = sqrt(sum(y_true2.^2)*sysInfo.dt/sysInfo.T);
y2 = y; %       = y_true2 + y_true_norm*nsr*randn(tn, 1);
b2 = b; %     = L_operator'*y2; 

%  estimate by regression and normalize at the end. NOT "correct": Ax=b not satisfied 
[x_reg,lambda_opt,cLcurveAll] = reg_Lcurve_3in1(A,b2,B,plotON,normType);
est_array = [cLcurveAll.creg_l2,cLcurveAll.creg_L2,cLcurveAll.creg_RKHS,f_true2];
est_array = est_array*diag(1./(dx*sum(est_array))); 
plot_estimators(est_array,lgnd,xgrid,rho); 
plot_estimators_image(est_array,tgrid,L_operator,y2,y_true2,lgnd);  

% augment to the transform, weak constraint, so that the function space is right >>> need to change rho?
con_type ='strong';
[A_con_bar,b_con_bar,rho_con] = LS_Ab_con_sum1(L_operator,y2,rho,con_type,dx);
B_con   = diag(rho_con);
[x_reg,lambda_opt,cLcurveAll] = reg_Lcurve_3in1(A_con_bar,b_con_bar,B_con,plotON,normType);  
estl2 = cLcurveAll.creg_l2;      estl2  = [estl2; 1/dx-sum(estl2)]; 
estLL2 = cLcurveAll.creg_L2;    estLL2  = [estLL2; 1/dx-sum(estLL2)]; 
estRKHS = cLcurveAll.creg_RKHS; estRKHS = [estRKHS; 1/dx-sum(estRKHS)]; 
est_array = [estl2,estLL2,estRKHS,f_true2];
plot_estimators(est_array,lgnd,xgrid,rho); 
plot_estimators_image(est_array,tgrid,L_operator,y2,y_true2,lgnd);

% augment to the transform, weak constraint, so that the function space is right >>> need to change rho? 
con_type ='weak';
[A_con_bar,b_con_bar,rho_con] = LS_Ab_con_sum1(L_operator,y2,rho,con_type,dx);
[x_reg,lambda_opt,cLcurveAll] = reg_Lcurve_3in1(A_con_bar,b_con_bar,B,plotON,normType);  
est_array = [cLcurveAll.creg_l2,cLcurveAll.creg_L2,cLcurveAll.creg_RKHS,f_true2];
est_array = est_array*diag(1./(dx*sum(est_array))); 
plot_estimators(est_array,lgnd,xgrid,rho); 
plot_estimators_image(est_array,tgrid,L_operator,y2,y_true2,lgnd);

% LSE with sum_j x_j*dx =1 added another row to A: weak constraint. Regu-norm should be changed? 
A_ext = [A; ones(1,length(b2))];  b_ext = [b2;dx];  c_ext = A_ext\b_ext; 
A_ext_bar = A_ext'*A_ext; b_ext_bar = A_ext'*b_ext; 
[x_reg,lambda_opt,cLcurveAll] = reg_Lcurve_3in1(A_ext_bar,b_ext_bar,B,plotON,normType);  
est_array = [cLcurveAll.creg_l2,cLcurveAll.creg_L2,cLcurveAll.creg_RKHS,f_true2];
% est_array = est_array*diag(1./(dx*sum(est_array))); 
plot_estimators(est_array,lgnd,xgrid,rho); 
plot_estimators_image(est_array,tgrid,L_operator,y2,y_true2,lgnd);


% LSE with sum_j x_j*dx =1  enforced to (n-1) vector: strong constraint to A,b 
n     = length(b2); 
A_con = zeros(n,n-1);  b_con = b2- A(:,n);  % change the regression matrix and vector: rectangular now!
for i=1:n
   A_con(i,:) =  A(i,1:n-1)- A(i,n); 
end
A_con_bar = A_con'*A_con;    % change to square matrix  ~~~~~ function space may mismatch here! 
b_con_bar = A_con'*b_con; 

rho_con = rho(1:n-1)+rho(n); rho_con = rho_con/(sum(rho_con)); 
B_con   = diag(rho_con);
[x_reg,lambda_opt,cLcurveAll] = reg_Lcurve_3in1(A_con_bar,b_con_bar,B_con,plotON,normType);
estl2 = cLcurveAll.creg_l2;      estl2  = [estl2; 1-sum(estl2)]; 
estLL2 = cLcurveAll.creg_L2;    estLL2  = [estLL2; 1-sum(estLL2)]; 
estRKHS = cLcurveAll.creg_RKHS; estRKHS = [estRKHS; 1-sum(estRKHS)]; 
est_array = [estl2,estLL2,estRKHS,f_true2];
plot_estimators(est_array,lgnd,xgrid,rho); 
plot_estimators_image(est_array,tgrid,L_operator,y2,y_true2,lgnd);


% using fmincon: not good 
normtype= 'RKHS'; % 'l2','L2','RKHS'
[f_reg,lambda_opt,loss_reg] = lse_by_opt_con(L_operator,y,B,dx, normtype);
est_array= [f_reg,f_true2]; lgnd= {normtype,'f-true'}; 
plot_estimators(est_array,lgnd,xgrid,rho); 

