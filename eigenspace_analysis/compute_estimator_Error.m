function [errors] = compute_estimator_Error(B,V_AB,eigAB,V_A,r,est_array,rho,xgrid)
% compute the errors of the estimator
%{
  - estimator errors: L2(rho), l2;  (rkhs can be infinity, so not computed)
  - projection errors: L2(rho),l2 , rkhs;   
%}


%% vector estimator: L2rho, l2  (L2 = l2* scalar)
[n_vec,n_est] = size(est_array);  
dx_array      = xgrid(2:end)-xgrid(1:end-1); dx_array = [dx_array;dx_array(end)]; 
true_val      = est_array(:,end);
diff_array    = est_array(:,1:end-1)- true_val;  
err_L2rho     = rho*diff_array.^2;         err_L2rho  = sqrt(err_L2rho); 
err_l2        = sum(diff_array.^2,1);       err_l2    = sqrt(err_l2);

errors.L2rho_l2  = [err_L2rho;err_l2]; 
errors.L2rho_l2_discription = {'L2rho errors: l2,L2,RKHS'; 'l2 errors: l2,L2,RKHS'}; 
%% projection of the estimators: in L2(rho), project to V_AB space 
%  errL2rho^2 = sum_i c_i^2 int_x phi_i(x)^2 rho(x) dx =   sum_i c_i^2, because V_AB basis are onb in L2rho
est_array_proj = V_AB'*B*est_array;  
est_array_proj = est_array_proj(1:r,:); 
true_val_proj  = est_array_proj(:,end);

diff_array     = est_array_proj(:,1:end-1)- true_val_proj; 
f_diff         = V_AB(:,1:r)*diff_array; 
err_L2rho_proj = sum(diff_array.^2,1);       err_L2rho_proj = sqrt(err_L2rho_proj); 
err_l2_proj    = sum(f_diff.^2,1);           err_l2_proj    = sqrt(err_l2_proj);
err_rkhs_proj  = sqrt(sum(diag(1./eigAB(1:r))*diff_array.^2, 1)); 
errors.L2rho_l2_projAB = [err_L2rho_proj;err_l2_proj; err_rkhs_proj];


errors.discription_L2rho_l2_projAB= {'L2rho errors AB-projected: l2,L2,RKHS '; 'l2 errors AB-projected: l2,L2,RKHS';'rkhs errors AB-projected: l2,L2,RKHS'}; 
%
%% projection of the estimators: in L2, project to V_A space
%  --- this is not compatiable with the function space of learning, L2rho. But it is in common-use 
%  errL2rho^2 = sum_i c_i^2 int_x phi_i(x)^2 rho(x) dx =   sum_i c_i^2, because V_A basis are onb in L2
est_array_proj = V_A'*eye(n_vec)*est_array;  
est_array_proj = est_array_proj(1:r,:); 
true_val_proj  = est_array_proj(:,end);

diff_array     = est_array_proj(:,1:end-1)- true_val_proj; 
f_diff         = V_A(:,1:r)*diff_array; 
err_L2rho_proj = rho*f_diff.^2;         err_L2rho_proj  = sqrt(err_L2rho_proj); 
err_l2_proj    = sum(f_diff.^2,1);       err_l2_proj    = sqrt(err_l2_proj);

errors.L2rho_l2_projA = [err_L2rho_proj; err_l2_proj];
errors.discription_L2rho_l2_projA = {'L2rho-A-projected: l2,L2,RKHS';'l2-A-projected: l2,L2,RKHS'}; 
%}

end