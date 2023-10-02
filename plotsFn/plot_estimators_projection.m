function plot_estimators_projection(B,V,r,est_array,lgnd,titl,new_figure)
% plot the eigenspace projection of: 
%       ftrue, estimator (regularized)
% Input: 
%     B          - the basis matrix 
%    V,eig_val   - eigenvector, eigenvalue 
%     r          - number of eigenvectors for projection (eig_val>1e-8)
%     Est_array: each colum is an estimator, with the last column = true


colors; % get linestyle and dark colors

[n_vec,n_est] = size(est_array);  
xi = 1:n_vec; 
if ~exist('new_figure','var'); new_figure=1; end
if new_figure==1; figure; end 

%% projection of the estimators 
est_array_proj = V'*B*est_array;  
true_val       = est_array_proj(:,end);


for n=1:n_est-1
    plot(xi(1:r),est_array_proj(1:r,n),linestyle{n},'Color',dred_do_db(n,:),'linewidth',1); hold on;
end
plot(xi(1:r),true_val(1:r),'linestyle',linestyle{n_est},'Color','k','linewidth',2); 

 
if min(true_val) >0 ;  lb = min(true_val)*0.8; 
else;                  lb = min(true_val)*1.2; 
end
if  max(true_val)>0;   ub = max(true_val)*1.2; 
else;                  ub = max(true_val)*0.8; 
end 
ylim([lb,ub]); ylabel('Coefficients'); xlabel('Coefficient index'); 
title(titl);

legend(lgnd); legend('location','best')



end