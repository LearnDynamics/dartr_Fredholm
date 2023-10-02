function full_estimator_err(err_cells,normType,fig_name,nsr_seq,xlabelstr)
% function [L2rho_err_full_est, l2_err_full_est]= extract_L2rho_err_full_estimator(err_cells,normType)
% extract the full estimator's L2 error 
[n_simuls,n_nsr] = size(err_cells); 

temp1 = zeros(n_simuls,length(normType),n_nsr);
temp2 = temp1; 
for n=1:n_simuls
    for j = 1:n_nsr
        err1 = err_cells{n,j};
        temp1(n,:,j)  = err1.L2rho_l2(1,:);
        temp2(n,:,j)  = err1.L2rho_l2(2,:);
    end
end
L2rho_err_full_est = temp1; 
l2_err_full_est    = temp2; 

fig_name1 = [fig_name,'_L2rho']; 

figure; 
subplot(121)
label_y = 'L2rho error'; 
plot_mean_std(nsr_seq,L2rho_err_full_est,0,normType,label_y,fig_name,xlabelstr);  
title('Full estimator L2rho error'); 

subplot(122)
label_y = 'l2 error';
plot_mean_std(nsr_seq,l2_err_full_est,0,normType,label_y,fig_name,xlabelstr);  
title('Full estimator l2 error'); 


end