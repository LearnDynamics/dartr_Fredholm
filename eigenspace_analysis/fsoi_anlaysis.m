function fsoi_anlaysis(B,V_A, eigA,V_L,eigL,r,f_true,f_est)
% analyse the project of f_true into the FSOI
xn  = length(f_true); 

figure; s= subplot(211); 
yyaxis left;
f_proj_eigenL = V_L'*B*f_true;
f_est_proj    = V_L'*B*f_est;
plot(1:xn,f_proj_eigenL,1:xn,f_est_proj,'linewidth',1); 
%plot(1:xn,f_proj_eigenL,'-.','linewidth',2); hold on
% plot(1:xn,f_est_proj,'--','linewidth',1); 
% symlog(s, 'y');

% super-impose eigenvalues
yyaxis right;  
semilogy(1:xn,abs(eigL),'k:','linewidth',2);
yticklabels('auto'); % h=gca;  set(h,'YColor',[0.7 0.7 0.7])
legend('f-true projection','Estimator projection','Eigenvalue');
title('Projection EigenAB'); 

subplot(212);  
yyaxis left;
f_proj_eigenA = V_A'*f_true;
f_est_proj    = V_A'*f_est;
plot(1:xn,f_proj_eigenA,1:xn,f_est_proj,'linewidth',1); 
% plot(1:xn,f_proj_eigenA,'-.','linewidth',2); hold on
% plot(1:xn,f_est_proj,'--','linewidth',1); 

yyaxis right;  
semilogy(1:xn,abs(eigA),'k:','linewidth',2);
yticklabels('auto'); % h=gca;  set(h,'YColor',[0.7 0.7 0.7])
legend('f-true projection','Estimator projection','Eigenvalue');
title('Projection EigenA'); 

end


function semilogy_neg(x,y,label_x,label_y)  % incorrect
% semilogy_neg(1:xn,f_proj_eigenL,'x'); 
% plot y in log-scale, allowing y to be either positive or negative 
ylog = sign(y).*log10(abs(y));
plot(x,ylog,'linewidth',1); 
if ~exist('label_x','var'); label_x = 'x'; end
if ~exist('label_y','var'); label_y = 'sign(y)log10(y)'; end
xlabel(label_x);ylabel(label_y);
end