function eigenspace_projection(B,V,eig_val,r,f_true,c,L_operator,y,y_true,titl)
% plot the eigenspace projection of: 
%       ftrue, estimator (regularized), truncated SVD for y and y_true
% Input: 
%     B          - the basis matrix 
%    V,eig_val   - eigenvector, eigenvalue 
%     r          - number of eigenvectors for projection (eig_val>1e-8)

coord_true = V'*B*f_true;
coord_est  = V'*B*c;
coord_naiv = pinv(diag(eig_val))*V'*L_operator'*y;
coord_naiv_true = pinv(diag(eig_val))*V'*L_operator'*y_true;


coord_true = coord_true(1:r);
coord_est = coord_est(1:r);
coord_naiv = coord_naiv(1:r);
coord_naiv_true = coord_naiv_true(1:r);


figure; subplot(121)
plot(coord_true,'linewidth',1);hold on;
plot(coord_est,'linewidth',1);
plot(coord_naiv, '-x', 'MarkerSize', 5)
plot(coord_naiv_true, '-o', 'MarkerSize', 5)
legend('coordinate of f true', 'f est', 'using y noisy, direct cut-off', 'using y true, r cut-off')
legend('location','best');
ylim([min(coord_true)-3, max(coord_true)+3]);
title([titl,'r-coeffs'])

project_true = V(:, 1:r)*coord_true;
project_est = V(:, 1:r)*coord_est;
project_naiv = V(:, 1:r)*coord_naiv;
project_naiv_true = V(:, 1:r)*coord_naiv_true;

% figure;
subplot(122)
plot(project_true,':','linewidth',3);hold on;
plot(project_est,'-','linewidth',1);
plot(f_true, 'k:', 'linewidth',3);
plot(c, '--', 'linewidth',2)
plot(project_naiv, '-.x','linewidth',1);
plot(project_naiv_true, '-.o','linewidth',1);
ylim([min(project_true)-3, max(project_true)+3]);
legend('f true project', 'f est project', 'f true', 'f est', 'using y noisy, r cut-off', 'using y true, r cut-off')
legend('location','best');
title([titl,'r-coef. estimator'])

end