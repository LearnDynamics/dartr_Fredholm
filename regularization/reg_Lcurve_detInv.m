function [x_reg,lambda_opt,cLcurveAll] = reg_Lcurve_detInv(A,b,B,plotON, normType, phi)
% regularizaiton using L-curve to select optimal lambda, considering three norms: l2, L2, L2-RKHS
% NOTE: L-curve RKHS may need a different range for the parameter.   
%{
normType  = {'l2','L2','RKHS'} or part of them; other types are possible. 
                  x_l = argmin_x  |Ax -b|^2 + lambda x'*B*x
                      = (A'*A+lambda*B)\(A'*b)
 Set lambda to damp the error in A\b while keeping the norm x'*B*x small
 Select lambda that at the "corner" of the L-curve (maximum curvature)
               ( |Ax_l -b|^2, x_l'*B*x_l)   = (xlambda,ylambda)
 See Hansen: the L-curve and its use in the numerical treatment of inverse problems
Input:  
   - A, b, B= basis matrix  
   - phi: - basis function vectors for H1 norm 
Regulariztaion methods:
    norm_type = 'l2':     |Ax -b|^2 + lambda x'* Id *x
    norm_type = 'L2':     |Ax -b|^2 + lambda x'* B *x,   with B being the L2 norm matrix
    norm_type = 'RKHS':   |Ax -b|^2 + lambda x'* B-rkhs *x,
Note: 
   1. A, b are NOT Abar = A' * A; bbar = A' * b; 
   2. The RKHS norm regularizatin is
      It is equivalent to (Abar*x + lambda*I)\ Abar*b
% Copyright(c): Fei Lu, 2022-09-26
%}

% tol = 1e-20;
%  B   = Lift_smallestEigen(B,tol);  % replace the small eigval by tol; to avoid numerical issues in the generalized eigenvalue problem

[m, n] = length(A); 
%% Loop for all normTypes 
for nn = 1:length(normType)
    norm_type = normType{nn};
    switch norm_type
        case 'l2'     %% 1. Tikhonov regulariztion, norm_type = 'l2';
            Bmat      = eye(length(b));
            titl      = ['L-curve with norm: ',norm_type];
            [x_reg,lambda_opt,loss_val] = Lcurve_with_Norm_pinv(A,b,Bmat,titl,plotON,norm_type);
            cLcurveAll.creg_l2             = x_reg;
            cLcurveAll.creg_l2_lambda_opt  = lambda_opt;
            cLcurveAll.creg_l2_loss_val    = loss_val; 
            %             cLcurveAll.inds_lambda_l2      = find(eigA>lambda_opt);
        case 'L2'     %% 2. Tikhonov regulariztion, norm_type = 'L2';
            Bmat      = B;
            titl      = ['L-curve with norm: ',norm_type];
            [x_reg,lambda_opt,loss_val] = Lcurve_with_Norm_pinv(A,b,Bmat,titl,plotON,norm_type);
            cLcurveAll.creg_L2             = x_reg;
            cLcurveAll.creg_L2_lambda_opt  = lambda_opt;
            cLcurveAll.creg_L2_loss_val    = loss_val; 
            %             cLcurveAll.inds_lambda_L2      = find(eigL>lambda_opt);
        case 'RKHS'   %% 3. Lcurve-rkhs regulariztion, norm_type = 'rkhs'
            titl      = ['L-curve with norm: ',norm_type]; 
            if 0  % to do later: reduce numerical error in B_rkhs
                [V, eigL]  = eig(A,B);  % generalized eigenvalue   A = B*V*eigL; V'*B*V =I;  V'*A*V = eigL; >>>  B_rkhs = inv(V*diag(eigL)*V')
                eigL       = real(diag(eigL));  [~,ind] = sort(eigL,'descend');
                eigL       = eigL(ind); V = V(:,ind);
                % B_rkhs     = pinv(V*diag(real(eigL))*V');   % = pinv(Abar);  B_rkhs= inv(Abar) when Bmat=Id;
                eig_diag_temp = pinv(diag(eigL.^(1/2)));
                temp2 = B*V*eig_diag_temp;
                B_rkhs = temp2*temp2'; % This is the correct way to preserve the symmetry of C
                [x_reg,lambda_opt,loss_val] = Lcurve_with_Norm_pinv(A,b,B_rkhs,titl,plotON,norm_type);
            else   % use lsq to sovlve x_lambda
                [x_reg,lambda_opt,loss_val] = Lcurve_sidaRKHS_lsq2(A,b,B,titl,plotON);
            end            
            cLcurveAll.creg_RKHS             = x_reg;
            cLcurveAll.creg_RKHS_lambda_opt  = lambda_opt;
            cLcurveAll.creg_RKHS_loss_val    = loss_val; 
        case 'H1'
            %%
            if ~exist('phi','var')
                % this L does not have full rank.
                % will cause small issue
                % L = [-eye(n-1, n-1), zeros(n-1, 1)] + [zeros(n-1, 1), eye(n-1, n-1)];
                L = eye(n) + [zeros(1, n-1), 0;-eye(n-1), zeros(n-1, 1)];
                BH1 = L'*L;
            else
                BH1 = zeros(n, n);
                for i = 1:n
                    for j = 1:n
                        BH1(i,j) = gradient(phi(:, i))'*gradient(phi(:, j));
                    end
                end
            end
            titl      = ['L-curve with norm: ',norm_type];
            [x_reg,lambda_opt,loss_val] = Lcurve_with_Norm_pinv(A,b,BH1,titl,plotON,norm_type);
            cLcurveAll.creg_H1             = x_reg;
            cLcurveAll.creg_H1_lambda_opt  = lambda_opt;
            
        case 'H2'
            %%
            if ~exist('phi','var')
                % this L doest not have full rank.
                % will cause small issue.
                %                 L = [-eye(n-1, n-1), zeros(n-1, 1)] + [zeros(n-1, 1), eye(n-1, n-1)];
                %                 LL = [-eye(n-2, n-2), zeros(n-2, 1)] + [zeros(n-2, 1), eye(n-2, n-2)];
                %                 B = L'*LL'*LL*L;
                L = 2*eye(n) + [zeros(1, n-1), 0;eye(n-1), zeros(n-1, 1)] +  [zeros(1, n-1), 0;-eye(n-1), zeros(n-1, 1)]';
                BH2 = L' * L;
            else
                BH2 = zeros(n, n);
                for i = 1:n
                    for j = 1:n
                        BH2(i,j) = gradient(gradient(phi(:, i)))'*gradient(gradient(phi(:, j)));
                    end
                end
            end
            titl      = ['L-curve with norm: ',norm_type];
            [x_reg,lambda_opt,loss_val] = Lcurve_with_Norm_pinv(A,b,BH2,titl,plotON,norm_type);
            cLcurveAll.creg_H2             = x_reg;
            cLcurveAll.creg_H2_lambda_opt  = lambda_opt;
        case 'gRKHS'  % ad hoc RKHS with a given reproducing kernel
            if ~exist('phi','var')
               B_grkhs = eye(n); 
            end
            [x_reg,lambda_opt,loss_val] = Lcurve_with_Norm_pinv(A,b,BH2,titl,plotON,norm_type);
            cLcurveAll.creg_gRKHS             = x_reg;
            cLcurveAll.creg_gRKHS_lambda_opt  = lambda_opt;
    end
end
end



