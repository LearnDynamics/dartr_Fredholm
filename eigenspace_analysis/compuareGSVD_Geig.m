function compuareGSVD_Geig(A,B)
% Conclusion: we should use eig(A,B) to get V'*B*V = I 
% Compare methods computing generalized eigenvalue problems
%   gsvd, eig(A,B), eig(A/B), svd(A/B)

fprintf('Compare methods computing the Generalized eigenvectors:\n'); 
fprintf('eig(A,B) is needed for such a decomposition.\n'); 

if nargin ==0 
    rng('default')
    n = 1000; 
    X = randn(n,50); B = X'*X/n;
    Y = rand(n,50);  A = Y'*Y/n;
    Y = [.1*Y(:,1:3), 100*Y(:,4), 1e3*X(:,5:end)]; A = Y'*Y;
end

% 1. gsvd is not good for G-eig
[Ugsvd,V,Xgsvd,C,S] = gsvd(A,B); % --- it is not eig(A,B) % unitary matrices U and V, a (usually) square matrix X, and nonnegative diagonal matrices
%   C and S so that 
%       A = U*C*X'
%       B = V*S*X'
%       C'*C + S'*S = I 
% When B is square and nonsingular, gsvd(A,B) corresponds to the ordinary singular values, svd(A/B)
% Note: U and V are Unitary, hence, neither is g-eigen-vector; 
%       A'*A = X* C.^2*X'  
%       B'*B = X* S.^2*X'

% X in gsvd does not diagonalize A or B 
figure; subplot(121);  imagesc(Xgsvd'*A'*A*Xgsvd);   title('X*A^t*A*X^t:   gsvd(A, B)');
        subplot(122);  imagesc(Xgsvd'*B'*B*Xgsvd);   title('X*B^t*B*X^t:   gsvd(A, B)');
Xgsvd_inv = pinv(Xgsvd*S); 
figure; subplot(131);  imagesc(Xgsvd_inv*A'*A*Xgsvd_inv');   title('Xinv*A^t*A*Xinv^t:   gsvd(A, B)');
        subplot(132);  imagesc(Xgsvd_inv*B'*B*Xgsvd_inv');   title('Xinv*B^t*B*Xinv^t:   gsvd(A, B)');
        subplot(133);  imagesc(Xgsvd_inv*B*Xgsvd_inv');   title('Xinv*B*Xinv^t:   gsvd(A, B)'); 


[Croot_gsvd, indgsvd] = sort(diag(C),'descend');  
Sroot_gsvd = diag(S); 
Ugsvd     = Ugsvd(:, indgsvd);        
Xgsvd     = Xgsvd(:, indgsvd);
eig_gsvd   = sort(diag(C)./diag(S),'descend');   

% check if the eigenvectors are B orthonormal
figure; subplot(121);  imagesc(Xgsvd'*Xgsvd);  
        subplot(122);  imagesc(Xgsvd'*B*Xgsvd); 

% % eig(A,B) is for G-eig
[V_eig, eigL0]  = eig(A, B); 
[eigAB, indAB]  = sort(diag(eigL0),'descend'); 
V_eig           = V_eig(:, indAB);

figure; subplot(131);  imagesc(V_eig'*V_eig);                     title('V^t*V:   eig(A, B)'); 
        subplot(132);  imagesc(V_eig'*B*V_eig);                   title('V^t*B*V=I?: eig(A, B)'); 
        subplot(133);  imagesc(A*V_eig - B*V_eig*diag(eigAB));    title('AV-BVS=0?: eig(A,B)'); colorbar;    



% % svd(A/B) NOT good for G-eig
ABinv      = A/B; 
[Usvd_AinvB, svdAinvB,V_svdAinvB] = svd(ABinv); 
[svdAinvB, indABinv2] = sort(diag(svdAinvB),'descend'); 
Usvd_AinvB = Usvd_AinvB(:, indABinv2);
V_svdAinvB = V_svdAinvB(:,indABinv2); 


figure; subplot(131);  imagesc(V_svdAinvB'*V_svdAinvB);    title('V^t*V:   svd(A/B)');       
        subplot(132);  imagesc(V_svdAinvB'*B*V_svdAinvB);  title('V^t*B*V: svd(A/B)');
        subplot(133);  imagesc(A*V_svdAinvB - B*V_svdAinvB*diag(svdAinvB));  title('AV-BVS=0?: svd(A/B)'); colorbar; 

figure; subplot(131);  imagesc(Usvd_AinvB'*Usvd_AinvB);    title('U^t*U:   svd(A/B)');       
        subplot(132);  imagesc(Usvd_AinvB'*B*Usvd_AinvB);  title('U^t*B*U=I?: svd(A/B)');
        subplot(133);  imagesc(A*Usvd_AinvB - B*Usvd_AinvB*diag(svdAinvB));    title('AU-BUS=0?: U svd(A/B)');  colorbar;  


% % % % eig(A/B) NOT good for G-eig
% % [Veig_AinvB, eigAinvB]  = eig(ABinv);   % A/B may not be symmetric; so the eigenvector may be complex 
% % [eigAinvB, indABinv] = sort(abs(diag(eigAinvB)),'descend'); 
% % Veig_AinvB = Veig_AinvB(:, indABinv);
% % 
% % figure; subplot(131);  imagesc(Veig_AinvB'*Veig_AinvB);                      title('V^t*V: eig(A/B)');       
% %         subplot(132);  imagesc(Veig_AinvB'*B*Veig_AinvB);                    title('V^t*B*V=I?: eig(A/B)');
% %         subplot(133);  imagesc(A*Veig_AinvB - B*Veig_AinvB*diag(eigAinvB));  title('AV-BVS=0?: eig(A/B)'); colorbar; 

% SVD(rootB*A*rootB) to compute G-eig(A,B): good eigv, not good V
rootB     = chol(B);    % B = R'*R. figure; imagesc(rootB'*rootB -B) 
A_rootB   = rootB*A*rootB; 
[U2,S2,~] = svd(A_rootB); eigL = diag(S2); 
V_L       = rootB\U2; 
figure; 
subplot(131);  imagesc(V_L'*V_L);            title('V^t*V:   svd(rootB*A*rootB)'); 
subplot(132);  imagesc(V_L'*B*V_L);          title('V^t*B*V=I?: svd(rootB*A*rootB)'); 
subplot(133);  imagesc(A*V_L - B*V_L*S2);    title('AV-BVS=0?: svd(rootB*A*rootB)'); colorbar; 
        

figure; 
        semilogy(abs(eigAB),'r-x','linewidth',1); hold on;
        % semilogy(eigAinvB,'b-*','linewidth',1);    % 'eigAinvB',
        semilogy(svdAinvB,'c-o','linewidth',1);     
        semilogy(eig_gsvd,'k:d','linewidth',1);
        semilogy(S2,'k:x','linewidth',2)
        legend('EigAB', 'svd(AinvB)','gsvd','svd(rootB*A*rootB)');
        title('Generalized eigenvalues and singular values')



figure;
        subplot(131); plot(V_eig(:,1:3),'linewidth',1); title('eig(A,B)')
        subplot(132); plot(Usvd_AinvB(:, 1:3),'linewidth',1); title('svd(A/B)')
        subplot(133); plot(Xgsvd(:, 1:3),'linewidth',1); title('gsvd(A,B)')
        sgtitle('Generalized eigenvectors and singular vectors')
%       subplot(142); plot(Veig_AinvB(:,1:3),'linewidth',1); title('eig(A/B)')
%         subplot(143); plot(Usvd_invBA(:, 1:3),'linewidth',1); title('svd(invBA)')
%         subplot(144); plot(Veig_invBA(:, 1:3),'linewidth',1); title('eig invBA)')


end