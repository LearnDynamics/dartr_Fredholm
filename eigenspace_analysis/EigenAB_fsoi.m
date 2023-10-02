function [V_A,eigA,V_L, eigL,r]= EigenAB_fsoi(A,B, plotON,method)
%% plot the eigenvalues of A and generalized (A,B)

if ~exist('method','var')
    method = 'svdA'; % 'svdA' not eigA; eig(A,B) not svd(A/B)
    %  using eig(A) leads to similar V, but SVD ensures >= eigenvalues;
    %  eig(A,B) is the right one to get A*V = B*V*D, V'*B*V = I 
end
switch method
    case 'svdA'  % eig(A,B) for G-eig, but uses SVD(A)
        % [V_A, eigA]   = eig(A);      % A*V = V*D.     V'*V = I 
         [Usvd,eigA,Vsvd] = svd(A);      [eigA, indA] = sort(diag(eigA),'descend');   V_A = Usvd(:, indA);

        [V_L, eigL0]  = eig(A, B);   % A*V = B*V*D.   V'*B*V = I;  % --- it may have negative small eigenvalues 
        [eigA, indA]  = sort(diag(eigA),'descend');
        [eigL, indAB] = sort(diag(eigL0),'descend');     
        V_A           = V_A(:, indA);
        V_L           = V_L(:, indAB);
    case 'svdAB'   % for G-eig(A,B) % A = U*S*V';  U is "same" as V_A (different sign); but S is nonnegative
        [Usvd,eigA,Vsvd] = svd(A);      [eigA, indA] = sort(diag(eigA),'descend');   V_A = Usvd(:, indA);

        % use SVD to compute eigen(A,B)
        rootB     = chol(B);    % B = R'*R.
        A_rootB   = rootB*A*rootB; 
        [U2,S2,~] = svd(A_rootB); 
        V_L2  = rootB\U2; 
        figure; 
        subplot(121);  imagesc(V_L2'*V_L2);      title('V^t*V:   svd(rootB*A*rootB)'); 
        subplot(122);  imagesc(V_L2'*B*V_L2);    title('V^t*B*V: svd(rootB*A*rootB)'); 
        
        figure; imagesc(A*V_L2- B*V_L2*S2); % large error in A equation; 

        [Veig_L, eigAB]  = eig(A, B);   % A*V = B*V*D.   V'*B*V = I;  % --- it may have negative small eigenvalues 
        [eigAB, indAB] = sort(diag(eigAB),'descend');          V_L           = Veig_L(:, indAB);
         figure; imagesc(A*V_L- B*V_L*diag(eigAB)); 

        V_L = V_L2; eigL = diag(S2); 
end

% plot the eigenvalues 
    figure; 
    eigA_abs = abs(eigA(eigA>0));  r_eigA = length(eigA_abs); 
    eigL_abs = abs(eigL(eigL>0));  r_eigL = length(eigL_abs); 
    r   = min(r_eigA,r_eigL);
    semilogy(1:r,eigA_abs(1:r),'r-x','linewidth',1); hold on;
    semilogy(1:r,eigL_abs(1:r),'b-o','linewidth',1);
    legend('Eigenvalue A','Eigenvalue (A,B)');% title('Eigenvalue of A and (A,B)')
    xlabel('i'); ylabel('\lambda_i');
    figname = ['eigen_val_vec',method];
    if exist('figname','var')
        set_positionFontsAll;   print([figname,'.pdf'],'-dpdf', '-bestfit');
    end  
    
% r = rank(A);  % the threshold is 1e-16
    r = length(find(abs(eigL)>1e-8)); 


if plotON ==1
    figure;
    subplot(131)
    semilogy(abs(eigA(eigA>0)),'r-x','linewidth',1); hold on;
    semilogy(abs(eigL(eigL>0)),'b-o','linewidth',1);
    legend('Eigenvalue A','Eigenvalue (A,B)');title('Eigenvalue of A and (A,B)')
    xlabel('i'); ylabel('\lambda_i');
    
    subplot(132)
    plot(V_A(:,1:min(r,6)),'linewidth',1); title('Eigenvectors A');
    xlabel('k'); % ylabel('V');
    subplot(133)
    plot(V_L(:, 1:min(r,6)),'linewidth',1); title('Eigenvectors (A,B)');
    xlabel('k'); % ylabel('V');
end

end 


function compuareSVD_eig(A,B)
% compare SVD with eig, just checking the numerical difference 
[U,S,V]      = svd(A); %  diagonal matrix S, with nonnegative diagonal elements in decreasing order, and unitary matrices U and V so that A = U*S*V'.
S = diag(S); 

[V_A, eigA]  = eig(A); 
[eigA, indA] = sort(diag(eigA),'descend');
V_A          = V_A(:, indA); 

    figure;
    subplot(131)
    semilogy(abs(eigA),'r-x','linewidth',1); hold on;
    semilogy(abs(S),'b-o','linewidth',1);      legend('Eig','SVD'); 
    
    subplot(132)  % the eigen-vectors are the same, but with different signs 
    plot(V_A(:,1:3),'linewidth',1); title('eig vect A')
    subplot(133)
    plot(-U(:, 1:3),'linewidth',1); title('eig vect SVD')
    
    
%    subplot(132);  imagesc(V_A); title('V in Eig'); subplot(133); imagesc(U); title('U in SVD'); 
    
[V_L, eigL0]  = eig(A, B);
[eigL, indAB] = sort(diag(eigL0),'descend'); 
V_L           = V_L(:, indAB);


[U_L, S_L, V_svd]  = svd(B\A);   % here B is diagnal, so B\A is sysmetric 
figure; subplot(131); semilogy(abs(eigL),'r-x','linewidth',1); hold on;
        semilogy(S_L,'b-o','linewidth',1);      legend('EigAB','svd((invB)A)');
        subplot(132); plot(V_L(:,1:3),'linewidth',1); title('eig(A,B)')
        subplot(133); plot(U_L(:, 1:3),'linewidth',1); title('svd((invB)A)')

    
end


function compare_Geig(A,B)
        [Usvd,eigA,Vsvd] = svd(A);      [eigA, indA] = sort(diag(eigA),'descend');   V_A = Usvd(:, indA);

        % use gsvd counterpart --- it has correct sv, but not eigenvector. V'*B*V = I. 
        [Usvd_AinvB, svdAinvB,Vsvd_AinvB]  = svd(A/B);     % A/B = U*S*V.          -  the eigenvectors are not B-onb    
        [svdAinvB, indABsvd] = sort(diag(svdAinvB),'descend');   Vsvd_AinvB = Vsvd_AinvB(:, indABsvd);

        
         [Veig_AinvB, eigAinvB]  = eig(A/B);    % not good: A/B not symmetric, may have negative/complex eigv
         [eigAinvB, indABinv] = sort(diag(eigAinvB),'descend'); 
         Veig_AinvB = Veig_AinvB(:, indABinv);

        % use SVD to compute eigen(A,B)
        rootB     = chol(B);    % B = R'*R.
        A_rootB   = rootB*A*rootB; 
        [U2,S2,~] = svd(A_rootB); eigL = diag(S2); 
        V_L2  = rootB\U2; 
        figure; 
        subplot(121);  imagesc(V_L2'*V_L2);      title('V^t*V:   svd(rootB*A*rootB)'); 
        subplot(122);  imagesc(V_L2'*B*V_L2);    title('V^t*B*V: svd(rootB*A*rootB)'); 

        [Veig_L, eigAB]  = eig(A, B);   % A*V = B*V*D.   V'*B*V = I;  % --- it may have negative small eigenvalues 
        [eigAB, indAB] = sort(diag(eigAB),'descend');          V_L           = Veig_L(:, indAB);

        figure; 
           semilogy(abs(eigA),'r-x','linewidth',1); hold on;
           semilogy(abs(eigAB),'c-x','linewidth',1); hold on;
           semilogy(abs(eigAinvB),'b-*','linewidth',1);    
           semilogy(svdAinvB,'c-o','linewidth',1);  
           semilogy(S2,'k:d','linewidth',1);
           legend('eigA','EigAB','eigAinvB', 'svd(AinvB)','svd(rootB*A*rootB)');
           title('Generalized eigenvalues and singular values')

        figure; imagesc(V_L -V_L2); % the eigenvectors are very different.  

        figure; 
        subplot(121); imagesc(A*V_L2- B*V_L2*S2); title('AV-BVS: svd(rootB*A*rootB)')% large error in A equation; 
        subplot(122); imagesc(A*Vsvd_AinvB- B*Vsvd_AinvB*svdAinvB); title('AV-BVS: svd(A/B)') % large error in A equation;
        % subplot(133); imagesc(A*Veig_AinvB- B*Veig_AinvB*eigAinvB); title('AV-BVS: eig(A/B)') % large error in A equation;
end

