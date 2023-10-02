function lambda_seq = lambda_seq_from_regMats(A_bar,B,N)
%% 1. get range of lambda from eigenvalues of A_bar or the generalized eigenvalue
tol = 1e-16; 
% l2 eigenvalue: eig(A_bar)
[UofA,eigA,~]    = svd(A_bar); eigA  = diag(eigA); 
[~,ind ]   = sort(eigA,'descend');   UofA = real(UofA(:,ind)); eigA = eigA(ind); % eigA sorted,  not necessary, just for convienience
EigenA      = eigA(eigA>tol);
minEigen    = max( min(EigenA), tol); maxEigen = min(max(EigenA(5:end)),0.5);   % since maxEigen\c_i, we set lambda <0.5
% L2 eigenvalue eig(A_bar,B)-- generalize eigenvalue 
 [V,eigL]    = eig(A_bar,B); eigL = diag(eigL); 
[~,ind ]   = sort(eigL,'descend');   V = real(V(:,ind)); eigL = eigL(ind); % eigL sorted, not necessary, just for convienience



minEigenL   = max(min(eigL), 1e-10);  maxEigenL = min(max(eigL(5:end)),0.1);    
minEigen    = min(minEigenL,minEigen); maxEigen = max([maxEigen,maxEigenL]);

temp        = linspace(log10(minEigen), log10(maxEigen),N);  
lambda_seq= [0 10.^temp];
end