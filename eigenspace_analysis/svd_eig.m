% SVD and eig are similar for symmetric positive definite matrix. 

% No singularities
rng('default')
X=randn(1000,5); A = X'*X; 
sv=svd(X);
eigv=sort(eig(A).^.5,'descend');
disp([sv eigv]) 

% With singularity 
rng('default')
X=randn(1000,5);
X = [.01*X(:,1:3), 1000*X(:,4), 1e6*X(:,5)]; A = X'*X; 
sv2=svd(X);
eigv2=sort(eig(A).^.5,'descend');
disp(log10([sv2 eigv2])) 



%%% GSVD. 

% No singularities
rng('default')
X = randn(1000,5); B = X'*X; 
Y = rand(1000,5);  A = Y'*Y; 
Y = [.01*Y(:,1:3), 1000*Y(:,4), 1e2*X(:,5)]; A = Y'*Y; 
[Ugsvd,V,Xgsvd,C,S]    = gsvd(A,B);     % 
%       A = U*C*X'
%       B = V*S*X'
%       C'*C + S'*S = I    >>> 
[Croot_gsvd, indgsvd] = sort(diag(C),'descend');   Ugsvd = Ugsvd(:, indgsvd);
Sroot_gsvd = diag(S); 
Xgsvd      = Xgsvd(:, indgsvd);

eig_gsvd2  = sort(sqrt(diag(C'*C)./diag(S'*S)),'descend');
[V_eig, eigAB]  = eig(A, B); 
[eigAB, indAB]  = sort(diag(eigAB),'descend'); 
V_eig           = V_eig(:, indAB);

disp(log10([eig_gsvd2 eigAB,Croot_gsvd,Sroot_gsvd])) 



% With singularity 
rng('default')
X=randn(1000,5);
X = [.01*X(:,1:3), 1000*X(:,4), 1e6*X(:,5)]; A = X'*X; 
sv2=svd(X);
eigv2=sort(eig(A).^.5,'descend');
disp(log10([sv2 eigv2])) 


