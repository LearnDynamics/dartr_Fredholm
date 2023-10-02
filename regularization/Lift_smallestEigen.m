
function B = Lift_smallestEigen(B,tol)
% lift the smallest eigenvalue of a matrix to a tolerannce
[U,eigB]            = eig(B); % faster than SVD: [U,eigB,V]  svd(B);   % B = U*S*V'; U'*U = I;   V'*V= I; U=V since B is symmetric 
eigBval             = diag(eigB); 
if min(eigBval) <tol
    eigBval(eigBval<tol)  = tol;
    fprintf('\n Warning: some eigenvalues of regularization matrix Bmat are less than %e. Set them to be it. \n', tol);
end
B = B;  fprintf('\n DO nonthing to Bmat'); 
% B = U*diag(eigBval)*U';

end

