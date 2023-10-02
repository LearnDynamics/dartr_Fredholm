function [A_con_bar,b_con_bar,rho_con] = LS_Ab_con_sum1(A0,b0,rho,con_type,dx)
% Add strong constraint sum_j x_j *dx = 1 to least squares problem A0 x= b0
%{
INPUT
  A0,b0      - original rectangular regression matrix and vector
  con_type   - 'strong' or 'weak'
Ouput
  A,b        - square regression matrix and vector with the constraint
%}

[m,n]    = size(A0); 
if ~exist('con_type','var'); con_type = 'strong'; end
switch con_type
    case 'strong'
        A_con = zeros(m,n-1);  b_con = b0- A0(:,n)/dx;  % change the regression matrix and vector: rectangular now!
        for i=1:m
            A_con(i,:) =  A0(i,1:n-1)- A0(i,n);
        end
        A_con_bar = A_con'*A_con;    % change to square matrix  ~~~~~ function space may mismatch here!
        b_con_bar = A_con'*b_con;
         % the exploration measure 
         rho_con = rho(1:n-1)+rho(n); rho_con = rho_con/(dx*sum(rho_con)); 
    case 'weak'
        A_con = [A0; ones(1,n)];  b_con = [b0;1*dx];  % change the regression matrix and vector: rectangular now!
        A_con_bar = A_con'*A_con;    % change to square matrix  ~~~~~ function space may mismatch here!
        b_con_bar = A_con'*b_con;
        rho_con = rho; 
end


end