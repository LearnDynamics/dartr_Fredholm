function sysInfo = system_settings()
% basic settings of the weighted transformation
%{
 Dicrete-Fredholm integral equation:
      \int_lb^rb K(t,x) f(x)dx + noise = y(t)
 Goal: Given y(t) at discrete-times, to estimate f
%}


lb = 1;     
rb = 5;      % Space Integration Interval  [lb,rb]
xn = 100;  


T = 5;      % Time integration interval  [0,T]
tn = 500;   % Time observation steps



% f_true_func = @(x) (sin(x-6)).^2-3;   % True f
% f_true_func = @(x) (x-6).^2+3;   % True f

% f_true_func = @(x) 15*((sin(x-6)).^2-3);   % True f
% f_true_func = @(x) 15*((sin(2*x-6)).^2-3);
% sigma = 0.001;    % noise level

%%
sysInfo.lb = lb;
sysInfo.rb = rb;
sysInfo.T = T;
sysInfo.tn = tn;
sysInfo.xn = xn;

%% select kernel types 
sysInfo.kernel_type = 'Vogel';  % 'exp'; 'poly'; 
switch sysInfo.kernel_type
    case 'exp'
        sysInfo.phi = @(t, x) x.^(-2).*exp(-x.*t); % This phi is the integral kernel K(t,x)
        % sysInfo.kernel_type = 'exp'; 
    case 'poly'
        sysInfo.phi = @(t, x) x.^(-1).*abs(sin(x.*t +1)); % This phi is the integral kernel K(t,x)
        % sysInfo.kernel_type = 'poly'; 
    case 'Vogel'      % Vogel's example: Vogel1996: non-convergence of the L-curve regularization parameter selection method
        fun_a = @(z) 30*abs(z).*(1-abs(z)); % 30*z.*(1-z); %  periodic, z in [0,1] % 30*abs(z).*(1-abs(z)); 
        sysInfo.phi = @(t, x) fun_a(t-x); % This phi is the integral kernel K(t,x)
       % sysInfo.kernel_type = 'Vogel';     
       sysInfo.lb = 0;       sysInfo.rb = 1; sysInfo.xn = tn; 
       sysInfo.T  = 1;       sysInfo.tn = tn; 
       % sysInfo.tn = 100;  sysInfo.xn = 100;
end 


%% 
% f_true_func = @(x) (sin(x-6)).^2-3;   % True f
% f_true_func = @(x) (x-6).^2-1;   % True f
% f_true_func = @(x) 15*((sin(x-6)).^2-3);   % True f
% f_true_func = @(x) 15*((sin(2*x-6)).^2-3);
% sigma = 0.1;    % noise level

% sysInfo.f_true_func = f_true_func;
% sysInfo.sigma = sigma;

sysInfo = update_system_settings(sysInfo);


end