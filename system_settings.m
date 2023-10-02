function sysInfo = system_settings()
% basic settings of the weighted transformation
%  

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
sysInfo.phi = @(t, x) x.^(-2).*exp(-x.*t);


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