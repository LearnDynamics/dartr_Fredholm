function [y_fine_noisy,L_operator_fine,sysInfoFine] = generate_fine_data(sysInfo,f_true,nsr)
%% get fine y-data to represent continuous solution 
% Then downsample in time for discrete data; 

% in discrete: sysInfo.tn = 500; for [0,5]
sysInfoFine    = sysInfo; 
sysInfoFine.tn = 10*sysInfo.tn;
sysInfoFine    = update_system_settings(sysInfoFine); 

L_operator_fine = sysInfoFine.L_operator; 
y_true = L_operator_fine*f_true; 
y_true_norm = sqrt(sum(y_true.^2)*sysInfo.dt/sysInfo.T); % to refine the noise setting later
y_fine_noisy      = y_true + y_true_norm*nsr*randn(sysInfoFine.tn, 1);   

