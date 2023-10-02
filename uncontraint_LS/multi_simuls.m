


function [err_L2rho_projAB,err_cells,err_l2_projA,loss_array] = multi_simuls(f_true,A,B,V_AB,V_A,r,xgrid,rho,eigAB,sysInfo,nsr_seq,normType,saveDIR)
if exist(saveDIR,'file')
    load(saveDIR); return; 
end
% multiple simulations: 
n_simuls = 100; 
% normType  = {'l2','L2','RKHS'}; 
err_L2rho_projAB = zeros(n_simuls,length(normType),length(nsr_seq)); 
err_cells = cell(n_simuls,length(nsr_seq)); 
plotON    =0; 

L_operator = sysInfo.L_operator;
%% nsr = 0: only 1 dataset since there is no noise;
nsr = 0; sqrt_deltat = sqrt(sysInfo.dt); 
tic
y_true = L_operator*f_true; 
y_true_norm = sqrt(sum(y_true.^2)*sysInfo.dt/sysInfo.T);
y      = y_true + y_true_norm*nsr*randn(sysInfo.tn, 1)*sqrt_deltat;
b      = L_operator'*y;

[err_L2rho_projAB_nrs0,errors_nrs0] = reguLSE_1dataset(A,b,B,f_true,V_AB,V_A,r,xgrid,rho,eigAB,plotON,normType);
t1  = toc; 
fprintf('nsr = 0 with 1 dataset. Time elapsed: %2.2f minutes \n', t1/60); 
fprintf('Total time to expect: %2.2f minutes / your worker in parallel \n ', t1/60*n_simuls* length(nsr_seq));


err_l2_projA   = err_L2rho_projAB; 
loss_array     = err_L2rho_projAB; 

%% nsr > 0: 100 datasets   
for m = 1:length(nsr_seq) 
    nsr = nsr_seq(m)*sqrt_deltat;  tn =sysInfo.tn;
    parfor n=1:n_simuls
        y      = y_true + y_true_norm*nsr*randn(tn, 1);
        b      = L_operator'*y;
        [err_L2rho_projAB(n,:,m), err_cells{n,m},err_l2_projA(n,:,m),loss_array(n,:,m)] = reguLSE_1dataset(A,b,B,f_true,V_AB,V_A,r,xgrid,rho,eigAB,plotON,normType);
    end
    fprintf('Progress %i out of %i \n ',  m,length(nsr_seq)); 
end
fprintf('DONE \n'); 

if exist('saveDIR','var') 
    save(saveDIR,'err_cells','err_L2rho_projAB','err_L2rho_projAB_nrs0','errors_nrs0','err_l2_projA','loss_array'); 
end

end



%% computing estimator and its errors 
function [err_L2rho,errors,err_l2,loss_array] = reguLSE_1dataset(A,b,B,f_true,V_AB,V_A,r,xgrid,rho,eigAB,plotON,normType)
if ~exist('plotON','var'); plotON=1;end
% plotON = 1; 
%% test different regularizations: l2, L2, RKHS, and may be H1? 
[~,~,cLcurveAll] = reg_Lcurve_3in1(A,b,B,plotON,normType); 
est_array        = [cLcurveAll.creg_l2,cLcurveAll.creg_L2,cLcurveAll.creg_RKHS,f_true];
% lgnd      = {'Est-l2','Est-L2','Est-RKHS','True'}; 
% plot_estimators(est_array,lgnd,xgrid,rho,1);
[errors]         = compute_estimator_Error(B,V_AB,eigAB,V_A,r,est_array,rho,xgrid); 
err_L2rho        = errors.L2rho_l2_projAB(1,:);
err_l2           = errors.L2rho_l2_projA(2,:);
loss_array       = [cLcurveAll.creg_l2_loss_val,cLcurveAll.creg_L2_loss_val,cLcurveAll.creg_RKHS_loss_val];
end


