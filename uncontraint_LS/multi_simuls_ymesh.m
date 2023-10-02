


function [err_L2rho_projAB,err_cells,err_l2_projA,loss_array,disc_est] = multi_simuls_ymesh(f_true,A,B,V_AB,V_A,r,xgrid,rho,eigAB,sysInfo,ymesh_seq,nsr,normType,saveDIR)
if exist(saveDIR,'file')
    load(saveDIR);    
    disc_est.err_L2rho_projAB = err_L2rho_projAB_disc;
    disc_est.err_cells    = err_cells_disc; 
    disc_est.err_l2_projA = err_l2_projA_disc; 
    disc_est.loss_array   = loss_array_disc; 
    return; 
end
% multiple simulations: 
n_simuls = 100; 
% normType  = {'l2','L2','RKHS'}; 
err_L2rho_projAB = zeros(n_simuls,length(normType),length(ymesh_seq)); 
err_l2_projA     =  0*err_L2rho_projAB;  
err_cells        = cell(n_simuls,length(ymesh_seq)); 
loss_array       = 0*err_L2rho_projAB; 


err_L2rho_projAB_disc = zeros(n_simuls,length(normType),length(ymesh_seq)); 
err_l2_projA_disc     =  0*err_L2rho_projAB;  
err_cells_disc        = cell(n_simuls,length(ymesh_seq)); 
loss_array_disc       =  0*err_L2rho_projAB; 
plotON    =0; 

L_operator = sysInfo.L_operator; 

y_true       = L_operator*f_true; 

y_true_norm  = sqrt(sum(y_true.^2)*sysInfo.dt/sysInfo.T); % to refine the noise setting later


for m=1:length(ymesh_seq)
    gap  = ymesh_seq(m);
    ymesh_ind  = 1:gap:sysInfo.tn;
    L_operator_downsample = L_operator(ymesh_ind,:);
    
    y_true1 = y_true(ymesh_ind); sqrt_deltat = sqrt(sysInfo.dt*gap); 
    parfor n=1:n_simuls
        y   = y_true1 + y_true_norm*nsr*randn(length(y_true1), 1)*sqrt_deltat;      
        b   = L_operator_downsample'*y;
        
        [err_L2rho_projAB(n,:,m), err_cells{n,m},err_l2_projA(n,:,m),loss_array(n,:,m)] = reguLSE_1dataset(A,b,B,f_true,V_AB,V_A,r,xgrid,rho,eigAB,plotON,normType);
        
        Abar = L_operator_downsample'*L_operator_downsample; 
        [err_L2rho_projAB_disc(n,:,m), err_cells_disc{n,m},err_l2_projA_disc(n,:,m),loss_array_disc(n,:,m)] = reguLSE_1dataset(Abar,b,B,f_true,V_AB,V_A,r,xgrid,rho,eigAB,plotON,normType);
    end
    fprintf('Progress %i out of %i \n ',  m,length(ymesh_seq));
end
fprintf('DONE \n'); 

disc_est.err_L2rho_projAB = err_L2rho_projAB_disc;
disc_est.err_cells    = err_cells_disc; 
disc_est.err_l2_projA = err_l2_projA_disc; 
disc_est.loss_array   = loss_array_disc; 
if exist('saveDIR','var') 
    save(saveDIR,'err_cells','err_L2rho_projAB','err_l2_projA','loss_array',...
        'err_cells_disc','err_L2rho_projAB_disc','err_l2_projA_disc','loss_array_disc','disc_est'); 
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


