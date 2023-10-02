%% Get regression vector for each f. 
nsr   = 2;    % noise to signal ratio    
           % -- issue when nsr=0: the optimal lambda=0, but numerical error in inversion prevents us from get to it.  
           % Solution: add lambda =0 estimator, and select between lambda_opt and 0 by min-loss ( a factor (1e2) to be robust) >> estimator

if ~exist('f_true','var')
   % f_true from functions
   % f_true_func = @(x) (sin(x-6)).^2+1;   % True f
  % f_true_func = @(x) 0.7*exp(-(x-2).^2/0.25)*sqrt(1/(2*.5*pi)) +.3*exp(-(x-4).^2/0.09)*sqrt(1/(2*0.3*pi));   % True f
   %  f_true_func = @(x) (x-6).^2+1;   % True f
   % f_true_func = @(x) 15*((sin(x-6)).^2-3);   % True f
   % f_true_func = @(x) 15*((sin(2*x-6)).^2-3);
    f_true_func = @(x) x.^2;   % True f
   f_true = f_true_func(xgrid);
  
   % f_true inside and outside of FSOI:
%   f_true = V_AB(:,10);   % f_true = V_AB(:,5);    % ind>5 ( i.e, eigAB<1e-9): rkhs not good, l2+L2 can slightly tolerate more, bc. not using rkhs inversion).
end

% f_true = f_true/sqrt(dx*sum(f_true.^2));
f_true_norm = sum(f_true.^2.*rho); 

sqrt_deltat = sqrt(sysInfo.dt); 
y_true = L_operator*f_true; 
y_true_norm = sqrt(sum(y_true.^2)*sysInfo.dt/sysInfo.T);
y      = y_true + y_true_norm*nsr*randn(tn, 1)*sqrt_deltat;
b      = L_operator'*y;


%% unconstained LSE with regularuzations: l2, L2, RKHS
plotON = 1;  
if exist('fig_dir','var') 
    figname = [fig_dir,file_str,'_estimator']; 
else 
    clear figname; 
end
%% test different regularizations: l2, L2, RKHS, and may be H1? 
normType = {'l2','L2','RKHS'}; 
[~,~,cLcurveAll] = reg_Lcurve_3in1(A,b,B,plotON,normType); 

if plotON==1
    lgnd      = {'l2','L2','RKHS','True'};
    est_array = [cLcurveAll.creg_l2,cLcurveAll.creg_L2,cLcurveAll.creg_RKHS,f_true];
    new_figure = 0;
    h= figure;
    subplot(121); plot_estimators(est_array,lgnd,xgrid,rho,new_figure);
    subplot(122); plot_estimators_image(est_array,tgrid,L_operator,y,y_true,lgnd,new_figure);   
    if exist('figname','var')
       figure(h); set_positionFontsAll;   print([figname,'.pdf'],'-dpdf', '-bestfit');
    end
    
    figure; 
    subplot(121);
    % projection to the V_AB
    titl = 'EigenAB projection';
    plot_estimators_projection(B,V_AB,2*r,est_array,lgnd,titl,new_figure);
    % project to V_A
    subplot(122);
    titl = 'EigenA projection';
    plot_estimators_projection(eye(xn),V_A,2*r,est_array,lgnd,titl,new_figure);
    if exist('figname','var')
        set_positionFontsAll;   print([figname,'projection.pdf'],'-dpdf', '-bestfit');
    end
end

[errors] = compute_estimator_Error(B,V_AB,eigAB,V_A,r,est_array,rho,xgrid); 

display(errors.L2rho_l2_discription)
display(errors.L2rho_l2)


display(errors.discription_L2rho_l2_projAB)
display(errors.L2rho_l2_projAB)

display(errors.discription_L2rho_l2_projA)
display(errors.L2rho_l2_projA)

fprintf('Loss values: l2,L2,rkhs \n')
loss_array  = [cLcurveAll.creg_l2_loss_val,cLcurveAll.creg_L2_loss_val,cLcurveAll.creg_RKHS_loss_val]

