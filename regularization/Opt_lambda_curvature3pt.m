function [lambda] = Opt_lambda_curvature3pt(zeta,ita,all_lambda,plotON,titl)
% find the maximal curvature of the curve (zeta,ita) using 3pt circle 
zeta   = log10(zeta)/2;
ita   = log10(ita)/2;

ind_loss = find(zeta>-20); 
ind_norm = find(ita(ind_loss)>-20); 
zeta1    = zeta(ind_loss);  zeta1 = zeta1(ind_norm); 
ita1     = ita(ind_loss);   ita1  = ita1(ind_norm); 
all_lambda = all_lambda(ind_loss); all_lambda = all_lambda(ind_norm); 

X = [zeta1, ita1];
[~,R2,K2] = curvature3pt(X);      % using 3-point circle radium for curvature 
sgn = -1 + 2 * (K2(:,1)>0).*(K2(:,2)>0);

[val, ind] = max(1./R2.*sgn);
lambda = all_lambda(ind);


if plotON
    figure;h_fig = subplot(121);
    h = plot(zeta1,ita1, '-o', 'LineWidth',2); grid on; axis equal;    %set(h,'marker','.');
    ylabel log_{10}(||x||_{B});  xlabel log_{10}(||Ax-b||);    title('2D curve with curvature vectors');   hold on
    quiver(zeta1,ita1,K2(:,1),K2(:,2),'LineWidth',1);   hold off
    if nargin<5; titl = 'L-curve and normal vector'; end 
    title(titl);
    legend('L-curve','normal vector')
    sgn = -1 + 2 * (K2(:,1)>0).*(K2(:,2)>0);

   %  set(h_fig, 'position', [0.13 0.1 0.3 0.8] ); set(gca,'FontSize',15 );
    
    h_fig = subplot(122);
    semilogx(all_lambda(ind), val, '*','Color', '#D95319','MarkerSize',12, 'LineWidth',5 );hold on;
    semilogx(all_lambda, 1./R2 .*sgn,'Color','#0072BD', 'LineWidth',2);
    xlabel('\lambda'); ylabel('curvature');legend(['\lambda_0 = ', num2str(all_lambda(ind))],'Location','best')
    title('Signed curvature')
    
%     sgtitle('Use L-curve to find the optimal \lambda_0','FontSize',25)
    % set(h_fig, 'position', [0.6 0.175 0.3 0.647] )     set(gca,'FontSize',25 );
    % fig = gcf;    fig.Units = 'inches';     fig.Position = [2 2   14 12];
    
end
end