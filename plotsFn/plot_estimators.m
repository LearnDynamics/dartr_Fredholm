function plot_estimators(est_array,lgnd,xi,rho,new_figure)
% plot the estimators 
% Input 
%     Est_array: each colum is an estimator, with the last column = true


colors; % get linestyle and dark colors

[n_vec,n_est] = size(est_array); 
true_val      = est_array(:,end);
if ~exist('new_figure','var'); new_figure=1; end
if new_figure==1; figure; end 
for n=1:n_est-1
    plot(xi,est_array(:,n),linestyle{n},'Color',dred_do_db(n,:),'linewidth',1); hold on;
end
plot(xi,true_val,'linestyle',linestyle{n_est},'Color','k','linewidth',2); 

 
if min(true_val) >0 ;  lb = min(true_val)*0.8; 
else;                  lb = min(true_val)*1.2; 
end
if  max(true_val)>0;   ub = max(true_val)*1.2; 
else;                  ub = max(true_val)*0.8; 
end 
yyaxis left; ylim([lb,ub]);
xlabel('u'); ylabel('f(u)'); 

% super-impose rho
yyaxis right;  
normalizer = 1; % max(abs(true_val))/2/max(rho); 
rho_length = length(rho); nx = length(xi); gap = ceil(nx/rho_length); indx = 1:gap:nx; xirho = xi(indx); % sparser index
area(xirho,rho(indx)*normalizer, 'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.3, 'EdgeAlpha', 0,'Displayname', '\rho');  
% area(xirho,rho(indx)*normalizer, 'FaceColor', 'c', 'FaceAlpha', 0.2,'EdgeAlpha', 0,'Displayname', '\rho');
yticklabels('auto'); h=gca;  set(h,'YColor',[0.7 0.7 0.7]); 
lgnd = [lgnd,'rho']; 
legend(lgnd,'Location','best'); 


end
