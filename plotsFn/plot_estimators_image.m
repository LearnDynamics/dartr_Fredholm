function plot_estimators_image(est_array,tgrid,L_operator,y,y_true,lgnd,new_figure)
% plot the image of the estimators

colors; % get linestyle and dark colors

if ~exist('new_figure','var'); new_figure=1; end
if new_figure==1; figure; end 
[n_vec,n_est] = size(est_array); 
for n=1:n_est-1
    f_est = est_array(:,n);
    plot(tgrid, L_operator*f_est,linestyle{n},'Color',dred_do_db(n,:),'linewidth',1);hold on;
end
plot(tgrid, y_true,'k:','linewidth',2);hold on;
scatter(tgrid,y,'+','MarkerEdgeColor',[0.4660 0.6740 0.1880],'MarkerEdgeAlpha',0.6)
% plot(tgrid, y,':','linewidth',1,'Color',[0.4660 0.6740 0.1880,0.9]); % [0.2 0.5 0.9 0.2]);


lgnd = [lgnd(1:end-1),'True','Observed',];
xlim([tgrid(1),tgrid(end)*1.01]); 
title('Estimated y');legend(lgnd); legend('location','best')
xlabel('t'); ylabel('y(t)');
end