
function plot_ensemble_quatile(ens_traj,tt,traj_true,ylabel_val)
% plot an ensemble of trajectories with quantile 
[n_t,n_ens] = size(ens_traj); 

% for nn =1:n_ens  % plot the ensemble lines 
%     h1= plot(tt,ens_traj(:,nn), 'c-','linewidth',1); hold on
% end

prctiles= [20 25 70 75 90 95];
%[lh, ph] = fanChart(tt, ens_traj, 'mean', prctiles); 
[lineh, bandsh] = fanChart(tt, ens_traj, 'mean', prctiles,'alpha', .7);
% [lineh, bandsh] = fanChart(tt, ens_traj, 'mean', prctiles, ...    'alpha', .2, 'colormap', {'shadesOfColor', [0 0 .8]});
% txt = strcat({'Percentile'}, cellstr(int2str((20:20:80)')));
% h2 = plot(tt,mean(ens_traj,2), 'r-.','linewidth',1.5); 
hold on; 
h3 = plot(tt,traj_true,'k:','linewidth',2); 
txt = strcat({'Percentile '}, cellstr(int2str([25 75 95]')));
% legend([bandsh;lineh], [txt;{'mean'}])
legend([bandsh;lineh;h3], [txt;{'Ensemble mean'};{'FOM'}],'location','best')

 xlabel('Time'); ylabel(ylabel_val);
% legend([ph,lh, h3],lgnd,'location','best'); 

end