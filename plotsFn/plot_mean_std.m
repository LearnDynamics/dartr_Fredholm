function h = plot_mean_std(x_ind,data_array,newfigure,lgnd,label_y,figname,xlabelstr) 
% plot the ensemble mean with std bar of samples in data_array:  X(samples,types,x_ind)
%{
    x_ind      - the index variable for the data 
    data_array -  n_samples x n_types x n_ind
%}
[n_samples,n_types,n_ind] = size(data_array); 
% [n_types,n_ind,n_samples] = size(data_array);
std_data    = zeros(n_types,n_ind);
mean_data   = zeros(n_types,n_ind);
for nn= 1:n_ind
    temp     = squeeze(data_array(:,:,nn));
    std_data(:,nn) = std(temp,0,1);
    mean_data(:,nn)= mean(temp,1);
    % boxplot(temp,'Notch','on');
end
fprintf('\n Mean of %i simulations \n ',n_samples);
disp(mean_data);


if max(max(mean_data))/min(min(mean_data)) > 1e2 || max(mean_data(2,:)./mean_data(3,:)) > 1e2  % use log-scale plot
    std_data_ub = log10(mean_data+std_data);
    y           = mean_data-std_data; 
    std_data_lb = sign(y).*log10(abs(y)); 
    mean_data   = log10(mean_data); 
    std_data    = max(std_data_ub-mean_data,std_data_lb-mean_data);
    label_y = ['log_{10} ', label_y];
end

if exist('newfigure','var') && newfigure ==1; figure; end 
colors;     % colors and linestyples 
for i= 1:n_types
    h = errorbar(x_ind,mean_data(i,:),std_data(i,:),linestyles{i},'linewidth',1.5,'Color',dred_do_db(i,:)); hold on;
end
lb_x = min(x_ind)- 0.05*(x_ind(end)-x_ind(1)); 
ub_x = max(x_ind)+0.05*(x_ind(end)-x_ind(1)); 
xlim([lb_x,ub_x ]); 

% special treatment to the 2nd row that has too large values 
% range_val = min(mean_data([1,3],:) -std_data([1,3],:));  % remove the 2nd row: it has too large values 
range_val = min(min(mean_data -std_data),min(mean_data));  % all rows
if min(range_val) >0 ;  lb = min(range_val)*0.8; 
else;                   lb = min(range_val)*1.2; 
end
% range_val = max(mean_data([1,3],:) +std_data([1,3],:)); 
range_val = max(max(mean_data +std_data),max(mean_data));  % all rows
if  max(range_val)>0;   ub = max(range_val)*1.2; 
else;                   ub = max(range_val)*0.8; 
end 
ylim([lb,ub]); xlabel(xlabelstr);

xticklabel= cell(1,length(x_ind)); 
for i= 1:length(x_ind)
   xticklabel{i}=  num2str(round(10^x_ind(i),3)); % sprintf('%04d',10^x_ind(i));
   % xticklabel{'0.125','0.25','0.5','1','2'}; 
end
xticks(x_ind); xticklabels(xticklabel); 

if exist('lgnd','var'); legend(lgnd,'Location','best');  end
if exist('label_y','var'); ylabel(label_y); yticklabels('auto'); end
if exist('figname','var') && ~isempty(figname)
    set_positionFontsAll;   print([figname,'.pdf'],'-dpdf', '-bestfit'); 
end
end


