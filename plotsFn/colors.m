%% Nice Colors for plots
% % .........................................................................
Color  = [zeros(3,1) eye(3)];
Color(:,1) = [85;170;170]/255;
Color(:,2) = [60;60;230]/255;
Color(:,3) = [170;0;170]/255;
Color(:,4) = [200;0;0]/255;
Color(:,5) = [225;125;0]/255;
Color(:,6) = [120;160;200]/255;
% .........................................................................

% matlab RGB 
matcolor = zeros(6,3);
matcolor(1,:) = [0 0.4470 0.7410];        % 'darkblue'
matcolor(2,:) = [0.8500 0.3250 0.0980];   % 'darkred'
matcolor(3,:) = [0.9290 0.6940 0.1250];   % 'darkyellow'
matcolor(4,:) = [0.4940 0.1840 0.5560];   % 'darkpurple'
matcolor(5,:) = [0.4660 0.6740 0.1880];   % 'darkgreen'
matcolor(6,:) = [0.3010 0.7450 0.9330];   % 'lightblue' 
code_str = {'#0072BD', '#D95319','#EDB120',	'#7E2F8E', '#77AC30','#4DBEEE','#D95319'};
name_str = {'darkblue','darkorange','darkyellow','darkpurple','darkgreen','lightblue','darkred'};

dark_color.matcolor = matcolor; 
dark_color.code_str = code_str; 
dark_color.name_str = name_str;

dred_do_db = matcolor([2,3,1],:); 


linestyle = {'-','-.','--',':'};  %  A,B,C=the best, true

linestyles = {'--x','-.o','--d',':'};