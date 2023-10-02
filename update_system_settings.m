function sysInfo = update_system_settings(sysInfo)

sysInfo.xgrid = linspace(sysInfo.lb, sysInfo.rb, sysInfo.xn)';
sysInfo.tgrid = linspace(0, sysInfo.T, sysInfo.tn)';
sysInfo.dt = sysInfo.tgrid(2) - sysInfo.tgrid(1);
sysInfo.dx = sysInfo.xgrid(2) - sysInfo.xgrid(1);


%% Generate Data
tn = sysInfo.tn;
xn = sysInfo.xn;

tgrid = sysInfo.tgrid;
xgrid = sysInfo.xgrid;
dx = sysInfo.dx;
phi = sysInfo.phi;


L_operator = zeros(tn, xn);
for i = 1:tn
    t = tgrid(i);
    L_operator(i, :) = phi(t, xgrid)*dx;
end
sysInfo.L_operator = L_operator;


end