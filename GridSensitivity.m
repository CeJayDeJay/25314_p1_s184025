dz_sens = [0.25 0.5 1 1.25 1.5];
figure;
for i = 1:5
param.dz = dz_sens(i);
dz = param.dz;
param.zm = 200 ; 
param.z = param.dz/2:param.dz:(param.zm-param.dz/2) ;
z = param.z;
param.nGrid = length(param.z); 
[t,Y] = model4(dz);
P = Y(:,[1:param.nGrid]);
plot(P(3000,:),-param.z,'LineWidth', 6/i/1.5)
hold on
end 
legend('\Deltaz = 0.25','\Deltaz = 0.5', '\Deltaz = 1', '\Deltaz = 1.25', '\Deltaz = 1.5','FontSize', 12)
title('Converged solution over different grid spacing, time = 3000 days','FontSize', 20)
xlabel('Phytoplankton [mmol Nitrogen m⁻3]','FontSize', 16)
ylabel('Depth [m]','FontSize', 16)

% figure;
% for i = 1:5
% param.dz = dz_sens(i);
% dz = param.dz;
% param.zm = 200 ; 
% param.z = param.dz/2:param.dz:(param.zm-param.dz/2) ;
% z = param.z;
% param.nGrid = length(param.z); 
% [t,Y] = model4(dz);
% N = Y(:,[param.nGrid+1:param.nGrid*2]);
% plot(N(3000,:),-param.z)
% hold on
% end 
% legend('\Deltaz = 0.25','\Deltaz = 0.5', '\Deltaz = 1', '\Deltaz = 1.25', '\Deltaz = 1.5','FontSize', 12)
% title('Converged solution over different grid spacing, time = 3000 days','FontSize', 20)
% xlabel('Nutrients [mmol Nitrogen m⁻3]','FontSize', 16)
% ylabel('Depth [m]','FontSize', 16)
% 
% figure;
% for i = 1:5
% param.dz = dz_sens(i);
% dz = param.dz;
% param.zm = 200 ; 
% param.z = param.dz/2:param.dz:(param.zm-param.dz/2) ;
% z = param.z;
% param.nGrid = length(param.z); 
% [t,Y] = model4(dz);
% D = Y(:,[param.nGrid*2+1:end]);
% plot(D(3000,:),-param.z)
% hold on
% end 
% legend('\Deltaz = 0.25','\Deltaz = 0.5', '\Deltaz = 1', '\Deltaz = 1.25', '\Deltaz = 1.5','FontSize', 12)
% title('Converged solution over different grid spacing, time = 3000 days','FontSize', 20)
% xlabel('Detritus [mmol Nitrogen m⁻3]','FontSize', 16)
% ylabel('Depth [m]','FontSize', 16)
