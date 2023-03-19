% model of NPZ including detritus 
% parameter values from "reference values" in article doi:10.1016/j.pocean.2007.09.002
% made 21-02-2023, last changed 19-03-2023 
%function [t,Y] = model4(dz)
% use it as a function to test grid sensitivity 

dz = 1;

%parameters 
param.zm = 200 ; % meters % water column depth
param.dz = dz ; % Grid spacing (m)
param.z = param.dz/2:param.dz:(param.zm-param.dz/2) ;
param.nGrid = length(param.z);  % no. of grid cells

param.HI = 20;%[umol photons m⁻2 s⁻1] 
param.HN = 0.3; %[mmol N/m³]

param.kp = 0.05; %[m² (mmol N)^-1] self shading
param.kw = 0.0375; %[m⁻1] &background turbidity
param.Av = 5*10^(-5)*(60*60*24); %[m^2 day^-1] 
param.loss = 0.03; %[day⁻1] 
param.gamma = 1.5; %[m^3 (mmol N)^-1 day^-1] grazing mortality 
param.zt = 20; %meters % thermocline depth in stratified waters 
param.I0 = 600; %[umol photons m⁻2 s⁻1] 
param.pmax = 0.5; %[day⁻1] 
param.NB = 30; %nutrient concentration at bottom [mmol N m^-3]
param.tau = 0.1; % [day^-1]
param.w = 15; % [m day^-1] sinking velocity

%Initial conditions
P0 = exp(-(param.z-param.zm/2).^2/5); %[mmol Nitrogen m⁻3]
N0 = repmat(param.NB,1,param.nGrid); %[mmol Nitrogen m^-3]
D0 = repmat(0,1,param.nGrid); %[mmol Nitrogen m^-3]
IC = [P0,N0,D0];

tSpan = [0:3500]; %[days]
[t, Y] = ode23(@dYmodel4, tSpan, IC, [], param);
P = Y(:,[1:param.nGrid]);
N = Y(:,[param.nGrid+1:param.nGrid*2]);
D = Y(:,[param.nGrid*2+1:end]);

%% sensitivity plots
% % % different grazing rates
% gamma = [0.5 1 1.5 2 2.5];
% figure; 
% for i = 1:5
% param.gamma=gamma(i);
% [t, Y] = ode23(@dYmodel4, tSpan, IC, [], param);
% P = Y(:,[1:param.nGrid]);
% N = Y(:,[param.nGrid+1:param.nGrid*2]);
% D = Y(:,[param.nGrid*2+1:end]);
% 
% plot(P(3000,:),-param.z, 'LineWidth', i/1.5)
% hold on
% end
% title('Phytoplankton over depth with different grazing rates, time = 3000 days','FontSize', 20)
% legend('Gamma=0.5 [m^3 (mmol N)^-1 day^-1]', 'Gamma=1 [m^3 (mmol N)^-1 day^-1]', 'Gamma=1.5 [m^3 (mmol N)^-1 day^-1]', 'Gamma=2 [m^3 (mmol N)^-1 day^-1]', 'Gamma=2.5 [m^3 (mmol N)^-1 day^-1]', 'FontSize', 12)
% xlabel('Phytoplankton [mmol Nitrogen m⁻3]','FontSize', 16)
% ylabel('Depth [m]','FontSize', 16)
% 
% 
% figure; 
% for i = 1:5
% param.gamma=gamma(i);
% [t, Y] = ode23(@dYmodel4, tSpan, IC, [], param);
% P = Y(:,[1:param.nGrid]);
% N = Y(:,[param.nGrid+1:param.nGrid*2]);
% D = Y(:,[param.nGrid*2+1:end]);
% 
% plot(N(3000,:),-param.z, 'LineWidth', i/1.5)
% hold on
% end
% title('Nutrients over depth with different grazing rates, time = 3000 days','FontSize', 20)
% legend('Gamma=0.5 [m^3 (mmol N)^-1 day^-1]', 'Gamma=1 [m^3 (mmol N)^-1 day^-1]', 'Gamma=1.5 [m^3 (mmol N)^-1 day^-1]', 'Gamma=2 [m^3 (mmol N)^-1 day^-1]', 'Gamma=2.5 [m^3 (mmol N)^-1 day^-1]','FontSize', 12)
% xlabel('Nutrients [mmol Nitrogen m⁻3]','FontSize', 16)
% ylabel('Depth [m]','FontSize', 16)
% % 
% figure;
% for i = 1:5
% param.gamma=gamma(i);
% [t, Y] = ode23(@dYmodel4, tSpan, IC, [], param);
% P = Y(:,[1:param.nGrid]);
% N = Y(:,[param.nGrid+1:param.nGrid*2]);
% D = Y(:,[param.nGrid*2+1:end]);
% 
% plot(D(3000,:),-param.z, 'LineWidth', i/1.5)
% hold on
% end
% title('Detritus over depth with different grazing rates, time = 3000 days','FontSize', 20)
% legend('Gamma=0.5 [m^3 (mmol N)^-1 day^-1]', 'Gamma=1 [m^3 (mmol N)^-1 day^-1]', 'Gamma=1.5 [m^3 (mmol N)^-1 day^-1]', 'Gamma=2 [m^3 (mmol N)^-1 day^-1]', 'Gamma=2.5 [m^3 (mmol N)^-1 day^-1]','FontSize', 12)
% xlabel('Detritus [mmol Nitrogen m⁻3]','FontSize', 16)
% ylabel('Depth [m]','FontSize', 16)
% 
% % % Different bottom nutrient values
% bottomN = [20 25 30 35 40];
% figure; 
% for i = 1:5
% param.NB=bottomN(i);
% [t, Y] = ode23(@dYmodel4, tSpan, IC, [], param);
% P = Y(:,[1:param.nGrid]);
% N = Y(:,[param.nGrid+1:param.nGrid*2]);
% D = Y(:,[param.nGrid*2+1:end]);
% 
% plot(P(3000,:),-param.z, 'LineWidth', i/1.5)
% hold on
% end
% title('Phytoplankton over depth with different nutrients concentrations at the bottom (NB), time = 3000 days','FontSize', 20)
% legend('NB=20 [mmol N m^-3]','NB=25 [mmol N m^-3]','NB=30 [mmol N m^-3]','NB=35 [mmol N m^-3]','NB=40 [mmol N m^-3]','FontSize', 12)
% xlabel('Phytoplankton [mmol Nitrogen m⁻3]','FontSize', 16)
% ylabel('Depth [m]','FontSize', 16)
% 
% figure; 
% for i = 1:5
% param.NB=bottomN(i);
% [t, Y] = ode23(@dYmodel4, tSpan, IC, [], param);
% P = Y(:,[1:param.nGrid]);
% N = Y(:,[param.nGrid+1:param.nGrid*2]);
% D = Y(:,[param.nGrid*2+1:end]);
% 
% plot(N(3000,:),-param.z, 'LineWidth', i/1.5)
% hold on
% end
% title('Nutrients over depth with different nutrients concentrations at the bottom (NB), time = 3000 days','FontSize', 20)
% legend('NB=20 [mmol N m^-3]','NB=25 [mmol N m^-3]','NB=30 [mmol N m^-3]','NB=35 [mmol N m^-3]','NB=40 [mmol N m^-3]','FontSize', 12)
% xlabel('Nutrients [mmol Nitrogen m⁻3]','FontSize', 16)
% ylabel('Depth [m]','FontSize', 16)
% 
% figure; 
% for i = 1:5
% param.NB=bottomN(i);
% [t, Y] = ode23(@dYmodel4, tSpan, IC, [], param);
% P = Y(:,[1:param.nGrid]);
% N = Y(:,[param.nGrid+1:param.nGrid*2]);
% D = Y(:,[param.nGrid*2+1:end]);
% 
% plot(D(3000,:),-param.z, 'LineWidth', i/1.5)
% hold on
% end
% title('Detritus over depth with different nutrients concentrations at the bottom (NB), time = 3000 days','FontSize', 20)
% legend('NB=20 [mmol N m^-3]','NB=25 [mmol N m^-3]','NB=30 [mmol N m^-3]','NB=35 [mmol N m^-3]','NB=40 [mmol N m^-3]','FontSize', 12)
% xlabel('Detritus [mmol Nitrogen m⁻3]','FontSize', 16)
% ylabel('Depth [m]','FontSize', 16)


%% plots 

% figure; 
% clf
% surface(t,-param.z,P')
% shading flat
% colorbar
% title('Phytoplankton growth over time [mmol Nitrogen m⁻3]','FontSize', 20)
% xlabel('Time [days]','FontSize', 16)
% ylabel('Depth [meters]','FontSize', 16)
% clim([0 0.3]);
% % 
% % 
% figure; 
% clf
% surface(t, -param.z, N')
% shading flat
% colorbar
% title('Nutrient concentration over time [mmol Nitrogen m⁻3]','FontSize', 20)
% xlabel('Time [days]','FontSize', 16)
% ylabel('Depth [meters]','FontSize', 16)
% 
% figure; 
% clf
% surface(t, -param.z, D')
% shading flat
% colorbar
% title('Detritus concentration over time [mmol Nitrogen m⁻3]','FontSize', 20)
% xlabel('Time [days]','FontSize', 16)
% ylabel('Depth [meters]','FontSize', 16)

% 
% figure; 
% yyaxis left 
% plot(-param.z,P(end,:))
% ylabel('Phytoplankton [cells m⁻3]')
% yyaxis right 
% plot(-param.z,N(end,:))
% ylabel('Nutrients [mmol nutrients m⁻3]')
% xlabel('Depth [meters]')

% figure; 
% yyaxis left 
% plot(-param.z,P(end,:))
% ylabel('Phytoplankton [cells m⁻3]')
% yyaxis right 
% plot(-param.z,Imodel(P(end,end),D(end,end),param))
% ylabel('Light intensity [umol photons m⁻2 s⁻1]')
% xlabel('Depth [meters]')
% 
% figure; 
% yyaxis left 
% plot(-param.z,N(end,:))
% ylabel('Nutrients [mmol nutrients m⁻3]')
% yyaxis right 
% plot(-param.z,Imodel(P(end,end), D(end,end),param))
% ylabel('Light intensity [umol photons m⁻2 s⁻1]')
% xlabel('Depth [meters]')

% figure; 
% yyaxis left 
% plot(-param.z,P(3000,:)*10^3,-param.z,D(3000,:)*10^3)
% ylabel('Concentrations [umol nitrogen m^-3]')
% yyaxis right 
% plot(-param.z,N(3000,:))
% legend('Phytoplankton','Detritus','Nutrients')
% ylabel('Nutrients [mmol nitrogen m^-3]')
% xlabel('Depth [m]')

% figure;
% x1 = N(3000,:);
% y1 = -param.z;
% x2 = P(3000,:)*10^3;
% x3 = D(3000,:)*10^3;
% tl = tiledlayout(1,1);
% ax1 = axes(tl);
% plot(ax1,x1,y1,'g','LineWidth', 3)
% ylabel('Depth [m]','FontSize',12)
% xlabel('mmol nitrogen m^-3, concentration scale for N','FontSize',12)
% ax1.XColor = 'g';
% ax1.YColor = 'k';
% legend('Nutrients','FontSize',12)
% ax2 = axes(tl);
% plot(ax2,x2,y1,'r',x3,y1,'b','LineWidth', 3)
% ylabel('Depth [m]','FontSize',12)
% xlabel('umol nitrogen m^-3, concentration scale for P & D','FontSize',12)
% ax2.XColor = 'k';
% ax2.YColor = 'k';
% ax2.XAxisLocation = 'top';
% ax2.YAxisLocation = 'right';
% ax2.Color = 'none';
% legend('Phytoplankton','Detritus','FontSize',12)
% title('A. Steady state solution for the state variables, time=3000 days','FontSize',12)

% II = Imodel(P(3000,:),D(3000,:),param);
% IItest = (II(:,:)./(param.HI+II(:,:)));
% figure;
% xx1 = (N(3000,:)./(N(3000,:)+param.HN));
% yy1 = -param.z;
% xx2 = (II(:,171)./(param.HI+II(:,171)));
% xx3 = (II(:,end)./(param.HI+II(:,end)));
% xx4 = (II(:,200)./(param.HI+II(:,200)));
% plot(xx1,yy1,xx2,yy1,xx3,yy1,xx4,yy1)
% legend('Nutrients','II1','II2','II3')
% xlabel('Limiting coefficient, 0=0%, 1=100%')
% ylabel('Depth [m]')
% title('Limiting factors at time=3000')

% figure; 
% semilogx(P(3000,:),-param.z,N(3000,:),-param.z,D(3000,:),-param.z)
% xlabel('Concentrations [mmol nutrients m⁻3]','FontSize',16)
% ylabel('Depth [m]','FontSize',16)
% title('Steady state solutions for the state variables, time=3000 days','FontSize',20)
% legend('Phytoplankton','Nutrients','Detritus','Fontsize',12)
% 
% figure; 
% semilogx(P(1500,:),-param.z,N(1500,:),-param.z,D(1500,:),-param.z)
% xlabel('Concentrations [mmol nutrients m⁻3]','FontSize',16)
% ylabel('Depth [m]','FontSize',16)
% title('Not steady state solutions for the state variables, time=1500 days','FontSize',20)
% legend('Phytoplankton','Nutrients','Detritus','Fontsize',12)

% figure; 
% semilogx(P(3500,:),-param.z,N(3500,:),-param.z,D(3500,:),-param.z)
% xlabel('Concentrations [mmol nutrients m⁻3]','FontSize',16)
% ylabel('Depth [m]','FontSize',16)
% title('Steady state solutions for the state variables, time=3500 days','FontSize',20)
% legend('Phytoplankton','Nutrients','Detritus','Fontsize',12)

%% functions
    function dYdt = dYmodel4(t, Y, param)  
        P = Y(1:param.nGrid);
        N = Y(param.nGrid+1:param.nGrid*2);
        D = Y(param.nGrid*2+1:param.nGrid*3);

        ix = 2:(param.nGrid);

        JPdiff(ix) = -param.Av*(P(ix)-P(ix-1))/param.dz;
        JPdiff(1) = 0; 
        JPdiff(param.nGrid+1) = 0;  
        JP = JPdiff;
        dPdt = -(JP(2:(param.nGrid+1))-JP(1:param.nGrid))/param.dz; 


        JNdiff(ix) = -param.Av*(N(ix)-N(ix-1))/param.dz;
        JNdiff(1) = 0; 
        JNdiff(param.nGrid+1) = -param.Av.*((param.NB-N(end))./param.dz);  
        %JNdiff(param.nGrid+1) = -param.Av.*((N(end)+param.NB)./param.dz); 
        JN = JNdiff;
        dNdt = -(JN(2:(param.nGrid+1))-JN(1:param.nGrid))/param.dz;


        JDadv(ix) = param.w*D(ix-1);
        JDadv(1) = 0;  
        JDadv(param.nGrid+1) = 0; 
        JDdiff(ix) = -param.Av*(D(ix)-D(ix-1))/param.dz;
        JDdiff(1) = 0; 
        JDdiff(param.nGrid+1) = param.w*D(param.nGrid);   
        JD = JDadv + JDdiff;
        dDdt = -(JD(2:(param.nGrid+1))-JD(1:param.nGrid))/param.dz; 


        I = Imodel(P,D,param);
        g = param.pmax.*(I./(param.HI+I)).*(N./(N+param.HN));

        dPdt = dPdt' + g.*P-param.loss*P-param.gamma*P.^2;

        dDdt = dDdt' + param.loss*P + param.gamma*P.^2 - param.tau*D; 

        dNdt = dNdt' - g.*P + param.tau*D; %- param.w*dDdt;

        dYdt = [dPdt;dNdt;dDdt];
       
    end 


    function I = Imodel(P,D,param)
        Q = param.kp*param.dz*(cumsum(P)-P/2+cumsum(D)-D/2) ;
        
        I = param.I0*exp(-param.kw*param.z'-Q) ;
    end 
%end