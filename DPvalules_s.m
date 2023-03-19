% model of NPZ including detritus and seasonal light intensity 
% parameter values from "D>P values" in article doi:10.1016/j.pocean.2007.09.002
% made 14-03-2023, last changed 19-03-2023 
dz = 1;
% 
% %parameters 
param.zm = 200 ; % meters % water column depth
param.dz = dz ; % Grid spacing (m)
param.z = param.dz/2:param.dz:(param.zm-param.dz/2) ;
param.nGrid = length(param.z);  % no. of grid cells

param.HI = 20;%[umol photons m⁻2 s⁻1] 
param.HN = 0.349; %[mmol N/m³]

param.kp = 0.033; %[m² (mmol N)^-1] self shading
param.kw = 0.045; %[m⁻1] &background turbidity
param.Av = 6.4*10^(-5)*(60*60*24); %[m^2 day^-1] 
param.loss = 0.023; %[day⁻1] 
param.gamma = 1.728; %[m^3 (mmol N)^-1 day^-1] grazing mortality 
param.zt = 20; %meters % thermocline depth in stratified waters 
param.I0 = 600; %[umol photons m⁻2 s⁻1] 
param.pmax = 0.5; %[day⁻1] 
param.NB = 31.84; %nutrient concentration at bottom [mmol N m^-3]
param.tau = 0.086; % [day^-1]
param.w = 5.54; % [m day^-1] sinking velocity
    
%Initial conditions
P0 = exp(-(param.z-param.zm/2).^2/5); %[mmol Nitrogen m⁻3]
N0 = repmat(param.NB,1,param.nGrid); %[mmol Nitrogen m^-3]
D0 = repmat(0,1,param.nGrid); %[mmol Nitrogen m^-3]
IC = [P0,N0,D0];

tSpan = [0:6000]; %[days]
[t, Y] = ode23(@dYmodel4, tSpan, IC, [], param);
P = Y(:,[1:param.nGrid]);
N = Y(:,[param.nGrid+1:param.nGrid*2]);
D = Y(:,[param.nGrid*2+1:end]);

%% plots 
figure;
x1 = N(5000,:);
y1 = -param.z;
x2 = P(5000,:)*10^3;
x3 = D(5000,:)*10^3;
tl = tiledlayout(1,1);
ax1 = axes(tl);
plot(ax1,x1,y1,'-g','LineWidth', 3)
ylabel('Depth [m]','FontSize',12)
xlabel('mmol nitrogen m^-3, concentration scale for N','FontSize',12)
ax1.XColor = 'g';
ax1.YColor = 'k';
legend('Nutrients','FontSize',12)
ax2 = axes(tl);
plot(ax2,x2,y1,'r',x3,y1,'b','LineWidth', 3)
ylabel('Depth [m]','FontSize',12)
xlabel('umol nitrogen m^-3, concentration scale for P & D','FontSize',12)
ax2.XColor = 'k';
ax2.YColor = 'k';
ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
ax2.Color = 'none';
legend('Phytoplankton','Detritus','FontSize',12)
title('B. SS-sol for the state variables with seasonality of I, time=5000 days','FontSize',12)

figure; 
clf
surface(t,-param.z,P')
shading flat
colorbar
title('Phytoplankton growth over time [mmol Nitrogen m⁻3]','FontSize', 20)
xlabel('Time [days]','FontSize', 16)
ylabel('Depth [meters]','FontSize', 16)
clim([0 0.25]);
% % 
% % 
figure; 
clf
surface(t, -param.z, N')
shading flat
colorbar
title('Nutrient concentration over time [mmol Nitrogen m⁻3]','FontSize', 20)
xlabel('Time [days]','FontSize', 16)
ylabel('Depth [meters]','FontSize', 16)

figure; 
clf
surface(t, -param.z, D')
shading flat
colorbar
title('Detritus concentration over time [mmol Nitrogen m⁻3]','FontSize', 20)
xlabel('Time [days]','FontSize', 16)
ylabel('Depth [meters]','FontSize', 16)


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


        I = Imodel(P,D,t,param);
        g = param.pmax.*(I./(param.HI+I)).*(N./(N+param.HN));

        dPdt = dPdt' + g.*P-param.loss*P-param.gamma*P.^2;

        dDdt = dDdt' + param.loss*P + param.gamma*P.^2 - param.tau*D; 

        dNdt = dNdt' - g.*P + param.tau*D; %- param.w*dDdt;

        dYdt = [dPdt;dNdt;dDdt];
       
    end 


    function I = Imodel(P,D,t,param)
        season = 0.4118*cos(2*pi/365*t)+1;       
        Q = param.kp*param.dz*(cumsum(P)-P/2+cumsum(D)-D/2) ;
        I = param.I0*exp(-param.kw*param.z'-Q)*season ;
    end 