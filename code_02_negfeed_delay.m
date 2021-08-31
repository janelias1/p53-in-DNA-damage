% Mathematical model of p53 
% For the Special Issue "The Functional Landscape of p53"
% in International Journal of Molecular Sciences
% (Eds. Andreas Prokesch and Jelena Krstic)
% By Jan Elias and Cicely K. Macnamara
% Date: 23.07.2021

% Simple p53-Mdm2 negative feedback with delay

global g_param;

T1 = 60; % for DDE
timeint2 = [0,T1];

% X = [p53], Y = [Mdm2]

% initial concentrations
X0 = 0.2;
Y0 = 0.1;


initdata = [X0;Y0];

% parameters
ks = 2;
k1 = 2; 
K1 = 0.1;
dx = 1; 
n = 2; 
k2 = 2; 
K2 = 1; 
dy = 1;
g_param = [ks, k1, K1, dx, n, k2, K2, dy];

options = odeset('RelTol',1e-8,'AbsTol',1e-8); % options for the ODE/DDE solver

% set the delay in DDE
lags = 3; % lag in Mdm2-dependent degradation of p53
history = initdata;

sol = dde23(@DDE_01,lags,history,timeint2,options);
tsoldde1 = sol.x'; udde1 = sol.y';

% plot results
plotsol(tsoldde1,udde1,lags);

% plot phase plane
plotpp(lags,timeint2,options);


% ========================  Nested Functions  ===========================

% ====================  Differential Equations  =========================

function du = DDE_01(t,u,Z)
% ODE set up

global g_param;

u2lag = Z;
% ks, k1, K1, dx, n, k2, K2, dy

ks = g_param(1); 
k1 = g_param(2); 
K1 = g_param(3); 
dx = g_param(4); 
n = g_param(5); 
k2 = g_param(6); 
K2 = g_param(7); 
dy = g_param(8);

du = zeros(2,1);

% equations
% X = [p53]
du(1)=ks - k1*u2lag(2)*u(1)./(K1+u(1)) - dx*u(1);
% y = [Mdm2]
du(2)= (k2*u(1).^n)./(K2^n+u(1).^n) - dy*u(2);
end


% ========================  Plot solution  ==============================


function plotsol(tsoldde1,udde1,lags)

% Define own colours
newcolors = [0    0.4470    0.7410
    0.8500    0.3250    0.0980];

global g_param;
% ks, k1, K1, dx, n, k2, K2, dy

ks = g_param(1); 
k1 = g_param(2); 
K1 = g_param(3); 
dx = g_param(4); 
n = g_param(5); 
k2 = g_param(6); 
K2 = g_param(7); 
dy = g_param(8);

figure;
axes('FontSize',22);
ax = subplot(1,1,1); hold on;
plot(ax,tsoldde1,udde1(:,1),'LineWidth',3,'Color',newcolors(1,:),'DisplayName','p53');
plot(ax,tsoldde1,udde1(:,2),'LineWidth',3,'Color',newcolors(2,:),'DisplayName','Mdm2');
xlabel('Time','FontSize',24,'interpreter','latex'); 
ylabel('Concentration','FontSize',24,'interpreter','latex');
legend('show','Interpreter','latex')
ax.TickLabelInterpreter = 'latex';
hold off;

% plot nullclines
t1 = find(tsoldde1==lags); t2 = length(tsoldde1); 
tint = tsoldde1(t1:t2); udde1d = interp1(tsoldde1, udde1,tint);
x1 = udde1d(:,1);
x2 = udde1(:,1);
xdde1=udde1(:,1); 
ydde1=udde1(:,2);

figure;
axes('FontSize',22);
ax = subplot(1,1,1); hold on;
plot(ax,x1,(ks-dx*x1).*((K1+x1)./(k1*x1)),'LineWidth',3,'Color',newcolors(1,:),'DisplayName','p53-nullcline');
plot(ax,x2,(k2/dy)*(x2.^n)./(K2^n+x2.^n),'LineWidth',3,'Color',newcolors(2,:),'DisplayName','Mdm2-nullcline');
plot(ax,xdde1(1:end),ydde1(1:end),'-','LineWidth',2,'Color','k','HandleVisibility','off');
ylim([0,2]);
set(findobj(gcf,'type','axes'),'FontSize',22);
xlabel('p53','FontSize',24,'interpreter','latex'); 
ylabel('Mdm2','FontSize',24,'interpreter','latex');
legend('show','Interpreter','latex')
ax.TickLabelInterpreter = 'latex';
hold off;

end

% =======================  plot phase plane  ============================

function plotpp(lags,timeint2,options)

L=2;
NX=5;
dx=0;

X0=zeros(1,4*NX-3); Y0=X0;

X0(1:NX) = linspace(dx,L-dx,NX);
Y0(1:NX) = dx;

X0(NX:2*NX-1) = L-dx;
Y0(NX:2*NX-1) = linspace(dx,L-dx,NX);

X0(2*NX-1:3*NX-2) = flip(linspace(dx,L-dx,NX));
Y0(2*NX-1:3*NX-2) = L-dx;

X0(3*NX-2:4*NX-4) = dx;
dxd = flip(linspace(dx,L-dx,NX)); Y0(3*NX-2:4*NX-4) = dxd(1:end-1);

X0(5) = 1.8;
Y0(1) = 0.2;
Y0(9) = 1.8;

X0(4*NX-3) = 0.5; Y0(4*NX-3)=0.5;
%X0(4*NX+2) = 0.5; Y0(4*NX+2)=1.5;


headWidth = 10;
headLength = 16;
LineLength = 1;


figure; axes('FontSize',22,'TickLabelInterpreter','latex'); hold on;
for i=1:4*NX-3
    
    history = [X0(i),Y0(i)];
    
    sol = dde23(@DDE_01,lags,history,timeint2,options);
    tsol = sol.x'; udde1 = sol.y';
    plot(udde1(:,1),udde1(:,2),'color','k','LineWidth',1);
    X = udde1(:,1); U = X(2:end); X = X(1:end-1); dX = U-X;
    Y = udde1(:,2); V = Y(2:end); Y = Y(1:end-1); dY = V-Y;

    % draw arrow heads #1
    tsol2 = find(abs(tsol - 0.25)<0.01);
    if length(tsol2) > 1
        tsol2 = tsol2(1);
    end
    X1 = X(tsol2); Y1 = Y(tsol2); dX1 = dX(tsol2); dY1 = dY(tsol2);   
    ah = annotation('arrow','color','k', ... 
        'headStyle','vback1','HeadLength',headLength,'HeadWidth',headWidth);
    set(ah,'parent',gca);
    set(ah,'position',[X1 Y1 LineLength*dX1 LineLength*dY1]);
    
    % draw arrow heads #2
    tsol2 = find(abs(tsol - 6)<0.01);
    if length(tsol2) > 1
        tsol2 = tsol2(1);
    end
    X1 = X(tsol2); Y1 = Y(tsol2); dX1 = dX(tsol2); dY1 = dY(tsol2);
    ah = annotation('arrow','color','k', ... 
        'headStyle','vback1','HeadLength',headLength,'HeadWidth',headWidth);
    set(ah,'parent',gca);
    set(ah,'position',[X1 Y1 LineLength*dX1 LineLength*dY1]);
    
    disp(i)
    
end
history = [X0(1),Y0(1)];
sol = dde23(@DDE_01,lags,history,timeint2*3,options);
tsol = sol.x'; udde1 = sol.y';
plot(udde1(12000:end,1),udde1(12000:end,2),'color','k','LineWidth',3);
axis([0 L 0 L]);
box on;
xlabel('p53','FontSize',24,'interpreter','latex'); 
ylabel('Mdm2','FontSize',24,'interpreter','latex');
hold off

end