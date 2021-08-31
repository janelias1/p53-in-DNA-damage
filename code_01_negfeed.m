% Mathematical model of p53 
% For the Special Issue "The Functional Landscape of p53"
% in International Journal of Molecular Sciences
% (Eds. Andreas Prokesch and Jelena Krstic)
% By Jan Elias and Cicely K. Macnamara
% Date: 23.07.2021

% Simple p53-Mdm2 negative feedback

global g_param;

T1 = 10;
timeint1 = [0,T1];

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

% solve ODE
options = odeset('RelTol',1e-8,'AbsTol',1e-8); % options for the ODE solver
[tsol1,sol1] = ode23(@ODE_01,timeint1,initdata,options);
[tsol2,sol2] = ode23(@ODE_01,timeint1,initdata*10,options);
[tsol3,sol3] = ode23(@ODE_01,timeint1,[0.5,1.7],options);

% plot results
plotsol(tsol1,sol1,tsol2,sol2,tsol3,sol3);

% plot phase plane
plotpp(timeint1,options)


% ========================  Nested Functions  ===========================

% ====================  Differential Equations  =========================

function du = ODE_01(t,u)
% ODE set up

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

du = zeros(2,1);

% equations
% X = [p53]
du(1)=ks - k1*u(2)*u(1)./(K1+u(1)) - dx*u(1);
% y = [Mdm2]
du(2)= (k2*u(1).^n)./(K2^n+u(1).^n) - dy*u(2);
end

% ========================  Plot solution  ==============================

function plotsol(tsol1,sol1,tsol2,sol2,tsol3,sol3)

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
plot(ax,tsol1,sol1(:,1),'LineWidth',3,'Color',newcolors(1,:),'DisplayName','p53');
plot(ax,tsol1,sol1(:,2),'LineWidth',3,'Color',newcolors(2,:),'DisplayName','Mdm2');
xlabel('Time','FontSize',24,'interpreter','latex');
ylabel('Concentration','FontSize',24,'interpreter','latex');
legend('show','Interpreter','latex')
ax.TickLabelInterpreter = 'latex';
hold off;

% plot nullclines
figure;
x1=sol1(:,1); x2=sol2(:,1); x3=sol3(:,1);
y1=sol1(:,2); y2=sol2(:,2); y3=sol3(:,2);
axes('FontSize',22);
ax = subplot(1,1,1); hold on;
plot(ax,x1,(ks-dx*x1).*((K1+x1)./(k1*x1)),'LineWidth',3,'Color',newcolors(1,:),'DisplayName','p53-nullcline');
plot(ax,x1,(k2/dy)*(x1.^n)./(K2^n+x1.^n),'LineWidth',3,'Color',newcolors(2,:),'DisplayName','Mdm2-nullcline');
plot(ax,x2,(ks-dx*x2).*((K1+x2)./(k1*x2)),'LineWidth',3,'Color',newcolors(1,:),'HandleVisibility','off');
plot(ax,x2,(k2/dy)*(x2.^n)./(K2^n+x2.^n),'LineWidth',3,'Color',newcolors(2,:),'HandleVisibility','off');
plot(ax,x1,y1,'-','LineWidth',2,'Color','k','HandleVisibility','off');
plot(ax,x2(20:end),y2(20:end),'--','LineWidth',2,'Color','k','HandleVisibility','off');
plot(ax,x3,y3,':','LineWidth',2,'Color','k','HandleVisibility','off');
set(findobj(gcf,'type','axes'),'FontSize',22);
xlabel('p53','FontSize',24,'interpreter','latex');
ylabel('Mdm2','FontSize',24,'interpreter','latex');
legend('show','Interpreter','latex')
ax.TickLabelInterpreter = 'latex';
hold off;

end

% =======================  plot phase plane  ============================

function plotpp(timeint1,options)

L=2;
NX=5;
dx=0;

X0=zeros(1,4*NX-4); Y0=X0;

X0(1:NX) = linspace(dx,L-dx,NX);
Y0(1:NX) = dx;

X0(NX:2*NX-1) = L-dx;
Y0(NX:2*NX-1) = linspace(dx,L-dx,NX);

X0(2*NX-1:3*NX-2) = flip(linspace(dx,L-dx,NX));
Y0(2*NX-1:3*NX-2) = L-dx;

X0(3*NX-2:4*NX-4) = dx;
dxd = flip(linspace(dx,L-dx,NX)); Y0(3*NX-2:4*NX-4) = dxd(1:end-1);

Y0(1) = 0.2;
Y0(9) = 1.8;

%col=hsv(4*NX);

headWidth = 10;
headLength = 16;
LineLength = 1;

figure; axes('FontSize',22,'TickLabelInterpreter','latex'); hold on;
for i=1:4*NX-4
    
    initdata = [X0(i),Y0(i)];
    
    [tsol,sol] = ode23(@ODE_01,timeint1,initdata,options);
    plot(sol(:,1),sol(:,2),'color','k','LineWidth',1);
    X = sol(:,1); U = X(2:end); X = X(1:end-1); dX = U-X;
    Y = sol(:,2); V = Y(2:end); Y = Y(1:end-1); dY = V-Y;
    
    % draw arrow heads #1
    tsol2 = find(abs(tsol - 0.2)<0.01);
    if length(tsol2) > 1
        tsol2 = tsol2(1);
    end
    X1 = X(tsol2); Y1 = Y(tsol2); dX1 = dX(tsol2); dY1 = dY(tsol2);
    
    ah = annotation('arrow','color','k', ... 
        'headStyle','vback1','HeadLength',headLength,'HeadWidth',headWidth);
    set(ah,'parent',gca);
    set(ah,'position',[X1 Y1 LineLength*dX1 LineLength*dY1]);

    
    % draw arrow heads #2
    tsol2 = find(abs(tsol - 1)<0.01);
    if length(tsol2) > 1
        tsol2 = tsol2(1);
    end
    X1 = X(tsol2); Y1 = Y(tsol2); dX1 = dX(tsol2); dY1 = dY(tsol2);
    ah = annotation('arrow','color','k', ... %newcolors(i,:), ...
        'headStyle','vback1','HeadLength',headLength,'HeadWidth',headWidth);
    set(ah,'parent',gca);
    set(ah,'position',[X1 Y1 LineLength*dX1 LineLength*dY1]);
    
    disp(i)
    
end

axis([0 L 0 L]);
box on;
xlabel('p53','FontSize',24,'interpreter','latex');
ylabel('Mdm2','FontSize',24,'interpreter','latex');
%set(gca,'FontSize',30,'fontweight','b','fontname','arial')
hold off

end