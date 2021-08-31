% Use this function to plot concentrations of p53, Mdm2 and Mdm2 mRNA as
% the integrated concentrations over the whole cell.
%
% This functions reads the data saved in the text file info_conc_cytnuc.txt
% created by the FF+ code code_04_negfeed_2Dspace.edp

% Define own colours
newcolors = [0    0.4470    0.7410
    0.8500    0.3250    0.0980
    0.9290    0.6940    0.1250];


% path to the file info_conc_cytnuc.txt, for example:
filename = '/Users/janelias/FreeFem++/info_conc_cytnuc.txt';
data = importdata(filename,',');


figure;
axes('FontSize',22);
ax = subplot(1,1,1); hold on;
plot(ax,data(:,1),data(:,2)+data(:,3),'LineWidth',3,'Color',newcolors(1,:),'DisplayName','p53');
plot(ax,data(:,1),data(:,4)+data(:,5),'LineWidth',3,'Color',newcolors(2,:),'DisplayName','Mdm2');
plot(ax,data(:,1),data(:,6)+data(:,7),'LineWidth',3,'Color',newcolors(3,:),'DisplayName','Mdm2 mRNA');
xlabel('Time','FontSize',24,'interpreter','latex'); 
ylabel('Concentration','FontSize',24,'interpreter','latex');
legend('show','Interpreter','latex')
ax.TickLabelInterpreter = 'latex';
hold off;

figure;
axes('FontSize',22);
ax = subplot(1,1,1); hold on;
plot(ax,data(:,1),data(:,2),'LineWidth',3,'Color',newcolors(1,:),'DisplayName','nuc p53');
plot(ax,data(:,1),data(:,4),'LineWidth',3,'Color',newcolors(2,:),'DisplayName','nuc Mdm2');
plot(ax,data(:,1),data(:,6),'LineWidth',3,'Color',newcolors(3,:),'DisplayName','nuc Mdm2 mRNA');
xlabel('Time','FontSize',24,'interpreter','latex'); 
ylabel('Concentration','FontSize',24,'interpreter','latex');
legend('show','Interpreter','latex')
ax.TickLabelInterpreter = 'latex';
hold off;

figure;
axes('FontSize',22);
ax = subplot(1,1,1); hold on;
plot(ax,data(:,1),data(:,3),'LineWidth',3,'Color',newcolors(1,:),'DisplayName','cyt p53');
plot(ax,data(:,1),data(:,5),'LineWidth',3,'Color',newcolors(2,:),'DisplayName','cyt Mdm2');
plot(ax,data(:,1),data(:,7),'LineWidth',3,'Color',newcolors(3,:),'DisplayName','cyt  Mdm2 mRNA');
xlabel('Time','FontSize',24,'interpreter','latex'); 
ylabel('Concentration','FontSize',24,'interpreter','latex');
legend('show','Interpreter','latex')
ax.TickLabelInterpreter = 'latex';
hold off;

% limit cycle
figure; axes('FontSize',22,'TickLabelInterpreter','latex'); hold on;
plot(data(:,2),data(:,4),'LineWidth',2,'Color','k');
hold off;
xlabel('nuclear p53','FontSize',24,'interpreter','latex'); 
ylabel('nuclear Mdm2','FontSize',24,'interpreter','latex');

