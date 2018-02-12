function [data]=RecData()
clc;
close all;
%% Call the ReadDataFile
data=ReadDataFile();
%% Plot the Data
figure;

subplot(4,1,1);
plot(data{1,1},data{1,2});
title('RAW DATA: CAM OUT','FontSize',18,'FontName','Arial Narrow')
xlabel('Time(s)','FontSize',14,'FontName','Arial');
ylabel('Voltage (V)','FontSize',14,'FontName','Arial');
grid on;

subplot(4,1,2);
plot(data{1,1},data{1,3});
title('RAW DATA: TRIG OUT','FontSize',18,'FontName','Arial Narrow')
xlabel('Time(s)','FontSize',14,'FontName','Arial');
ylabel('Voltage (V)','FontSize',14,'FontName','Arial');
grid on;

subplot(4,1,3);
plot(data{1,1},data{1,4});
title('Raw Data: Incident Gauge(G_f=2.04;R=1k ohm)','FontSize',18,'FontName','Arial Narrow')
xlabel('Time(s)','FontSize',14,'FontName','Arial');
ylabel('Voltage (V)','FontSize',14,'FontName','Arial');
grid on;

subplot(4,1,4);
plot(data{1,1},data{1,5});
title('Raw Data: Transmitted Gauge(G_f=2.04;R=1k ohm)','FontSize',18,'FontName','Arial Narrow')
xlabel('Time(s)','FontSize',14,'FontName','Arial');
ylabel('Voltage (V)','FontSize',14,'FontName','Arial');
grid on;

hold off;
%% Calculate strain from Voltage data
% Time=data{1,1};
% Inc_strain=HalfBridgeCkt_process(data{1,4},2.05,30);
% Trans_strain=HalfBridgeCkt_process(data{1,5},2.05,30);
% %% Plot Bar Strain Data
% figure;
% CKB_strain1=subplot(3,1,1);
% plot(Time,data{1,2});
% title('Raw Data: Optical Gate o/p (V)','FontSize',18,'FontName','Arial Narrow')
% xlabel('Time(s)','FontSize',14,'FontName','Arial');
% ylabel('Voltage (V)','FontSize',14,'FontName','Arial');
% grid on;
% 
% hold on;
% CKB_strain2=subplot(3,1,2);
% plot(Time,Inc_strain*(10^6));
% title('Raw Data: Incident Bar Strain','FontSize',18,'FontName','Arial Narrow')
% xlabel('Time(s)','FontSize',14,'FontName','Arial');
% ylabel('Strain x 10^{-6} ','FontSize',14,'FontName','Arial');
% grid on;
% 
% CKB_strain3=subplot(3,1,3);
% plot(Time,Trans_strain*(10^6));
% title('Raw Data: Transmitted Bar Strain','FontSize',18,'FontName','Arial Narrow')
% xlabel('Time(s)','FontSize',14,'FontName','Arial');
% ylabel('Strain x 10^{-6} ','FontSize',14,'FontName','Arial');
% grid on;
% 
% linkaxes([CKB_strain1,CKB_strain2,CKB_strain3],'x');
