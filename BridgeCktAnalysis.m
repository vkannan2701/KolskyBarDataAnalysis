%Input RawData and convert it to bar strain and bar stress
function [camout,BarStrain_INC,BarStrain_TR,ti,tf]=BridgeCktAnalysis()
%% INPUT VARIABLES
Gf1=2.04;
Gf2=170;
Vii=30;
Vit=15;
Gain_inc=10;
Gain_tr=0.1;
%% Read Data File: Column 3/9=Incident Bar Gauge, Column 4/9= Transmitted Bar Strain
[inc,tr,~,camout]=ReadDataFile();
%% Plot Data
Fig2=figure;
set(Fig2,'defaulttextinterpreter','tex');
plot(inc(:,1)*10^6,inc(:,2),'b','Linewidth',2);
T2={'Raw Strain gauge data'};
title(T2,'FontSize',20);
xlabel('Time($\mu$s)','FontSize',20,'Interpreter','latex');
ylabel('$Oscilloscope\:o/p (V)$','FontSize',20,'Interpreter','latex');
hold on;
plot(tr(:,1)*10^6,tr(:,2),'r','Linewidth',2);
%% Prompt for ti and tf => start and end points of interest in the signal
% ti=input('Start time (s) : ');
% tf=input('End time (s) : ');
ti=-65e-6;
tf=105e-6;
%% Compute strains from half bridge calc; epsilon=2*Vo/Gf*Vi
BarStrain_INC(:,1)=inc(:,1);
BarStrain_INC(:,2)=inc(:,2)*2.101/(Gf1*Vii*Gain_inc);
BarStrain_TR(:,1)=tr(:,1);
BarStrain_TR(:,2)=tr(:,2)*2.11/(Gf2*Vit*Gain_tr);
%% Plot Strain Data
Fig3=figure;
set(Fig3,'defaulttextinterpreter','latex');
plot(BarStrain_INC(:,1),BarStrain_INC(:,2)*(10^6),'b','Linewidth',2);
T3={'$BAR\:STRAINS$'};
title(T3,'FontSize',20);
xlabel('Time($\mu$s)','FontSize',20,'Interpreter','latex');
ylabel('$Bar\:Strain\:(\mu\epsilon)$','FontSize',20,'Interpreter','latex');
hold on;
plot(BarStrain_TR(:,1),BarStrain_TR(:,2)*(10^6),'r','Linewidth',2);
end
