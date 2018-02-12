% Reading Perception Data files
function [incident,trans,camtrig,camout]=ReadDataFile()
clc;
close all;
%% Choose file directory
dir_path=uigetdir;
% search_txt=strcat(dir_path,'\*.txt');
search_csv=strcat(dir_path,'\*.csv');
%% Incident Strain Gauge Data
% Syntax : M=csvread(<filename>,6,0)
[inc_file,inc_path]=uigetfile(search_csv);
inc_full_path=strcat(inc_path,inc_file);
incident=csvread(inc_full_path,6,0);
%% Plot raw data: Incident
Fig1=figure;
set(Fig1,'defaulttextinterpreter','latex');
plot(incident(:,1)*10^6,incident(:,2),'b','Linewidth',2);
T1={'$RAW\:DATA$'};
title(T1,'FontSize',20);
xlabel('Time($\mu$s)','FontSize',20,'Interpreter','latex');
ylabel('$Oscilloscope\:o/p$','FontSize',20,'Interpreter','latex');
hold on;
%% Transmitted Strain Gauge
[trans_file,trans_path]=uigetfile(search_csv);
trans_full_path=strcat(trans_path,trans_file);
trans=csvread(trans_full_path,6,0);
%% Plot transmitted pulse
plot(trans(:,1)*10^6,trans(:,2),'r','LineWidth',2);
%% Cam trig
[camtrig_file,camtrig_path]=uigetfile(search_csv);
camtrig_full_path=strcat(camtrig_path,camtrig_file);
camtrig=csvread(camtrig_full_path,6,0);
%% Plot
plot(camtrig(:,1)*10^6,camtrig(:,2)/50,'Color',[0 0 0],'LineWidth',0.5);
%% Cam output
[camout_file,camout_path]=uigetfile(search_csv);
camout_full_path=strcat(camout_path,camout_file);
camout=csvread(camout_full_path,6,0);
%% Plot
plot(camout(:,1)*10^6,camout(:,2)/50,'Color',[0 0.5 0],'LineWidth',2);
end