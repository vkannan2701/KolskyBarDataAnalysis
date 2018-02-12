% Modified to input the data array instead of the filename as in the
% original code
function [Signalf]=Filter_Kannan_realsignal(Time,Signal,Fs)
clc;
close all;
% %% Read data from file and copy into vectors Signal and Time
% [Data]=ReadDataFile(filename);
% Time=Data{1,1}; %Time vector
% Signal=Data{1,3}; %Incident Pulse
%% Input
% Fs=2.5e9;%Sampling frequency
% Ts=1/Fs; %Inter-sample time
Np_arr=size(Signal);
Np=Np_arr(1); %number of data points
freq=zeros(Np,1);
cutoff=2e6; %cutoff frequency (filter)
%% Plot Signal
figure;
subplot(3,1,1);
plot(Time,Signal);
title('Raw Data (Voltage Vs Time)','FontSize',18,'FontName','Arial Narrow')
xlabel('Time(s)','FontSize',14,'FontName','Arial');
ylabel('Voltage','FontSize',14,'FontName','Arial');
grid on;
%% Fourier Transform using in-built Matlab function
%NFFT=2^nextpow2(Np);
Y=fft(Signal);
for i=1:Np
    freq(i)=(i-1)*(Fs/Np);
end
%% Plot FFT
subplot(3,1,2);
plot(abs(Y(1:Np/2)));
title('FFT output','FontSize',18,'FontName','Arial Narrow');
xlabel('No. of points','FontSize',14,'FontName','Arial');
ylabel('FFT output','FontSize',14,'FontName','Arial');
grid on;

subplot(3,1,3);
plot(freq(1:Np/2),abs(Y(1:Np/2)));
title('Frequency Content','FontSize',18,'FontName','Arial Narrow');
xlabel('Frequency (Hz)','FontSize',14,'FontName','Arial');
ylabel('FFT output','FontSize',14,'FontName','Arial');
grid on;
%% Filter
for i=1:Np/2+1 %There is a catch here! I need to remove the frequency corresponding to the aliased part of my frequency spectrum as well.
    if freq(i)>=cutoff
        Y(i)=0;
    end
end
for i=2:Np/2+1 %There is a catch here! I need to remove the frequency corresponding to the aliased part of my frequency spectrum as well.
    if freq(i)>=cutoff
        Y(Np-i+2)=Y(i);
    end
end
Signalf=ifft(Y);
Yf=fft(Signalf);
%% Plot filtered Signal
figure;
subplot(3,1,1)
plot(Time,Signalf,'r');
title('Filtered Data','FontSize',18,'FontName','Arial Narrow');
xlabel('Time(s)','FontSize',14,'FontName','Arial');
ylabel('Voltage','FontSize',14,'FontName','Arial');
grid on;

subplot(3,1,2);
hold on;
plot(abs(Y(1:Np/2)),'color','blue');
plot(abs(Yf(1:Np/2)),'color','red');
title('FFT output','FontSize',18,'FontName','Arial Narrow');
xlabel('No. of points','FontSize',14,'FontName','Arial');
ylabel('FFT output','FontSize',14,'FontName','Arial');
grid on;
hold off;

subplot(3,1,3);
hold on;
plot(freq(1:Np/2),abs(Y(1:Np/2)));
plot(freq(1:Np/2),abs(Yf(1:Np/2)),'color','red');
title('Frequency Content','FontSize',18,'FontName','Arial Narrow');
xlabel('Frequency (Hz)','FontSize',14,'FontName','Arial');
ylabel('FFT output','FontSize',14,'FontName','Arial');
grid on;
hold off;
end