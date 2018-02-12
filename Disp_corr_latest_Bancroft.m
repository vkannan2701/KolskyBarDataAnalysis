% Dispersion Correction of linear elastic longitudinal pulses on a Kolsky
% bar apparatus
% Written by, Vignesh Kannan, Ramesh Lab, Latrobe 026
% plot syntax: plot(x,y,'--gs',...
%    'LineWidth',2,...
%    'MarkerSize',10,...
%    'MarkerEdgeColor','b',...
%    'MarkerFaceColor',[0.5,0.5,0.5])
function [Signal_corr,Time]=Disp_corr_latest_Bancroft(Signal,Time,z,Fs)
clc;
%% input parameters
% Fs=100e6;
% Size of time vector
Tnum=size(Time);
Tnumrec=Tnum(1);
% Input parameters
D=12.7e-3;
co=5000;
%% Time vector
Np=Tnumrec;
Np=int64(Np);
Time(:)=Time(:)-Time(1);
figure;
plot(Time,Signal,'-b');
%% Create Signal2
% A=[0;0;0.006;0.006];
% omega=[2*pi/200e-6;2*pi*10e6;2*pi*1e6];
% Signal(:,1)=A(1)*sin(omega(1)*Time(:,1))+A(2)*sin(omega(2)*Time(:,1))+A(3)*sin(omega(3)*Time(:,1))+A(4);
% figure;
% plot(Time,Signal,'-b');
%% Calculate the Fourier Transform
%% Using matlab's fft function
F1=fft(Signal);
size(F1)
%Frequency vector
freq=zeros(Np,1);
for i=1:Np
    freq(i,1)=(i-1)*(Fs/Np);
end
figure;
plot(freq(1:Np),abs(F1(1:Np)),'-ob');
%% Reconstruct the signal
Signal_reconst=ifft(F1(1:Np),'symmetric');
figure;
plot(Time(1:Np),Signal_reconst(1:Np),'-b','LineWidth',2);
title('Reconstructed Signal from modified frequency spectrum','FontSize',18,'FontName','Arial Narrow');
xlabel('Time (s)','FontSize',14,'FontName','Arial');
ylabel('Signal (arb)','FontSize',14,'FontName','Arial');
grid on;
%% Dispersion Relation for a cylindrical bar
%    %% Use Pochammer-Chree solution
%    omega=2*pi*freq;
%    [wn,c_p]=Solve_PochchammerChree_latest(omega);
%% Use Bancroft's data to plot an empirical function of phase velocity as a function of frequency
% Bancroft Data
d=zeros(27,1);
d(1,1)=0;
for i=2:21
    d(i,1)=d(i-1,1)+0.05;
end
d(22:27,1)=[1.2;1.4;1.6;1.8;2.0;10];
Cno=[1;0.99944;0.99774;0.99482;0.99054;0.98466;0.97691;0.96688;0.95410;0.93810;0.91854;0.89549;0.86964;0.84222;0.81466;0.78818;0.76357;0.74125;0.72130;
    0.70365;0.68814;0.64321;0.61687;0.60111;0.59139;0.58524;0.57516];
%% plot data
figure;
plot(d,Cno,'-ob','LineWidth',2);
hold on;
%    plot(D*wn/(2*pi),c_p/co);
title('C_n/C_o Vs d/\lambda','FontSize',18,'FontName','Arial Narrow');
xlabel('d/\lambda','FontSize',14,'FontName','Arial');
ylabel('C_n/C_o','FontSize',14,'FontName','Arial');
legend({'Bancroft(1941)','P-C solution'},'FontSize',18,'FontWeight','bold');
grid on;

% Replot Bancroft's data as Cn/Co Vs d\nu/C_o
dnu=d.*Cno;
%% Curve fit (Empirical function) from Bancroft's data(1941):Cn/Co Vs d\nu/C_o
% Form of equation given by Felice(1986)
syms x y;
y=0.575+0.4189/(63.04*x^4-62.43*x^3+27.76*x^2-7.731*x^1.5+1);
%% plot bancroft's data and the curve fit
figure;
plot(dnu,Cno,'ob','LineWidth',2);
hold on;
ezplot(y,[0,6,0.5,1]);
hold on;
%    size(omega)
%    size(c_p)
%    plot(D*omega/(2*pi*co),c_p/co);
title('(\nu=0.3) : C_n/C_o Vs d\nu/C_o','FontSize',18,'FontName','Arial Narrow');
xlabel('d\nu/C_o','FontSize',18,'FontName','Arial');
ylabel('C_n/C_o','FontSize',18,'FontName','Arial');
legend({'Bancroft(1941)','Curve fit (Felice(1986)'},'FontSize',18,'FontWeight','bold');
grid on;
%% Dispersion Correction
% Calculate phase shift as a function of frequency and apply this to the
% inverse fft
%% Calculate values of Cn for the data
cn=zeros(Np,1);
a=D*freq(1,1)/co;
cn(1,1)=co*(0.575+0.4189/(63.04*a^4-62.43*a^3+27.76*a^2-7.731*a^1.5+1));
for i=2:Np
    if i<=((Np-1)/2)+2
        a=D*freq(i,1)/co;
    else
        a=D*freq(Np-i+1,1)/co;
    end
    cn(i,1)=co*(0.575+0.4189/(63.04*a^4-62.43*a^3+27.76*a^2-7.731*a^1.5+1));
end
figure;
plot(freq,cn,'-b','LineWidth',2);
title('Cn Vs freq','FontSize',18,'FontName','Arial Narrow');
xlabel('\nu (s^{-1})','FontSize',18,'FontName','Arial');
ylabel('Phase velocity: C_n (m/s)','FontSize',18,'FontName','Arial');
grid on;
%% Phase shift: an array of the same size as the fft array
phi=zeros(Np,1);
phi(1)=2*pi*freq(1)*z*(1/co-1/cn(1,1));
for i=2:Np
    if i<=(Np-1)/2+2
        phi(i)=2*pi*freq(i)*z*(1/co-1/cn(i,1));
%         phi(i)=0;
    else
        phi(i)=2*pi*freq(Np-i+1)*z*(1/co-1/cn(i,1));
%         phi(i)=0;
    end
end
figure;
plot(freq,phi,'-b','LineWidth',2);
title('Cn Vs Phase shift','FontSize',18,'FontName','Arial Narrow');
xlabel('\nu (s^{-1})','FontSize',18,'FontName','Arial');
ylabel('Phase shift: \phi (rad)','FontSize',18,'FontName','Arial');
grid on;
%% Plot the phase shift and the phase velocities
figure;
hold on;
[X_1,Y1_1,Y2_1]=plotyy(freq,cn,freq,phi);
title('C_n and \phi as function of \nu','FontSize',18,'FontName','Arial Narrow');
xlabel('freq: \nu (s^{-1})','FontSize',14,'FontName','Arial');
ylabel(X_1(1),'Phase Velocity (m/s) ','FontSize',14,'FontName','Arial');
ylabel(X_1(2),'Phase shift (rad)','FontSize',14,'FontName','Arial');
set(Y1_1,'linewidth',0.5,'color','blue','LineWidth',2);
set(Y2_1,'linewidth',0.5,'color','red','LineWidth',2);
grid on;
hold off;
%% Dispersion Correction
%% Correction is done for the fourier coefficients of the transform
% Matlab's fft: X(k)=sum_j(x(j)exp(-i 2pi/N (j-1) (k-1)))
% x=sum(X(k)exp(-i 2pi freq_n t))
% adding a phase shift to the arguement of the above function,
% x_new=sum(X_new(k) exp(-i 2pi freq_n t)); X_new(k)=X(cos(phi(k))+isin(phi(k)))
F_new=zeros(Np,1);
size(phi)
size(F1)
for i=1:Np
    if freq(i)<=2e6
        F_new(i,1)=F1(i,1)*complex((cos(phi(i))),(sin(phi(i))));
    end
%     check(i)=phi(i)+2*pi*freq(i)*Time(i);
end
% F_new(((Np-1)/2)+2)=F_new(((Np-1)/2)+1);
figure;
plot(freq,real(F1),'-^b','LineWidth',0.5);
hold on;
plot(freq,real(F_new),'-xr','LineWidth',0.5);
title('Original and corrected fourier coefficients','FontSize',18,'FontName','Arial Narrow')
xlabel('freq: \nu (s^{-1})','FontSize',14,'FontName','Arial');
ylabel('Fourier Coefficients','FontSize',14,'FontName','Arial');
hold off;

figure;
plot(freq,imag(F1),'-^b','LineWidth',0.5);
hold on;
plot(freq,imag(F_new),'-xr','LineWidth',0.5);
title('Original and corrected fourier coefficients','FontSize',18,'FontName','Arial Narrow')
xlabel('freq: \nu (s^{-1})','FontSize',14,'FontName','Arial');
ylabel('Fourier Coefficients','FontSize',14,'FontName','Arial');
hold off;
%% Reconstruct the corrected signal
Signal_corr=ifft(F_new(1:Np),'symmetric');
figure;
plot(Time(1:Np),Signal_reconst(1:Np),'-b','LineWidth',0.2);
hold on;
plot(Time(1:Np),Signal_corr(1:Np),'-r','LineWidth',0.2);
hold off;
end