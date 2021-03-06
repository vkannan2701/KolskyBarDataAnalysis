function[wavenumber,phasev]=Solve_PochchammerChree_latest(omega)
clc;
close all;
%% Define input parameters for the dispersion relation
% Bar radius, a(m)
a=3.5e-3;
% Longitudinal wave speed,cd(m/s)
cd=5091;
% Shear wave speed, cs(m/s)
cs=3083;
% %Poisson's Ratio, nu
nu=0.3;
sigma=3.8317;
% Number of iterations
N=50;
%size of frequency vector
numberoffreq=size(omega)
Nw=numberoffreq(1);
% number of initial guesses
nsol=1;
%% Input the frequency spectrum
% omega=zeros(Nw,1);
% for i=1:Nw
%     omega(i)=2*pi*i*0.5e4;
% %     omega(i)=1;
% end
%% Solve for the full range of omega
init_guess=zeros(Nw,1);
sol=zeros(Nw,3);
for i=1:Nw
    %Use the Rayleigh approximation for init guess
    sol(i,:)=real(roots([-(nu^2*sigma^2/4)*(cd/cs) 0 (cd/cs) -omega(i)*a/(sigma*cs)]));
    for j=1:nsol
        init_guess(i,j)=sigma*abs(sol(i,3))/a;
    end
end
% Solve and rerun loop if the solution does not converge
i=1;
while i<=Nw
%     count=1;
    for j=1:nsol
        [K{i,j},~,C_P{i,j}]=NewtonIterativeScheme(omega(i),init_guess(i,j),N,1e-4);
        p=sqrt((omega(i)^2/cd^2)-K{i,j}(N,1)^2);
        q=sqrt((omega(i)^2/cs^2)-K{i,j}(N,1)^2);
        F(i,j)=(2/a)*p*(omega(i)^2/cs^2)*besselj(1,p*a)*besselj(1,q*a)-(4*K{i,j}(N,1)^2*p*q)*besselj(0,q*a)*besselj(1,p*a)-((2*K{i,j}(N,1)^2-(omega(i)^2/cs^2))^2)*besselj(0,p*a)*besselj(1,q*a);
%         if F(i,j)>1
%             if omega(i)<2*pi*2e6
%             disp('Loop re-run');
%             init_guess(i,j)=init_guess(i,j)+(0.5*pi/(a));
%             i=i-1;
% %             count=count+1;
%             end
%         end
    end
    i=i+1;
end       
%% Plot the solutions
figure;
hold on;
xlabel('ak/\pi','FontSize',14,'FontWeight','bold','Color','b');
ylabel('\omega a/\pi c_t','FontSize',14,'FontWeight','bold','Color','b');
for i=1:Nw
    for j=1:nsol
        plot((3.5e-3)*K{i,j}(N,1)/(pi),(omega(i)*(3.5e-3)/(pi*3000)),'o','MarkerEdgeColor','b','MarkerFaceColor','b');
        hold on;
        plot(a*init_guess(i,j)/(pi),(omega(i)*(3.5e-3)/(pi*3000)),'o','MarkerEdgeColor','r','MarkerFaceColor','r');
    end
end
[filename,path]=uigetfile();
filename=strcat(path,filename);
achenbach_data=xlsread(filename);
plot(achenbach_data(:,1),achenbach_data(:,2),'d','MarkerEdgeColor','g','MarkerFaceColor','g');
% figure;
% for i=1:1
%     for j=1:10
%         plot((7e-3)*omega(i)/(2*pi*5000),C_P{i,j}(20,1)/5000,'o');
%         hold on;
%     end
% end
% figure;
% for i=1:1
%     for j=1:10
%         plot((3.5e-3)*K{i,j}(20,1)/(2*pi),C_P{i,j}(20,1)/5000,'o');
%         hold on;
%     end
% end
%% Replot the P-C equation using the numerical solution obtained
%% Define the function and the derivative
% f=k^2-7*k+1;
% df=diff(f,k);
%% Plot
for i=1:Nw
    for j=1:nsol
        p=sqrt((omega(i)^2/cd^2)-K{i,j}(N,1)^2);
        q=sqrt((omega(i)^2/cs^2)-K{i,j}(N,1)^2);
        F(i,j)=(2/a)*p*(omega(i)^2/cs^2)*besselj(1,p*a)*besselj(1,q*a)-(4*K{i,j}(N,1)^2*p*q)*besselj(0,q*a)*besselj(1,p*a)-((2*K{i,j}(N,1)^2-(omega(i)^2/cs^2))^2)*besselj(0,p*a)*besselj(1,q*a);    
    end
end
figure;
for j=1:nsol
    hold on;
    plot(omega(:)*(3.5e-3)/(pi*3000),real(F(:,j)),'-o');
    xlabel('\omega a/\pi c_t','FontSize',14,'FontWeight','bold','Color','b');
    ylabel('F(\omega,k)- Value of the P-C Equation','FontSize',14,'FontWeight','bold','Color','b');
end
%% Save output variables
wavenumber=zeros(Nw,1);
phasev=zeros(Nw,1);
for i=1:Nw
    wavenumber(i)=K{i,nsol}(N-1);
    phasev(i)=C_P{i,nsol}(N-1);
end