function [x,err,Cp]=NewtonIterativeScheme(w,x0,niter,etol)
%% FUNCTION TO SOLVE THE POCHCHAMMER-CHREE DISPERSION RELATION NUMERICALLY
% THE NUMERICAL SCHEME AS USED IN A JMPS PAPER BY ZHAO AND GARY(1995)
% THE SOLUTION TO THE P-C RELATION IS A COMPLEX VALUE OF WAVENUMBER, K
% RE(k) GIVES US w VS Cp; Im(k): w Vs ATTENUATION COEFFICIENT
% USE A NON-LINEAR 2D OPTIMIZATION PROBLEM/ NEWTON's ITERATIVE PROCEDURE
% k_(n+1)=k_(n)-f(k_(n))/f'(k_(n))
% The fact that we use the first derivative means we assume (x-x0) is small
%% Define input parameters
% omega=w;
% init=x0;
% N=niter;
% tol=etol;
%% Define input parameters for the dispersion relation
% Bar radius, a(m)
a=3.5e-3;
% Longitudinal wave speed,cd(m/s)
cd=5091;
% Shear wave speed, cs(m/s)
cs=3083;
% %Poisson's Ratio, nu
% nu=0.3;
% delta
delta=1500;
%% Define the function and the derivative
syms k;
p=sqrt((w^2/cd^2)-k^2);
q=sqrt((w^2/cs^2)-k^2);
f=(2/a)*p*(w^2/cs^2)*besselj(1,p*a)*besselj(1,q*a)-(4*k^2*p*q)*besselj(0,q*a)*besselj(1,p*a)-((2*k^2-(w^2/cs^2))^2)*besselj(0,p*a)*besselj(1,q*a);
% f=k^2-7*k+1;
df=diff(f,k);
%% Newton's iterative solver
% Initial guess=x0
Cp=zeros(niter,1);
err=zeros(niter,1);
x=zeros(niter,1);
x(1)=x0-(subs(f,x0)/subs(df,x0));
err(1)=x(1)-x0;
i=1;
while i<niter
    x(i+1)=x(i)-(subs(f,x(i))/subs(df,x(i)));
    err(i)=x(i+1)-x(i);
    i=i+1;
end
if err(niter-1)>etol
    disp('No convergence!');
end 
%% Stability of the numerical scheme
stability=subs(f,x(niter-1));
disp(real(stability));
if abs(double(stability))>1
    disp('Possible numerical instability!');
%     w=w+delta;
%     [wn,e,c]=NewtonIterativeScheme(omega,init,N,tol);
%     x=wn;
%     err=e;
%     Cp=c;
end
%% Calculate Cp-phase velocity
for i=1:niter
    Cp(i)=w/x(i);
end
end
    
    
