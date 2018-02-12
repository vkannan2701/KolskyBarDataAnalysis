function [sol]=Rayleigh_Approx()
%% Rayleigh approximation for the dispersion relation for propagation of waves down and infinite bar of known cross-section
%% Input variables
cb=5000;
cs=3000;
nu=0.3;
sigma=3.8317;
a=3.5e-3;
%% Frequancy domain
w=zeros(200,1);
%% Solve the cubic 
sol=zeros(200,3);
for i=1:200
    w(i)=2*pi*i*3e4;
    sol(i,:)=roots([-(nu^2*sigma^2/4)*(cb/cs) 0 (cb/cs) -w(i)*a/(sigma*cs)]);
end
% figure;
% plot(sigma*abs(sol(:,1))/pi,a*w(:)/(pi*cs));
% hold on;
% plot(sigma*abs(sol(:,2))/pi,a*w(:)/(pi*cs));
% hold on;
plot(sigma*abs(sol(:,3))/pi,a*w(:)/(pi*cs));
