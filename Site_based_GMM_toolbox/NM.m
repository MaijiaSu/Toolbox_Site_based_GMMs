function [ u, du, ddu ] = NM(k,m,zeta,p,dt,u0,du0,gamma,beta)
% This funtion uses Newmark's method to calculate displacement,
% velocity and acceleration at each time step 
%   k,m,zeta: Properties of the system: Stiffness, mass and damping ratio
%   p,dt: Loading conditions: load and time step
%   u0,du0: Initial Conditions: Initial Displacement and Initial Velocity
%   gamma, beta: parameters that define variation of acceleration
%   du,ddu: Velocity and Acceleration 

% Initialization
u(1)=u0;
du(1)=du0;
wn=sqrt(k/m); %Natural frequency
c=2*m*wn*zeta; %Damping coefficient
ddu(1)=(p(1)-c*du(1)-k*u(1))/m;
a1=(1/(beta*dt*dt))*m+(gamma/(beta*dt))*c;
a2=m/(beta*dt)+(gamma/beta-1)*c;
a3=(1/(2*beta)-1)*m+dt*(gamma/(2*beta)-1);
kh=k+a1;
P=zeros(size(p)); %Preallocation
u=zeros(size(p)); 
du=zeros(size(p));
ddu=zeros(size(p));

%Calculation for each step
for i=1:length(p)-1
    P(i+1)=p(i+1)+a1*u(i)+a2*du(i)+a3*ddu(i);
    u(i+1)=P(i+1)/kh;
    du(i+1)=(gamma/(beta*dt))*(u(i+1)-u(i))+(1-gamma/beta)*du(i)+dt*(1-gamma/(2*beta))*ddu(i);
    ddu(i+1)=(1/(beta*dt^2))*(u(i+1)-u(i))-(1/(dt*beta)*du(i))-(1/(2*beta)-1)*ddu(i);
end
end


