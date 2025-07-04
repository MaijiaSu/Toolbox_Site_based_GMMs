function [ Sd, Sv, Sa,maxT,u] = ResponseSpectraNM(acc,dt,g,T,zeta,u0,du0,gamma,beta)
% This function calculates the Response spectrum using Newmark's method 

%% Input 
%   acc= Acceleration vector Unit [g] 
%   dt = Time step
%    g = gravity
%    T = Period vector
% zeta = Damping ratio unit [absolute], therefore 5% is equivalent to 0.05
%   u0 = Initial Displacement  typically [0]
%  du0 = Initial Velocity      typically [0]

%% Output 
% Sd displacement unit [cm]
% Sv pseudo-velocity   [cm*rad/s]
% Sa pseudo-acc        [cm*rad^2/s^2]

Sd=zeros(length(zeta),length(T));
Sv=zeros(length(zeta),length(T));
Sa=zeros(length(zeta),length(T));
m=1;
for i=1:length(zeta)
    for j=1:length(T)
        wn=2*pi/T(j); %Natural frequency
        k=wn^2; %Stiffness
        u(j,:) = NM(k,m,zeta(i),-acc*g',dt,u0,du0,gamma,beta);
        [Sd(i,j),id] = max(abs(u(j,:)));
        Sd(i,j) = Sd(i,j)*100;
        maxT(i,j) = (id-1)*dt;
%         if Sd(j) > 10000 % too large means discovegence
%            Sd(j) = 0;
%            Sa(j) = max(abs(acc));
%         end
        Sv(i,j)=Sd(i,j)*wn;
        Sa(i,j)=Sd(i,j)*wn^2/g/100;
    end
end
end

