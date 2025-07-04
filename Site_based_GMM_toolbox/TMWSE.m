 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%______________________________ETH_____________________________________________
%_____________________________B E R K E L E Y__________________________________
%______________Simulation of Ground Motion in Freqeuncy Domain ________________
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%______________________________________________________________________________
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function to solve the estimate a EPSD with te  Short-Time Thompson's multiple-window spectrum  
% Paper of reference: Fully Nonstationary Analytical Earthquake Ground-Motion Model 
% Authiors: P. Conte and F. Peng 

%Function wrote written by Marco Broccardo 
function [PHI, PHIn, sgt_2, S,w]=TMWSE(eq,EP_T)
% where
% PHI（t_i）: the estimaed EPSD sacled by sgt_2 at t_i
% PHIn(t_i): EPSD with nomalization that sacles to unit-variacne at t_i
% sgt_2 (t_i): the varince of the weight-signal by time-window at t_i
% wk:the discret point at frequency domain
% S: S(t_i,1:TMWSEPar.L) is the local time series centered at t_i
% w: the time-moving windons, a unit-integral hann-model is used
% L: the length of the Hann window

%% Choose Hanning window 

L=EP_T.L;   % Number of discrete point for the windows
Lt=EP_T.Lt; % Number of discrete point for time avaraging
K=EP_T.K;   % Number of leakage-resistant orthogonal windows 
dt=EP_T.dt; % Time integration step 

%% Preallocation 
S(1:numel(eq),1:L)=0;       % local time series centered at time t_i 
Sf(1:L,1:K)=zeros(L,K);     % Fourier transform of S
Sh(1:(L+1)/2)=0;            % Temporal mediated Fourier transform
SH(1:numel(eq),(L+1)/2)=0;  % Temporal allocation for the unormalized SH

W=(K+1)/(2*L);              % Bandwith of the discrete-time Fourier transforms of the first K DPSSs
w=hann(L)/norm(hann(L));    % weighth variance  
wt=hann(Lt)/sum(hann(Lt));  % weight time avaraging 
sgt_2(1:numel(eq))=0;       % estimated variance of the process 


%% Compute the empirical EPSD
% 1. Construct Leakage-resistance orthogonal windows, i.e.,K Discte prolate
% spherodial sequences (DPSS)
% 1.1. Define the toeplitz matrix, i.e.,T
pT=sin(2.*pi.*W.*((1:L)-1))./(pi*((1:L)-1)); %construction of the Toeplix matrix 
pT(1)=2*W;
T=toeplitz(pT);
% 1.2 Perform eigvetor decompostion to T
[V,D]=eig(T); % eigenvalues
% where V(:,L+1-k) is the kth of the DPSS
% and D(L+1-k,L+1-k) is the corresponding eigenvaule

sh=(L-1)/2;                     % Shift for padding the time series
sht=(Lt-1)/2;                   % Shift for padding the time series
eq_pad=padarray(eq,[sh 0]);     % padding add zero 


% 2. Compute the empirical EPSD
for i=1:numel(eq)
    % 2.1. Compute the K mediated Fourier transform, Sf
    % Peform DFT to the product of local-time series and k-th DPSS
    % the local time series,i,e, S(i,1:L) at t_i, is defines by a
    % time-moving window w.
    S(i,1:L)=w.*eq_pad(i:((i-1)+L));
    Ss=S(i,1:L)';
    sgt_2(i)=norm(Ss)^2;
    for k=1:K
        Sf(:,k)=abs( fft(Ss.*(V(:,L+1-k))) ).^2;     
    end

    % 2.2 A recursive procedure to compute the final empirical EPSD
    do=zeros(K,1);
    for j=1:((L+1)/2)
        dn=ones(K,1);
        while abs(dn-do)>10^-8
            dnn=abs(dn)/norm(dn);
            Sht=dnn'*Sf(j,1:K)';  %equation 27
            %Sht=(abs(dn).^2)'*Sf(j,1:K)'/norm(abs(dn))^2;
            do=dn;
            for k=1:K
                lm=D(L+1-k,L+1-k);
                dn(k)=sqrt(lm)*Sht/(lm*Sht+(1-lm)*sgt_2(i));
            end
        end
        Sh(j)=(abs(dn).^2)'*Sf(j,1:K)'/norm(abs(dn))^2;
    end
    SH(i,:)=Sh;
end


%% Time avareging 
SHt=padarray(SH',[0 sht]);
WT=padarray(wt',[0 (size(SHt,2)-numel(wt))],'post');
WTT=gallery('circul',WT)';
PHItt=SHt*WTT;
PHIt=PHItt(:,(sht+1):(end-sht));

%% Normalize according the area 
Ar=2*trapz(PHIt(:,:))*((2*pi)/(L*dt)); 
Arn=Ar+10^-8;     % spurious value to avoid inf %
PHIn=PHIt.*(ones((L+1)/2,1)*(1./Arn));
PHI=PHIn.*(ones((L+1)/2,1)*(sgt_2));

end 




