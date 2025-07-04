function [t_q_pad,qmod,padn] = TModulate_Estimator(Par,dt,L)
%% 
% Input: 
% Output:#1 A discrete temporal modulation,i.e., {(ti,qmod_i),i=1,...,N}
%        #2 Length of zeros padding

% DefaultPar = inputParser;
% DefaultPar.KeepUnmatched = 1;
% DefaultPar.addParameter('type','spline');
% parse(DefaultPar,Par);
% Par = DefaultPar.Results;
% if ~isfield(Par,'type')
%     Par.type = 'spline';
% end

if ~isfield(Par,'s_time')
    Par.s_time = 0;
end

%% ------------------------------------------------------------------------
%                   Method#1： Spline interpolation
%% ------------------------------------------------------------------------
if strcmp(Par.type,'spline')
    %% Par
    AI_p = Par.AI(end)*Par.p;
    t_AI_p = Par.t_AI_p;
    if ~all(diff(t_AI_p)>0)
        disp('Pars of Spline function are not increased monotonicaly')
    end

    s_time = Par.s_time;
    teq = (0:L-1)*dt;

    % Hann window
    Lt1 = ceil((s_time/dt)/2)*2+1;
    wt = hann(Lt1)/sum(hann(Lt1));  % weight time avaraging
    padn = (Lt1-1)/2;

    % Hermite Interpolation
    AI_t = pchip(t_AI_p,AI_p,teq);
    AI_t = AI_t/AI_t(end);
    AI_t_pad = [zeros(1,padn*2) AI_t ones(1,padn*2)]; % zeros-padding
    L_pad = L+padn*2;

    % Smooth AI(t) with hann window
    AI_t_smooth = zeros(1,L_pad);
    for k = 1:L_pad
        kk = padn+k;
        AI_t_smooth(k) = AI_t_pad(kk-padn:kk+padn)*wt;
    end

    % finite difference
    AI_t_smooth = AI_t_smooth/AI_t_smooth(end);
    q_2_t = [0 diff(AI_t_smooth)/dt];
    t_q_pad = ((1:numel(q_2_t))-1)*(dt);
    q_2 = q_2_t/trapz(t_q_pad,q_2_t)*Par.AI(end);

    % Time-modulation function
    qmod = sqrt(q_2);

end

%% ------------------------------------------------------------------------
%                     Method#2： Gamma function
%% ------------------------------------------------------------------------
if strcmp(Par.type,'Gamma')
   %% Input
   %  q(t)=theta_1*t.^(theta_2-1).*exp(-t/theta_3)
   AI      = Par.AI;
   D5_95  = Par.D5_95;
   tmid = Par.tmid;
   if L*dt<tmid
       error('Wrong Gamma modulating Pars: tf<tmid ')
   end
   theta_2 =(1/0.3)^2;
   theta_3 = tmid/((1/0.3)^2);
   theta_0 = [theta_2, theta_3];
   % Perfrom optimization
   f_option.Display = 'off';
   theta_t = fsolve(@(x) opt_q_para(x,D5_95,tmid),theta_0,f_option);
%    theta_t = fsolve(@opt_q_para,theta_0,[],D5_95,tmid);
   % update the initial guess
   theta_2 = theta_t(1);
   theta_3 = theta_t(2);

   % other Par   
   teq = (0:L-1)*dt;

   % Modulating function before smooth
   q_t = gampdf(teq,2*theta_2-1,theta_3/2);    %.*filt;

   % With Padding and Smoothness
   s_time = Par.s_time;
   Lt1    = ceil((s_time/dt)/2)*2+1;
   wt     = hann(Lt1)/sum(hann(Lt1));  % weight time avaraging
   yy     = [0 cumsum(q_t(2:end))]/sum(sum(q_t(2:end)));
   padn=(numel(wt)-1)/2;
   yy_pad = [zeros(1,padn*2) yy ones(1,padn*2)];
   yyf_t = zeros(1,numel(yy));
   for k = 1:(numel(yy)+padn*2)
       kk = padn+k;
       yyf_t(k) = yy_pad(kk-padn:kk+padn)*wt;
   end
   %yyf = yyf_t(padn:end-padn);
   yyf = yyf_t;
   yyf = yyf/yyf(end);
   q_2_t = [0 diff(yyf)/dt];
   t_q_pad = ((1:numel(q_2_t))-1)*(dt);
   q_2 = q_2_t/trapz(t_q_pad,q_2_t)*AI;

   % Time-modulation function
    qmod = sqrt(q_2);
end



end


%% ------------------------------------------------------------------------
%                     Defination of sub  function
%% ------------------------------------------------------------------------
function F=opt_q_para(theta,D5_95,t45)
theta_2=theta(1);
theta_3=theta(2);
F(1)=D5_95-gaminv(0.95,2*theta_2-1,theta_3/2)+gaminv(0.05,2*theta_2-1,theta_3/2);
F(2)=t45-gaminv(0.45,2*theta_2-1,theta_3/2);
end
