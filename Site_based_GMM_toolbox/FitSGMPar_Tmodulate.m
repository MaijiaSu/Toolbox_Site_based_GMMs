function  SGMPar = FitSGMPar_Tmodulate(GMi,SGMPar)
% input: GM, given model type of time-modulating function
% output: Estimated Pars for the given model

%% Data Preparation
% Compute the time instant 'GM.t_AI_p' that reaches p% of the total Arial Intensity
AI_t = cumtrapz(GMi.teq,GMi.eq.^2);
AI_n = AI_t/AI_t(end)*100;
[~,ID_repeat] = unique(AI_n);
GMi.AI_p = [0.001,0.05,1:99,99.5,99.99];
GMi.t_AI_p = interp1(AI_n(ID_repeat),GMi.teq(ID_repeat),GMi.AI_p);
% Complement the starting time and the ending time
GMi.AI_p = [0,GMi.AI_p,100];
GMi.AI_t = GMi.AI_p/100*AI_t(end);
GMi.t_AI_p = [0,GMi.t_AI_p,GMi.teq(end)];

% Compte the ID in GM.eq that reaches p% of the total Arial Intensity
% for nn = 1:numel(GMi.AI_p)
%     GMi.ID_p(nn) = find(AI_n>=GMi.AI_p(nn),1,'first');
% end

%% Modulate function#1: Spline Interpolation
if strcmp(SGMPar.Tmodulate.type,'spline') 
    % save the parameter
    Par.type = 'spline';
    Par.p = [0 5 30 45 75 95 100];
    Par.t_AI_p = interp1(GMi.AI_p,GMi.t_AI_p,[0 5 30 45 75 95 100]);
    Par.AI = AI_t(end);
    Par.s_time = SGMPar.Tmodulate.s_time;
end

%% Modulate function#2: Gamma function
if strcmp(SGMPar.Tmodulate.type,'Gamma')
    % Key features of AI(t)
    t_5 = GMi.t_AI_p(find(GMi.AI_p==5));
    t_95 = GMi.t_AI_p(find(GMi.AI_p==95));
    t_45 = GMi.t_AI_p(find(GMi.AI_p==45));
    D5_95  = t_95-t_5;
    tmid = t_45;
      
    % save the parameter
    Par.type = 'Gamma';
    Par.D5_95 = D5_95;
    Par.tmid = tmid;
    Par.AI = AI_t(end);
    Par.s_time = SGMPar.Tmodulate.s_time;
end

%% Modulate function#3: ....


%% save data
SGMPar.Tmodulate.Par = Par;

end


%% ------------------------------------------------------------------------
%                     Defination of sub  function
%% ------------------------------------------------------------------------
function F=opt_q_para(theta,D,t45)
theta_2=theta(1);
theta_3=theta(2);
F(1)=D-gaminv(0.95,2*theta_2-1,theta_3/2)+gaminv(0.05,2*theta_2-1,theta_3/2);
F(2)=t45-gaminv(0.45,2*theta_2-1,theta_3/2);
end