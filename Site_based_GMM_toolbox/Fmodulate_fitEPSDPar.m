function Filter_opt_par = Fmodulate_fitEPSDPar(Paramtericmodel,...
    FittedFilterPar,GM,Tmodulate)
%% Use a Paramteric model to catch the trend of extimated pars (BN) of analytical EPSD
% Input: i. Type of Parameteric models for EPSD Pars
%        ii. Estimated values of BN
%        iii. Other Par
% Output: Parameters of Parameteric models for analytical EPSD model Pars

%% Data Preparation
% Compute the time instant 'GM.t_AI_p' that reaches p% of the total Arial Intensity
AI = cumtrapz(GM.teq,GM.eq.^2);
AI_n = AI/AI(end)*100;
[~,ID_repeat] = unique(AI_n);
GM.AI_p = [0.001,0.05,1:99,99.5,99.99];
GM.t_AI_p = interp1(AI_n(ID_repeat),GM.teq(ID_repeat),GM.AI_p);
% Complement the starting time and the ending time
GM.AI_p = [0,GM.AI_p,100];
GM.AI = GM.AI_p/100*AI(end);
GM.t_AI_p = [0,GM.t_AI_p,GM.teq(end)];
% Compte the ID in GM.eq that reaches p% of the total Arial Intensity
for nn = 1:numel(GM.AI_p)
    GM.ID_p(nn) = find(AI_n>=GM.AI_p(nn),1,'first');
end

% AI_p = [0 5 30 45 75 95 100];
% t_AI_p = interp1(GM.AI_p,GM.t_AI_p,[0 5 30 45 75 95 100]);
ID_p = floor(interp1(GM.AI_p,GM.ID_p,[0 5 30 45 75 95 100]));
% AI = interp1(GM.AI_p,GM.AI,[0 5 30 45 75 95 100]);
% dt = GM.dt;
FilterNofPar = size(FittedFilterPar,2);

%% Define the functional type of Trend
% TrendType(nn) = 1,'costant'; = 2, 'linear';=3,'polyline'
if strcmp(Paramtericmodel,'Constantmodel')
    TrendType = ones(1,FilterNofPar)*1;
elseif strcmp(Paramtericmodel,'Linearmodel')
    TrendType = ones(1,FilterNofPar)*2;
elseif strcmp(Paramtericmodel,'ReducedLinearmodel')
    if FilterNofPar == 3
        TrendType = [2,1,nan];
    elseif FilterNofPar == 5
        TrendType = [2,1,2,1];
    elseif FilterNofPar == 6
        TrendType = [2,1,2,1,2];
    end
elseif strcmp(Paramtericmodel,'Polymodel')
    TrendType = ones(1,FilterNofPar)*3;
else 
    TrendType = Paramtericmodel;
end

%% Loop to estimate pars of Trend
for nn = 1:FilterNofPar-1
    switch TrendType(nn)

        case 1 % Costant (One Par)
            % % Prepare the linear fitting
            % compute the ID
            id_tmid = find(GM.AI_p==45);      % ID of time reaching 45% AI
            id_t5 = GM.ID_p(find(GM.AI_p==5));  % ID of time reaching 5% AI
            id_t95 = GM.ID_p(find(GM.AI_p==95)); % ID of time reaching 95% AI
            % compute the weighting coefficient
            [~,qmod,~] = TModulate_Estimator(Tmodulate.Par,GM.dt,GM.L);
            w_coef = qmod(id_t5:id_t95)'/sum(qmod(id_t5:id_t95));
            % construct the function
            tmid = GM.t_AI_p(id_tmid);
            modelfun_weight = @(b,t) sqrt(interp1(GM.teq(id_t5:id_t95),w_coef,t))...
                .*(b(1)+b(2)*(t-tmid));
            modelfun = @(b,t) (b(1)+b(2)*(t-tmid));

            % % Start to fit
            Neq_t = numel(GM.teq(ID_p(2):ID_p(end-1)));
            BN_smoth= smooth(FittedFilterPar(:,nn),floor(((Neq_t)/16)*2));
            % starting points
            beta0 = [BN_smoth(id_tmid),(BN_smoth(id_t5)-BN_smoth(id_t95))/(GM.teq(id_t95)-GM.teq(id_t5))];
            % data with weighting
            x_data = GM.teq(id_t5:id_t95)';
            y_data = sqrt(w_coef).*BN_smoth(id_t5:id_t95);
            % Peform weighting least square fitting
            options = optimset('Display','off');
            bn = lsqcurvefit(modelfun_weight,beta0,x_data,y_data,[],[],options);
            %         BW0 = modelfun(bn,GM.teq(id_t5:id_t95));
            %         BW0 = [ones(1,id_t5-1)*BW0(1),BW0,ones(1,GM.L-id_t95)*BW0(end)];

            %         BW(:,par_no) = BW0;
            Filter_opt_par(:,nn) = bn(1);
            Filter_opt_par(2,nn) = 0;

        case 2 % LinearModel (Two Pars)
            % % Prepare the linear fitting
            % compute the ID
            id_tmid = find(GM.AI_p==45);      % ID of time reaching 45% AI
            id_t5 = GM.ID_p(find(GM.AI_p==5));  % ID of time reaching 5% AI
            id_t95 = GM.ID_p(find(GM.AI_p==95)); % ID of time reaching 95% AI
            % compute the weighting coefficient
            [~,qmod,~] = TModulate_Estimator(Tmodulate.Par,GM.dt,GM.L);
            w_coef = qmod(id_t5:id_t95)'/sum(qmod(id_t5:id_t95));
            % construct the function
            tmid = GM.t_AI_p(id_tmid);
            modelfun_weight = @(b,t) sqrt(interp1(GM.teq(id_t5:id_t95),w_coef,t))...
                .*(b(1)+b(2)*(t-tmid));
            modelfun = @(b,t) (b(1)+b(2)*(t-tmid));

            % % Start to fit
            Neq_t = numel(GM.teq(ID_p(2):ID_p(end-1)));
            BN_smoth= smooth(FittedFilterPar(:,nn),floor(((Neq_t)/16)*2));
            % starting points
            beta0 = [BN_smoth(id_tmid),(BN_smoth(id_t5)-BN_smoth(id_t95))/(GM.teq(id_t95)-GM.teq(id_t5))];
            % data with weighting
            x_data = GM.teq(id_t5:id_t95)';
            y_data = sqrt(w_coef).*BN_smoth(id_t5:id_t95);
            % Peform weighting least square fitting
            options = optimset('Display','off');
            bn = lsqcurvefit(modelfun_weight,beta0,x_data,y_data,[],[],options);
            %         BW0 = modelfun(bn,GM.teq(id_t5:id_t95));
            %         BW0 = [ones(1,id_t5-1)*BW0(1),BW0,ones(1,GM.L-id_t95)*BW0(end)];

            %         BW(:,par_no) = BW0;
            Filter_opt_par(:,nn) = bn';

        case 3 % Polyline (Five Pars)
            Neq_t = numel(GM.teq(ID_p(2):ID_p(end-1)));
            y = smooth(FittedFilterPar(:,nn),floor(((Neq_t)/16)*2));
            % smooth the result for about 1/16 of the legth
            pi_5      = y(ID_p(2));
            pi_30     = y(ID_p(3));
            pi_45     = y(ID_p(4));
            pi_75     = y(ID_p(5));
            pi_95     = y(ID_p(6));
            PI_val_t  = [pi_5 pi_30 pi_45 pi_75 pi_95];
            %
            Filter_opt_par(:,nn)   = PI_val_t';      
    end
end

%%  polyline model
% if strcmp(Paramtericmodel,'Polymodel')
%     Filter_opt_par = zeros(5,FilterNofPar);
%     Neq_t = numel(GM.teq(ID_p(2):ID_p(end-1)));
%     x = GM.teq(ID_p(2):ID_p(end-1))';
%     X = [ones(size(x)) x];
%     for i = 1:FilterNofPar
%         y = smooth(FittedFilterPar(:,i),floor(((Neq_t)/16)*2));
%         smooth the result for about 1/16 of the legth
%         pi_5      = y(ID_p(2));
%         pi_30     = y(ID_p(3));
%         pi_45     = y(ID_p(4));
%         pi_75     = y(ID_p(5));
%         pi_95     = y(ID_p(6));
%         PI_val_t  = [pi_5 pi_30 pi_45 pi_75 pi_95];
% 
%         Filter_opt_par(:,i)   = PI_val_t';
%     end
%     wg_p = opt_par(:,1);
%     zeta_p = opt_par(:,2);
%     Par.wg_p = Filter_opt_par(:,1);
%     Par.zeta_p = Filter_opt_par(:,2);
% end

%%  Linear model
% if strcmp(Paramtericmodel,'Linearmodel')||strcmp(Paramtericmodel,'ReducedLinearmodel')||strcmp(Paramtericmodel,'Constantmodel')
%     Filter_opt_par = zeros(2,FilterNofPar);
%     % compute the ID
%     id_tmid = find(GM.AI_p==45);      % ID of time reaching 45% AI
%     id_t5 = GM.ID_p(find(GM.AI_p==5));  % ID of time reaching 5% AI
%     id_t95 = GM.ID_p(find(GM.AI_p==95)); % ID of time reaching 95% AI
%     % compute the weighting coefficient
%     [~,qmod,~] = TModulate_Estimator(Tmodulate.Par,GM.dt,GM.L);
%     w_coef = qmod(id_t5:id_t95)'/sum(qmod(id_t5:id_t95));
%     % construct the function
%     tmid = GM.t_AI_p(id_tmid);
%     modelfun_weight = @(b,t) sqrt(interp1(GM.teq(id_t5:id_t95),w_coef,t))...
%         .*(b(1)+b(2)*(t-tmid));
%     modelfun = @(b,t) (b(1)+b(2)*(t-tmid));
% 
%     for par_no = 1:FilterNofPar % loops for seperatly optimizing wn and zeta
%         Neq_t = numel(GM.teq(ID_p(2):ID_p(end-1)));
%         BN_smoth(:,par_no) = smooth(FittedFilterPar(:,par_no),floor(((Neq_t)/16)*2));
%         % starting points
%         beta0 = [BN_smoth(id_tmid,par_no),(BN_smoth(id_t5,par_no)-BN_smoth(id_t95,par_no))/(GM.teq(id_t95)-GM.teq(id_t5))];
%         % data with weighting
%         x_data = GM.teq(id_t5:id_t95)';
%         y_data = sqrt(w_coef).*BN_smoth(id_t5:id_t95,par_no);
%         % Peform weighting least square fitting
%         options = optimset('Display','off');
%         bn = lsqcurvefit(modelfun_weight,beta0,x_data,y_data,[],[],options);
%         %         BW0 = modelfun(bn,GM.teq(id_t5:id_t95));
%         %         BW0 = [ones(1,id_t5-1)*BW0(1),BW0,ones(1,GM.L-id_t95)*BW0(end)];
% 
%         %         BW(:,par_no) = BW0;
%         Filter_opt_par(:,par_no) = bn';
%     end
% 
%     % Linear model Par
%     Par.wg_tmid   = Filter_opt_par(1,1);
%     Par.wg_rate   = Filter_opt_par(2,1);
%     Par.zeta_tmid = Filter_opt_par(1,2);
%     Par.zeta_rate = Filter_opt_par(2,2);
%     %     Par.tmid = tmid;
%     %     Par.id_t5 = id_t5;
%     %     Par.id_t95 = id_t95;
% end

%%  ReducedLinearmodel model
% % % % % A linear model for filter frquency and a constant model for bandwidth
% % % % if strcmp(Paramtericmodel,'ReducedLinearmodel')
% % % %     % Linear model Par
% % % %     Par.wg_tmid   = Filter_opt_par(1,1);
% % % %     Par.wg_rate   = Filter_opt_par(2,1);
% % % %     Par.zeta_tmid = Filter_opt_par(1,2);
% % % %     Par.zeta_rate = 0;Filter_opt_par(2,2)=0;
% % % % 
% % % %     % Linear model Par
% % % %     if FilterNofPar == 6
% % % %         Filter_opt_par(2,2)=0;
% % % %         Filter_opt_par(2,4)=0;
% % % %         %         Filter_opt_par(2,5)=0;
% % % %     elseif  FilterNofPar == 5
% % % %         Filter_opt_par(2,2)=0;
% % % %         Filter_opt_par(2,4)=0;
% % % %     end
% % % % 
% % % % end

%% Constant model
% % % if strcmp(Paramtericmodel,'Constantmodel')
% % %     % Linear model Par
% % %     Par.wg_tmid   = Filter_opt_par(1,1);
% % %     Par.wg_rate   = 0;Filter_opt_par(2,1)=0;
% % %     Par.zeta_tmid = Filter_opt_par(1,2);
% % %     Par.zeta_rate = 0;Filter_opt_par(2,2)=0;
% % % 
% % %     % Linear model Par
% % %     if FilterNofPar == 6
% % %         Filter_opt_par(2,1)=0;
% % %         Filter_opt_par(2,2)=0;
% % %         Filter_opt_par(2,3)=0;
% % %         Filter_opt_par(2,4)=0;
% % %         Filter_opt_par(2,5)=0;
% % %     elseif  FilterNofPar == 5
% % %         Filter_opt_par(2,1)=0;
% % %         Filter_opt_par(2,2)=0;
% % %         Filter_opt_par(2,3)=0;
% % %         Filter_opt_par(2,4)=0;
% % %     end
% % % end

%% Discard
% if strcmp(Paramtericmodel,'NonParameter')
%     BW1 = BN(:,1);
%     BW2 = BN(:,2);
%     % Bound
%     BW1(BW1<=25/500*2*pi) = 25/500*2*pi;
%     BW2(BW2<=0.02)=0.02;
%     BW2(BW2>=0.98)=0.98;
%
%     % Smooth the parameters
%     s_time = 2;
%     wg_ti = smooth(BW1,s_time/Par.dt);
%     zeta_ti = smooth(BW2,s_time/Par.dt);
%
% end