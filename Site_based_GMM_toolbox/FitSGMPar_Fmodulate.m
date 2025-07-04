function [SGMPar,FitResult] = FitSGMPar_Fmodulate(GM,SGMPar,eEPSD)
%% ------------------------------------------------------------------------ 
%             2-I. Compute Spectral Paramters at each time ti
%% ------------------------------------------------------------------------

%% 2.0 Define a filter model
[GMFilter,FilterNofPar] = SetFilterModel(SGMPar.Fmodulate.type);

%% 2.1 Pad the data for Hann-window smoothing
% Unit-variance empirical EPSD 
PHIn = eEPSD.PHIn;
% Hann window
Time_smoothing = 3; % Equation 11 in ICOSSAR PAPER 
Lt = ceil((Time_smoothing/GM.dt)/2)*2+1; % Number of discrete point for time avaraging
wt = hann(Lt)/sum(hann(Lt));            
WT = ones(size(PHIn,1),1)*wt';
% Padding 
mid_p = (Lt-1)/2;
PHIn_ext = zeros(size(PHIn,1),size(PHIn,2)+(Lt-1));
PHIn_ext(:,1:mid_p) = PHIn(:,1)*ones(1,mid_p);
PHIn_ext(:,(mid_p+1):size(PHIn,2)+mid_p) = PHIn;
PHIn_ext(:,size(PHIn,2)+mid_p+1:end) = PHIn(:,end)*ones(1,mid_p);
 
%% 2.2 Start a loop for optimizing the pars
% Determine the fiting range
id1 = find(GM.AI_p==5);
id2 = find(GM.AI_p==95);
ind_tp1 = GM.ID_p(id1);         % index of time at p1% of Arias Intensity
ind_tp2 = GM.ID_p(id2); % index of time at p2% of Arias Intensity

FittedFilterPar = [];
Fitscalse = 'xylinear'; %'xylinear','ylog'

%% 2.2 Choice I -- linear Sacle 
if strcmp(Fitscalse,'xylinear')
    for i = ind_tp1:ind_tp2
        % ID of PHIn_ext
        ii = i+mid_p;
        % 2.2.1 Compute Wighting empirical PSD with a hann-window averaging at ti
        ePSD_ti = sum(PHIn_ext(:,(ii-mid_p):(ii+mid_p)).*WT,2);

        if sum(ePSD_ti)>0
            % 2.2.2 Determine the range to perform the fitting, i.e., optimize
            % up to 99% of the cumulative specrtum density
            X = eEPSD.wk;
            y = ePSD_ti/(trapz(X,ePSD_ti));
            y_cum = cumtrapz(X,y);
            [~, in_op]= min(abs(y_cum-0.99));
            x_opt_t = X(1:in_op);
            y_opt_t = y(1:in_op);
            % Data used for fiting filter model
            x_opt  = 0:X(2)/10:X(in_op);
            y_opt  = spline(x_opt_t ,y_opt_t ,x_opt);

            % 2.2.3 Pars of the optimization algorithm
            %
%             y=PHIn(:,ind_t05)';
            y = ePSD_ti';
            % the Thompson normalized EPSD at t=0;
            s_w_spect=20;
            % Ration of the smoothing size of the moving average for a single
            % spectrum to the total bandwidth
            window_size =((eEPSD.L-1)/2)/s_w_spect;
            % windows size for smoothing the spectrum at each time.
            ysm= smooth(y',window_size);
            % Smooth version of y
            pks=findpeaks(ysm);
            % find a suitable starting point
            yt=sort(pks,'descend');
            if isempty(yt)
                yt = [ysm(51),ysm(51)];
            elseif numel(yt)==1
                yt = [yt,yt];
            end
            xl(1)=find(ysm==yt(1));
            xl(2)=find(ysm==yt(2));
            xl=sort(xl);

            type = SGMPar.Fmodulate.type;
            if strcmp(type,'IIorderFilter')||strcmp(type, 'KanaiTajimi')
                % starting point
                beta0 = [X(xl(1)),0.05,1];
                phi_0 = trapz(x_opt,y_opt)/trapz(x_opt,GMFilter(beta0,x_opt));
                beta0(3) = phi_0;
                % bounds parameters
                lb = [X(2),0.01 0];
                ub = [1/GM.dt*pi,0.99,inf];
                % Constraints
                A = []; b =[]; Aeq = [];beq=[];
            elseif strcmp(type, 'CloughPenzien')||strcmp(type, 'CascadeII-2')||strcmp(type,'Cascade_2-KT')
                % starting point
                beta0 = [X(xl(1)),0.05,X(xl(2)),0.05,1];
                phi_0 = trapz(x_opt,y_opt)/trapz(x_opt,GMFilter(beta0,x_opt));
                beta0(5) = phi_0;
                % bounds parameters
                lb = [X(2),0.01,X(2),0.01,0];
                ub = [x_opt_t(end),0.99,x_opt_t(end),0.99,inf];
                % Constraints
                A = [1,0,-1,0,0]; b=0; Aeq = []; beq=[];
            elseif strcmp(type, 'CloughPenzien_hf')
                % starting point
                beta0 = [X(xl(1)),0.05,X(xl(2)),1];
                phi_0 = trapz(x_opt,y_opt)/trapz(x_opt,GMFilter(beta0,x_opt));
                beta0(4) = phi_0;
                % bounds parameters
                lb = [X(2),0.01,X(2),0];
                ub = [1/GM.dt*pi,0.99,10*2*pi,inf];
                % Constraints
                A = [-1,0,1,0]; b=0; Aeq = []; beq=[];

            elseif strcmp(type, 'ConvexII-2')||strcmp(type,'Convex_2-KT')
                % starting point
                alpha_0 = y(xl(1))/(y(xl(1))+y(xl(2)));
                beta0 = [X(xl(1)),0.05,X(xl(2)),0.05,alpha_0,1];
                phi_0 = trapz(x_opt,y_opt)/trapz(x_opt,GMFilter(beta0,x_opt));
                beta0(6) = phi_0;
                % bounds parameters
                lb = [X(2),0.01,X(2),0.01,0,0];
                ub = [x_opt_t(end),0.99,x_opt_t(end),0.99,1,inf];
                % Constraints
                A = [1,0,-1,0,0,0];b = 0;Aeq = [];beq=[];
            elseif strcmp(type,'Convex-CP_II')
                % starting point
                alpha_0 = y(xl(1))/(y(xl(1))+y(xl(2)));
                beta0 = [X(xl(1)),0.05,0.01*2*pi,X(xl(2)),0.05,alpha_0,1]; % third par is like the fc
                phi_0 = trapz(x_opt,y_opt)/trapz(x_opt,GMFilter(beta0,x_opt));
                beta0(7) = phi_0;
                % bounds parameters
                lb = [X(2),0.01,01*2*pi,X(2),0.01,0,0];
                ub = [1/GM.dt*pi,0.99,10*2*pi,1/GM.dt*pi,0.99,1,inf];
                % Constraints
                A = [-1,0,1,0,0,0,0];b = 0;Aeq = [];beq=[];
            end

            % % % % %         if FilterNofPar == 3
            % % % % %             % starting point
            % % % % %             beta0=[X(xl(1)),0.2,1];
            % % % % %             % bounds parameters
            % % % % %             lb = [0 0 0];
            % % % % %             ub = [1/GM.dt*pi 1 inf];
            % % % % %         elseif FilterNofPar == 5
            % % % % %             % starting point
            % % % % %             beta0=[X(xl(1)),0.2,X(xl(1)),0.2,1];
            % % % % %             % bounds parameters
            % % % % %             lb = [0,0,0,0,0];
            % % % % %             ub = [1/GM.dt*pi,1, 1/GM.dt*pi,1,inf];
            % % % % %         elseif FilterNofPar == 6
            % % % % %             % starting point
            % % % % %             alpha_0 = y(xl(1))/(y(xl(1))+y(xl(2)));
            % % % % %             beta0=[X(xl(1)),0.05,X(xl(2)),0.05,alpha_0,1];
            % % % % %             phi_0 = trapz(x_opt,y_opt)/trapz(x_opt,GMFilter(beta0,x_opt));
            % % % % %             beta0(6)=phi_0;
            % % % % %             % bounds parameters
            % % % % %             lb = [X(2),0.05,X(2),0.05,0,0];
            % % % % %             ub = [1/GM.dt*pi,0.95, 1/GM.dt*pi,0.95,1,inf];
            % % % % %             A = [1,0,-1,0,0,0]; b =0; Aeq = [];beq=[];
            % % % % %         end

            % 2.2.4. strat the optimization
            % -----------------------------------------------------------------
            % ---- opt. fit method#1
            %         bn = lsqcurvefit(GMFilter,beta0,x_opt,y_opt,lb,ub,options);
            % ---- opt fit. method#2
            %         opts = statset('Display','off');
            %         mdl = fitnlm(x_opt,y_opt,GMFilter,beta0,'Options',opts);
            %         bn = mdl.Coefficients.Estimate;
            % ---- opt fit. method#3
            options = optimset('Display','off');
            objfun = @(b,x) sum((GMFilter(b,x_opt)-y_opt).^2);
%             objfun = @(b,x) sum(y_opt.*(GMFilter(b,x_opt)-y_opt).^2);
%             objfun = @(b,x) sum(y_opt.^2.*(GMFilter(b,x_opt)-y_opt).^2);
            bn = fmincon(objfun,beta0,A,b,Aeq,beq,lb,ub,[],options);
            % ---- opt fit. method#4
            %         options = optimoptions('particleswarm', 'Display', 'iter', 'SwarmSize', 100);
            %         objfun = @(x) sum((GMFilter(x,x_opt)-y_opt).^2);
            %         [bn, fval] = particleswarm(objfun,numel(beta0),lb,ub,options);
            % -----------------------------------------------------------------
%                             figure;
%                             plot(x_opt/2/pi,y_opt);hold on
%                             plot(x_opt/2/pi,GMFilter(beta0,x_opt))
%                             set(gca,'Xscale','linear','Yscale','linear')
%                             legend('data','model')

            % Save data
            FitResult.x_opt{i} = x_opt; FitResult.y_opt{i} = y_opt;
            FitResult.phi{i} = GMFilter(bn,x_opt);
            FittedFilterPar(i,:) = bn;

        else % ePSD_ti is constantly zeros, no needed to fit
            FittedFilterPar(i,:) = FittedFilterPar(i-1,:);
        end
    end

%% 2.2 Choice II -- log-y Sacle 
elseif strcmp(Fitscalse,'ylog')
    for i = ind_tp1:ind_tp2
        % ID of PHIn_ext
        ii = i+mid_p;
        % 2.2.1 Compute Wighting empirical PSD with a hann-window averaging at ti
        ePSD_ti = sum(PHIn_ext(:,(ii-mid_p):(ii+mid_p)).*WT,2);

        if sum(ePSD_ti)>0
            % 2.2.2 Determine the range to perform the fitting, i.e., optimize
            % up to 99% of the cumulative specrtum density
            X = eEPSD.wk;
            y = ePSD_ti/(trapz(X,ePSD_ti));
            y_cum = cumtrapz(X,y);
            [~, in_op]= min(abs(y_cum-0.99));
%             in_op =  find(X>=10*2*pi);
            x_opt_t = X(1:in_op);
            y_opt_t = y(1:in_op);
            % Data used for fiting filter model
            x_opt  = 1*2*pi:X(2)/10:X(in_op);
            y_opt  = spline(x_opt_t ,y_opt_t ,x_opt);
            logy_opt = log10(abs(y_opt));

            % 2.2.3 Pars of the optimization algorithm
            %
            y=log10(PHIn(:,ind_tp1))';
            % the Thompson normalized EPSD at t=0;
            s_w_spect=20;
            % Ration of the smoothing size of the moving average for a single
            % spectrum to the total bandwidth
            window_size =((eEPSD.L-1)/2)/s_w_spect;
            % windows size for smoothing the spectrum at each time.
            ysm= smooth(y',window_size);
            % Smooth version of y
            pks=findpeaks(ysm);
            % find a suitable starting point
            yt=sort(pks,'descend');
            if isempty(yt)
                xl = 51;
            end
            xl(1)=find(ysm==yt(1));
            xl(2)=find(ysm==yt(2));
            xl=sort(xl);

            GMFilter_logy = @(b,x) log10(GMFilter(b,x));
            type = SGMPar.Fmodulate.type;
            if strcmp(type,'IIorderFilter')||strcmp(type, 'KanaiTajimi')
                % starting point
                beta0 = [X(xl(1)),0.05,1];
                phi_0 = trapz(x_opt,y_opt)/trapz(x_opt,GMFilter(beta0,x_opt));
                beta0(3) = phi_0;
                % bounds parameters
                lb = [X(2),0.01 0];
                ub = [1/GM.dt*pi,0.95,inf];
                % Constraints
                A = []; b =[]; Aeq = [];beq=[];
            elseif strcmp(type, 'CloughPenzien')||strcmp(type, 'CascadeII-2')
                % starting point
                beta0 = [X(xl(1)),0.5,X(xl(2)),0.05,1];
                phi_0 = trapz(x_opt,y_opt)/trapz(x_opt,GMFilter(beta0,x_opt));
                beta0(5) = phi_0;
                % bounds parameters
                lb = [X(2),0.01,X(2),0.01,0];
                ub = [1/GM.dt*pi,0.99,1/GM.dt*pi,0.99,inf];
                % Constraints
                A = [1,0,-1,0,0]; b=0; Aeq = []; beq=[];

            elseif strcmp(type, 'ConvexII-2')
                % starting point
                alpha_0 = y(xl(1))/(y(xl(1))+y(xl(2)));
                beta0 = [X(xl(1)),0.05,X(xl(2)),0.05,alpha_0,1];
                phi_0 = trapz(x_opt,y_opt)/trapz(x_opt,GMFilter(beta0,x_opt));
                beta0(6) = phi_0;
                % bounds parameters
                lb = [X(2),0.05,X(2),0.05,0,0];
                ub = [x_opt_t(end),0.99, x_opt_t(end),0.99,1,inf];
                % Constraints
                A = [1,0,-1,0,0,0];b = 0;Aeq = [];beq=[];
            elseif strcmp(type,'Convex-CP_II_v1') 
                % starting point
                alpha_0 = y(xl(1))/(y(xl(1))+y(xl(2)));
                beta0 = [X(xl(1)),0.05,0.01*2*pi,X(xl(2)),0.05,alpha_0,1]; % third par is like the fc
                phi_0 = trapz(x_opt,y_opt)/trapz(x_opt,GMFilter(beta0,x_opt));
                beta0(7) = phi_0;
                % bounds parameters
                lb = [X(2),0.05,01*2*pi,X(2),0.05,0,0];
                ub = [1/GM.dt*pi,0.95,10*2*pi,1/GM.dt*pi,0.95,1,inf];
                % Constraints
                A = [-1,0,1,0,0,0,0];b = 0;Aeq = [];beq=[];
            elseif strcmp(type,'Convex-CP_II')
                % starting point
                alpha_0 = abs(y(xl(1))/(y(xl(1))+y(xl(2))));
                beta0 = [0.01*2*pi,0.05,X(xl(1)),0.05,X(xl(2)),0.05,alpha_0,1]; 
                phi_0 = trapz(x_opt,y_opt)/trapz(x_opt,GMFilter(beta0,x_opt));
                beta0(7) = phi_0;
                % bounds parameters
                lb = [X(2),0.01,X(2),0.01,X(2),0.01,0,0];
                ub = [1/GM.dt*pi,0.99,1/GM.dt*pi,0.99,1/GM.dt*pi,0.99,1,inf];
                % Constraints
                A = [1,0,-1,0,0,0,0,0
                     0,0,1,0,-1,0,0,0];
                b = [0;0];
                Aeq = [];beq=[];
            end       


      
            % 2.2.4. strat the optimization
            % -----------------------------------------------------------------
            % ---- opt. fit method#1
            %         bn = lsqcurvefit(GMFilter,beta0,x_opt,y_opt,lb,ub,options);
            % ---- opt fit. method#2
            %         opts = statset('Display','off');
            %         mdl = fitnlm(x_opt,y_opt,GMFilter,beta0,'Options',opts);
            %         bn = mdl.Coefficients.Estimate;
            % ---- opt fit. method#3
            options = optimset('Display','off');
            objfun = @(b,x) sum(y_opt.*(GMFilter_logy(b,x_opt)-logy_opt).^2 );
            bn = fmincon(objfun,beta0,A,b,Aeq,beq,lb,ub,[],options);
%             figure
%             plot(X/2/pi,ysm)
%             hold on
%             plot(x_opt/2/pi,logy_opt)
% %             hold on
% % %             plot(x_opt/2/pi,GMFilter_logy(beta0,x_opt)) 
% %             hold on
%             plot(x_opt/2/pi,GMFilter_logy(bn,x_opt),'--','LineWidth',3)
%             legend('data','truncated data','fitted')

            % ---- opt fit. method#4
            %         options = optimoptions('particleswarm', 'Display', 'iter', 'SwarmSize', 100);
            %         objfun = @(x) sum((GMFilter(x,x_opt)-y_opt).^2);
            %         [bn, fval] = particleswarm(objfun,numel(beta0),lb,ub,options);
            % -----------------------------------------------------------------
            %                 figure;
            %                 plot(x_opt,y_opt);hold on
            %                 plot(x_opt,GMFilter(bn,x_opt))
            %                 legend('data','model')

            % Save data
            FitResult.x_opt{i} = x_opt; FitResult.y_opt{i} = y_opt;
            FitResult.phi{i} = GMFilter(bn,x_opt);
            FittedFilterPar(i,:) = bn;

        else % ePSD_ti is constantly zeros, no needed to fit
            FittedFilterPar(i,:) = FittedFilterPar(i-1,:);
        end
    end
end


%%
% 2.2.5. Complement the BN
FittedFilterPar(1:ind_tp1,:) = ones(ind_tp1,1)*FittedFilterPar(ind_tp1,:);
FittedFilterPar(ind_tp2:GM.L,:) = ones(numel(ind_tp2:GM.L),1)*FittedFilterPar(ind_tp2,:);

% figure(Position=[100,100,800,300])
% plot(GM.teq,FittedFilterPar(:,1))
% hold on
% plot(GM.teq,FittedFilterPar(:,3))
% xlabel('$t$ [s]')
% ylabel('$\omega_g(t)$')
% grid on

%% ------------------------------------------------------------------------ 
%             2-II. Fit trend functions for estimated value BN
%% ------------------------------------------------------------------------
Paramtericmodel = SGMPar.Fmodulate.Paramtericmodel;
Filter_opt_par = Fmodulate_fitEPSDPar(Paramtericmodel,FittedFilterPar,GM,SGMPar.Tmodulate);


%% ------------------------------------------------------------------------ 
%                             Save  Data
%% ------------------------------------------------------------------------
SGMPar.Fmodulate.BN = FittedFilterPar;
SGMPar.Fmodulate.opt_par = Filter_opt_par;


end
