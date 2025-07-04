function   [EPSD_unit_var,fitted_opt_par_smooth] =  FModulateComputeEPSD(Tmodulate,Fmodulate,L,dt,wk)
%% ------------------------------------------------------------------------
%                  Parametric model for EPSD Pars
%% ------------------------------------------------------------------------
%  polymodel input: Par.wg_p and Par.zeta_p
%  Linear input: Par.wg_tmid, Par.wg_rate, Par.zeta_tmid, Par.zeta_rate    
%  Constant input: Par.wg_tmid, Par.zeta_tmid

Paramtericmodel = Fmodulate.Paramtericmodel;
if isfield(Fmodulate,'opt_par')
    opt_par = Fmodulate.opt_par;
else
    Par = Fmodulate.Par;
    opt_par = [Par.wg_tmid,Par.zeta_tmid;
               Par.wg_rate,Par.zeta_rate];
end
% Par = Fmodulate.Par;
teq = (0:L-1)*dt;

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
[GMFilter,FilterNofPar] = SetFilterModel(Fmodulate.type);

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

%% Predict the Trend of filter Parameters
for nn = 1:FilterNofPar-1
    switch TrendType(nn)

        case 1 % Costant (One Par)
          fitted_opt_par(:,nn) = ones(L,1)*opt_par(1,nn);

        case 2 % LinearModel (Two Pars)
            if strcmp(Tmodulate.type,'spline')
                id = find(Tmodulate.Par.p==45);
                tmid = Tmodulate.Par.t_AI_p(id);
            elseif strcmp(Tmodulate.type,'Gamma')
                tmid = Tmodulate.Par.tmid;
            end
            Linearmodelfun = @(b,t) (b(1)+b(2)*(t-tmid));     
            PI_val1 = opt_par(1:2,nn);
            fitted_opt_par(:,nn) = Linearmodelfun(PI_val1,teq)';
            

        case 3 % Polyline (Five Pars)
            if strcmp(Tmodulate.type,'spline')
                t_AI_p = Tmodulate.Par.t_AI_p;
            else
                error('Polymodel is not recomended when Spline model is unused in time-modulating mdoel')
            end
            PI_val1 = opt_par(:,nn);
            PI_val1 = [PI_val1(1),PI_val1',PI_val1(end)];
            fitted_opt_par(:,nn) = interp1(t_AI_p,PI_val1,teq);
    end
end


%% Old-code
% % if strcmp(Paramtericmodel,'Polymodel')
% %     if strcmp(Tmodulate.type,'spline')
% %         t_AI_p = Tmodulate.Par.t_AI_p;
% %     else
% %         error('Polymodel is not recomended when Spline model is unused in time-modulating mdoel')
% %     end
% % 
% %     if FilterNofPar == 3 % i.e. II-order filter,KT filter,...
% %         Linear interpolation
% %         interpolate the frequency
% %         PI_val1 = Par.wg_p; PI_val1 = [PI_val1(1),PI_val1',PI_val1(end)];
% %         BW1 = interp1(t_AI_p,PI_val1,teq);
% %         interpolate the dampering ratio
% %         PI_val2 = Par.zeta_p; PI_val2 = [PI_val2(1),PI_val2',PI_val2(end)];
% %         BW2 = interp1(t_AI_p,PI_val2,teq);
% %         Bound the parameters
% %         fitted_opt_par = [BW1(:),BW2(:)];
% %     elseif FilterNofPar == 6||FilterNofPar == 5
% %         for nn = 1:FilterNofPar-1
% %             PI_val1 = opt_par(:,nn); PI_val1 = [PI_val1(1),PI_val1',PI_val1(end)];
% %             fitted_opt_par(:,nn) = interp1(t_AI_p,PI_val1,teq);
% %             PI_val1 = opt_par(:,nn);
% %             fitted_opt_par(:,nn) = Linearmodelfun(PI_val1,teq)';
% %         end
% %     end
% %  
% % elseif strcmp(Paramtericmodel,'Linearmodel')||strcmp(Paramtericmodel,'ReducedLinearmodel')
% %     
% %     if strcmp(Tmodulate.type,'spline')
% %         id = find(Tmodulate.Par.p==45);
% %         tmid = Tmodulate.Par.t_AI_p(id);
% %     elseif strcmp(Tmodulate.type,'Gamma')
% %         tmid = Tmodulate.Par.tmid;
% %     end
% %     Linearmodelfun = @(b,t) (b(1)+b(2)*(t-tmid)); 
% %     for nn = 1:FilterNofPar-1
% %         PI_val1 = opt_par(:,nn);
% %         fitted_opt_par(:,nn) = Linearmodelfun(PI_val1,teq)';
% %     end
% % 
% % elseif strcmp(Paramtericmodel,'Constantmodel')
% % 
% %     if FilterNofPar == 3
% %         wg_ti = ones(1,L)*Par.wg_tmid;
% %         zeta_ti = ones(1,L)*Par.zeta_tmid;
% %         fitted_opt_par = [wg_ti(:),zeta_ti(:)];
% %     elseif FilterNofPar == 6||FilterNofPar == 5
% %         fitted_opt_par = ones(L,1)*opt_par(1,1:end-1);
% %     end
% % 
% % end



%% ------------------------------------------------------------------------
%                  Bound and smooth fitted_opt_par
%% ------------------------------------------------------------------------
% % Truncate t<t5 and t>t95
[t_q_pad,qmod,padn] = TModulate_Estimator(Tmodulate.Par,dt,L);
AI_q = cumtrapz(t_q_pad,qmod.^2);
AI_n = AI_q/AI_q(end)*100;
[~,id_t5] = min(abs(AI_n-5));
[~,id_t95] = min(abs(AI_n-95));
fitted_opt_par(1:id_t5,:) =  repmat(fitted_opt_par(id_t5,:),id_t5,1);
fitted_opt_par(id_t95:end,:) = repmat(fitted_opt_par(id_t95,:),L+1-id_t95,1);

% % Bound
if FilterNofPar == 3
    % Bound wg_ti
    fitted_opt_par(fitted_opt_par(:,1)<=25/500*2*pi,1) = 25/500*2*pi;  
    % Bound zeta_ti
    fitted_opt_par(fitted_opt_par(:,2)<=0.02,2) = 0.02;
    fitted_opt_par(fitted_opt_par(:,2)>=0.98,2) = 0.98;
elseif FilterNofPar == 6||FilterNofPar == 5
    % Bound wg_ti
    fitted_opt_par(fitted_opt_par(:,1)<=25/500*2*pi,1) = 25/500*2*pi;
    fitted_opt_par(fitted_opt_par(:,3)<=25/500*2*pi,3) = 25/500*2*pi;
    % Bound zeta_ti
    fitted_opt_par(fitted_opt_par(:,2)<=0.02,2) = 0.02;
    fitted_opt_par(fitted_opt_par(:,2)>=0.98,2) = 0.98;
    fitted_opt_par(fitted_opt_par(:,4)<=0.02,4) = 0.01;
    fitted_opt_par(fitted_opt_par(:,4)>=0.98,4) = 1;
    if  FilterNofPar == 6
        % Bound the paratameter of alpha
        fitted_opt_par(fitted_opt_par(:,5)>=1,5) = 1-1e-6;
        fitted_opt_par(fitted_opt_par(:,5)<=0,5) = 1e-6;
    end
end

% Smooth the parameters
s_time = 2;
for nn = 1:FilterNofPar-1
    fitted_opt_par_smooth(:,nn) = smooth(fitted_opt_par(:,nn),s_time/dt);
end

%% ------------------------------------------------------------------------
%                            Estimated EPSD
%% ------------------------------------------------------------------------
% modelfun = @(b,x)((b(1).^4)./...
%     ((b(1).^2.-x.^2).^2.+4.*b(2).^2.*b(1).^2.*x.^2)).*b(3);

% wk = SGMPar.Fmodulate.wk;
% L = L;
EPSD_unit_var(1:numel(wk),1:L) = 0;

for ii = 1:L
%     fitted_opt_par = 1
%     fitted_opt_par = [wg_ti(ii),zeta_ti(ii),1];
    Phitemp = GMFilter([fitted_opt_par_smooth(ii,:),1],wk)';
    EPSD_unit_var(:,ii) = Phitemp/(2*trapz(wk,Phitemp));
end

dw = wk(2)-wk(1);
NormFactor = 2*sum(EPSD_unit_var)*dw;
EPSD_unit_var = EPSD_unit_var./(repmat(NormFactor,numel(wk),1));
% NormFactor = 2*sum(EPSD_unit_var)*dw;

%   temp = sum(EPSD_unit_var/dw,1);
% temp = sum(sqrt(EPSD_unit_var_pad),1);
% temp = sum(EPSD_unit_var*dw,1)

end