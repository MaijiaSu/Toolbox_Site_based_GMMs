% Define the Rectangle class
classdef PlotGM
    properties
        GMi          % 
        SGMPari
        FFT
        eEPSD
        SA
        NLSA
        BWResult
    end
    
    methods

        % Constructor
        function obj = PlotGM(GMi)
            % 1. GM time history 
            % discrete time instant
            GMi.teq = (0:GMi.L-1)*GMi.dt;
            % Compute the time instant 'GM.t_AI_p' that reaches p% of the total Arial Intensity
            AI = cumtrapz(GMi.teq,GMi.eq.^2);
            AI_n = AI/AI(end)*100;
            [~,ID_repeat] = unique(AI_n);
            GMi.AI_p = [0.001,0.05,1:99,99.5,99.99];
            GMi.t_AI_p = interp1(AI_n(ID_repeat),GMi.teq(ID_repeat),GMi.AI_p);
            % Complement the starting time and the ending time
            GMi.AI_p = [0,GMi.AI_p,100];
            GMi.AI = GMi.AI_p/100*AI(end);
            GMi.t_AI_p = [0,GMi.t_AI_p,GMi.teq(end)];
            % Compte the ID in GM.eq that reaches p% of the total Arial Intensity
            for nn = 1:numel(GMi.AI_p)
                GMi.ID_p(nn) = find(AI_n>=GMi.AI_p(nn),1,'first');
            end
            obj.GMi = GMi;   
        end
        
        % 1. Plot Real GM time series
        function obj = PlotTime(obj)
            dt = obj.GMi.dt; L = obj.GMi.L;
            eq = obj.GMi.eq; teq = (0:L-1)*dt;  
%             figure(Position=[100,100,800,300])
            plot(teq,eq,'k','LineWidth',1)
            grid on
            xlim([0,max(teq)])
            xlabel('$t$ [s]')
            ylabel('$a_g(t)$')
        end

        % 1-2：Plot Both real and simulated time series
        function obj = PlotSimTimeSeries(obj,simGMi)
            GMi = obj.GMi;
            dt = GMi.dt;
            L = numel(GMi.eq);
            teq = (0:L-1)*dt;
            ag_real = GMi.eq;
            Vg_real = cumtrapz(teq,ag_real);
            Dg_real = cumtrapz(teq,Vg_real);

            % Simulated
            dt = simGMi.dt;
            L = numel(simGMi(1).eq);
            teq_sim = (0:L-1)*dt;
            ag_sim = simGMi(1).eq;
            Vg_sim = cumtrapz(teq_sim,ag_sim);
            Dg_sim = cumtrapz(teq_sim,Vg_sim);


            figure(OuterPosition=[100,100,1000,300])
            %         subplot(3,1,1)
            plot(teq_sim,ag_sim,'LineWidth',1)
            hold on
            plot(teq,ag_real,'LineWidth',1)
            ylabel('$A(t)$');xlabel('$t$')
            grid on
            legend('sim.','real','Location','best')
             % Compute time-modulating function
%             [t_q_pad,qmod,padn] = TModulate_Estimator(SGMPar_ir.Tmodulate.Par,GMi.dt,GMi.L);
%             hold on
%             plot(t_q_pad,qmod,'--g','LineWidth',2)


            %         subplot(3,1,2)
            %         plot(t0+teq_sim,Vg_sim,'LineWidth',1)
            %         hold on
            %         plot(teq,Vg_real,'LineWidth',1)
            %         ylabel('$V(t)$');xlabel('$t$')
            %         grid on
            %         wg_mid  = SGMPar.Fmodulate.Par.wg_tmid;
            %         wg_rate = SGMPar.Fmodulate.Par.wg_rate;
            %         zeta_mid  = SGMPar.Fmodulate.Par.zeta_tmid;
            %         zeta_rate = SGMPar.Fmodulate.Par.zeta_rate;
            %
            %         subplot(3,1,3)
            %         plot(t0+teq_sim,Dg_sim,'LineWidth',1)
            %         hold on
            %         plot(teq,Dg_real,'LineWidth',1)
            %         ylabel('$D(t)$');xlabel('$t$')
        end

        % 2. Plot Sa(T)
        function obj = PlotSa(obj,varargin)

            if size(varargin,2) == 0
                T_t = 0.01:0.01:1.5;T_t = [T_t 1.6:0.1:10.1];
                zeta = 0.05;
            elseif size(varargin,2) == 1
                T_t = varargin{1}; zeta = varargin{2};
            else
                error('wromg PlotSa input')
            end        
            GMi = obj.GMi;g = 9.816;
            [Sd,Sv,Sa]  = ResponseSpectraNM(GMi.eq,GMi.dt,g,T_t,zeta,0,0,1/2,1/6);
            obj.SA.T_t = T_t; obj.SA.zeta = zeta; 
            obj.SA.Sa = Sd; obj.SA.Sv = Sv; obj.SA.Sa = Sa; 

%             figure(Position=[100,100,500,400])
            plot(T_t,Sa,'-','Linewidth',1.5,'Color','r');
            set(gca,'XScale','log')
            set(gca,'YScale','Linear')
            xlabel('$T[s]$')
            ylabel('$S_a(T)$')
            xlim([0.05,10])
            xticks([0.05,0.1,1,10])
%             ylim([0,1.2])
            grid on                   
        end

         % 2-1. Plot Sa(T)
        function obj = PlotlogSa(obj,varargin)
            if nargin == 1
                T_t = 0.01:0.01:1.5;T_t = [T_t 1.6:0.1:10.1];
                zeta = 0.05;
            elseif nargin == 2
                T_t = varargin{1}; zeta = varargin{2};
            else
                error('wromg PlotSa input')
            end        
            GMi = obj.GMi;g = 9.816;
            [Sd,Sv,Sa]  = ResponseSpectraNM(GMi.eq,GMi.dt,g,T_t,zeta,0,0,1/2,1/6);
            obj.SA.T_t = T_t; obj.SA.zeta = zeta; 
            obj.SA.Sa = Sd; obj.SA.Sv = Sv; obj.SA.Sa = Sa; 

%             figure(Position=[100,100,500,400])
            plot(T_t,log(Sa),'-','Linewidth',1.5,'Color','k');
%             set(gca,'XScale','log')
            set(gca,'YScale','Linear')
            xlabel('$T[s]$')
            ylabel('$S_a(T)$')
            xlim([0.05,10])
            xticks([0.05,0.1,1,10])
            % ylim([0,1.2])
            grid on                   
        end

        % 3. Plot eEPSD
        function obj =  ploteEPSD(obj,SGMPar)
            % eEPSD
            GMi = obj.GMi;
            TMWSEPar.TIME_Windows_f = 4;TMWSEPar.TIME_Windows_t = 4;                                   
            TMWSEPar.L = ceil((TMWSEPar.TIME_Windows_f/GMi.dt)/2)*2+1;                                  
            TMWSEPar.Lt = ceil((TMWSEPar.TIME_Windows_t/GMi.dt)/2)*2+1;      
            TMWSEPar.K = 5;  TMWSEPar.dt = GMi.dt;
            [PHI, PHIn, sgt_2, S,w] = TMWSE(GMi.eq,TMWSEPar);      
            eEPSD.PHI = PHI;eEPSD.sgt_2 = sgt_2;eEPSD.PHIn = PHIn;
            eEPSD.wk = (0:1:(TMWSEPar.L-1)/2)*2*pi/(GMi.dt*TMWSEPar.L);
            eEPSD.S = S;eEPSD.w = w;eEPSD.L = numel(w);
            eEPSD.TMWSEPar = TMWSEPar;
            obj.eEPSD = eEPSD;

            teq = (0:obj.GMi.L-1)*obj.GMi.dt; 
            fk = obj.eEPSD.wk/2/pi;
            [X1,X2] = meshgrid(teq,fk);
%             figure(Position=[100,100,500,400])
            contour(X1,X2,(abs(obj.eEPSD.PHI)),50)
            colormap('jet')
            shading interp; alpha 0.5;  box on
            xlabel('$t$ [s]');ylabel('$f$ [Hz]')
            xlim([0,max(teq)]);
            PP = (abs(obj.eEPSD.PHI));PP = sum(PP,2); cumPP = cumsum(PP);
            fkk =  fk(cumPP>0.99*cumPP(end));
            ylim([0,min(fkk)])
            eEPSD.fkk = fkk;
            obj.eEPSD = eEPSD; 
            grid on      

%             hold on
%             wk = SGMPar.Fmodulate.wk;
%             [~,fitted_opt_par] =  FModulateComputeEPSD(SGMPar.Tmodulate,...
%                 SGMPar.Fmodulate,GMi.L,GMi.dt,wk);
%             hold on
%             plot(teq,fitted_opt_par(:,1)/2/pi,'LineWidth',5,'Color',[6,6,6]/10)
 

        end

         % 3. Plot fitted eEPSD
        function obj =  plotFittedEPSD(obj,SGMPar)        
            GMi = obj.GMi;
            teq = (0:GMi.L-1)*GMi.dt; 
            wk = SGMPar.Fmodulate.wk;
            fk = wk/2/pi;
            [X1,X2] = meshgrid(teq,fk);
           
            [EPSD_unit_var,fitted_opt_par] =  FModulateComputeEPSD(SGMPar.Tmodulate,...
                SGMPar.Fmodulate,GMi.L,GMi.dt,wk);
            [~,qmod,~] = TModulate_Estimator(SGMPar.Tmodulate.Par,GMi.dt,GMi.L);
             EPSD_q = EPSD_unit_var.*repmat(qmod.^2,numel(wk),1);
            % M1
%             surfc(X1,X2,EPSD_unit_var);
%             polarmap([cmap('red',1,0,0);flipud(cmap('red',20,30,5))])
%             polarmap(300,0.9)
%             shading interp
%             alpha 0.7
%             view(0,90)

            % M2
            contour(X1,X2,EPSD_q,50);
            colormap('jet')
            shading interp; alpha 0.5;  box on

            xlabel('$t(s)$');ylabel('$f$(Hz)')
            xlim([0,max(teq)]);
            PP = sum(EPSD_q,2); cumPP = cumsum(PP);
            fkk =  fk(cumPP>0.99*cumPP(end));
            ylim([0,min(obj.eEPSD.fkk)])
%             grid on      
            hold on
            if size(fitted_opt_par,2) == 2
                p1 = plot(teq,fitted_opt_par(:,1)/2/pi,'LineWidth',3,'Color',[7,7,7]/10);
                legend(p1,'$\omega_{g}$')    
            elseif size(fitted_opt_par,2) >=3
               p1 = plot(teq,fitted_opt_par(:,1)/2/pi,'LineWidth',3,'Color',[4,4,4]/10);
                hold on
               p2 = plot(teq,fitted_opt_par(:,3)/2/pi,'--','LineWidth',3,'Color',[2,2,2]/10);
                legend([p1,p2],'$\omega_{g,1}$','$\omega_{g,2}$')      
            end

        end

        % 4. Plot FFT
        function obj =  plotFFT(obj)
%             figure(Position=[100,100,500,400])
            % 2. FFT of real GM
            GMi = obj.GMi;
            fs = 1/GMi.dt;
            Y = fft(GMi.eq);f_real = fs*(0:(GMi.L/2))/GMi.L;
            P2 = abs(Y/GMi.L);P1_real = P2(1:GMi.L/2+1);
            P1_real(2:end-1) = 2*P1_real(2:end-1);
            obj.FFT.P1 = P1_real;
            obj.FFT.f = f_real;
           
%             E_f = cumtrapz(f_real.f,P1_real.^2);
%             E_f_n = E_f/E_f(end)*100;
%             ff = [0.01:0.01:0.3,0.4:0.1:10];
%             E_f_n_feq1(ir,:) = interp1(GM(ir).f,E_f_n,ff);

            plot(obj.FFT.f,obj.FFT.P1,'-','Linewidth',1.5,'Color','k');
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            xlabel('$f$ [Hz]')
            ylabel('$A(f)$')
            xlim([0.05,max( obj.FFT.f)])
            grid on
        end


        % Plot FFT V1
        function obj =  ComputeFFT(obj)
            %             figure(Position=[100,100,500,400])
            % 2. FFT of real GM
            GMi = obj.GMi;
            fs = 1/GMi.dt;
            Y = fft(GMi.eq);f_real = fs*(0:(GMi.L/2))/GMi.L;
            P2 = abs(Y/GMi.L);P1_real = P2(1:GMi.L/2+1);
            P1_real(2:end-1) = 2*P1_real(2:end-1);
            obj.FFT.P1 = P1_real;
            obj.FFT.f = f_real;

            %             E_f = cumtrapz(f_real.f,P1_real.^2);
            %             E_f_n = E_f/E_f(end)*100;
            %             ff = [0.01:0.01:0.3,0.4:0.1:10];
            %             E_f_n_feq1(ir,:) = interp1(GM(ir).f,E_f_n,ff);

%             plot(obj.FFT.f,obj.FFT.P1,'-','Linewidth',1.5,'Color','k');
%             set(gca,'XScale','log')
%             set(gca,'YScale','log')
%             xlabel('$f$ [Hz]')
%             ylabel('$A(f)$')
%             xlim([0.05,max( obj.FFT.f)])
%             grid on
        end

      % 4. Plot FFT_v2
        function obj =  plotFFTv2(obj)
%             figure(Position=[100,100,500,400])
            % 2. FFT of real GM
            GMi = obj.GMi;
            fs = 1/GMi.dt;
            Y = fft(GMi.eq);f_real = fs*(0:(GMi.L/2))/GMi.L;
            P2 = abs(Y/GMi.L);P1_real = P2(1:GMi.L/2+1);
            P1_real(2:end-1) = 2*P1_real(2:end-1);
            obj.FFT.P1 = P1_real;
            obj.FFT.f = f_real;
           
            subplot(3,1,1)
            plot(obj.FFT.f,obj.FFT.P1,'-','Linewidth',1.5,'Color','k');
            set(gca,'XScale','log'); set(gca,'YScale','log');xlabel('$f$ [Hz]');ylabel('$A(f)$')
            grid on

            subplot(3,1,2)
            plot(obj.FFT.f,obj.FFT.P1./(2*pi*obj.FFT.f)','-','Linewidth',1.5,'Color','k');
            set(gca,'XScale','log'); set(gca,'YScale','log');xlabel('$f$ [Hz]');ylabel('$V(f)$')
            grid on

            subplot(3,1,3)
            plot(obj.FFT.f,obj.FFT.P1./(2*pi*obj.FFT.f)'.^2,'-','Linewidth',1.5,'Color','k');
            set(gca,'XScale','log'); set(gca,'YScale','log');xlabel('$f$ [Hz]');ylabel('$D(f)$')
            grid on     
%             xlim([0.05,max( obj.FFT.f)])            
        end

        % 4. Plot FFT
        function obj =  plotFFT_Sim(obj,SGMPar)
%             figure(Position=[100,100,500,400])
            % 2. FFT of real GM
            GMi = obj.GMi;
            fs = 1/GMi.dt;
            Y = fft(GMi.eq);f_real = fs*(0:(GMi.L/2))/GMi.L;
            P2 = abs(Y/GMi.L);P1_real = P2(1:GMi.L/2+1);
            P1_real(2:end-1) = 2*P1_real(2:end-1);
            obj.FFT.P1 = P1_real;
            obj.FFT.f = f_real;
           


            % FAS of Simulated GM
%             Option.T_t =T_t;
            Option.NofSGM = 100;
            Option.ComputeSa = 'off';
            dt = 0.02;
            [simGM,~] = ...
                SGMgeneration(SGMPar,dt,[],Option);
            for nn = 1:Option.NofSGM
                ag_sim_T(nn,:) = (simGM(nn).eq);
            end
            L = numel(simGM(nn).eq);
            fs = 1/dt;
            t = (0:L-1)*dt;              % Time vector
            f = fs*(0:(L/2))/L;
            Y = fft(ag_sim_T,[],2);
            P2 = abs(Y/L);
            P1_sim = P2(:,1:L/2+1);
            P1_sim(:,2:end-1) = 2*P1_sim(:,2:end-1);
            hold on
            p1 = plot(f,P1_sim,'-','Color',[9,9,9]/10);
            mean_P1_sim = mean(P1_sim);
            p2 = plot(f,mean_P1_sim,'-','Linewidth',1.5,'Color','b');
            % Real
            hold on
            p3 = plot(obj.FFT.f,obj.FFT.P1,'-','Linewidth',1.5,'Color','k');
            hold on

            % Plot setting
            set(gca,'XScale','log')
            set(gca,'YScale','log')
            xlabel('$f$ [Hz]')
            ylabel('$A(f)$')
            xlim([0.05,max( obj.FFT.f)])
            legend([p1(1),p2,p3],'single sim.','sim. mean','real', ...
                'Location','best')
            grid on
        end
       


        % 5. AnimateEPSD
        function obj = animateEPSD(obj,Type,filename)
            if ~isfield(obj,'eEPSD')
                GMi = obj.GMi;
                TMWSEPar.TIME_Windows_f = 4;TMWSEPar.TIME_Windows_t = 4;
                TMWSEPar.L = ceil((TMWSEPar.TIME_Windows_f/GMi.dt)/2)*2+1;
                TMWSEPar.Lt = ceil((TMWSEPar.TIME_Windows_t/GMi.dt)/2)*2+1;
                TMWSEPar.K = 5;  TMWSEPar.dt = GMi.dt;
                [PHI, PHIn, sgt_2, S,w] = TMWSE(GMi.eq,TMWSEPar);
                eEPSD.PHI = PHI;eEPSD.sgt_2 = sgt_2;eEPSD.PHIn = PHIn;
                eEPSD.wk = (0:1:(TMWSEPar.L-1)/2)*2*pi/(GMi.dt*TMWSEPar.L);
                eEPSD.S = S;eEPSD.w = w;eEPSD.L = numel(w);
                eEPSD.TMWSEPar = TMWSEPar;
                obj.eEPSD = eEPSD;
            end

            % % prepare data for real EPSD
            fs_win = (0:obj.eEPSD.TMWSEPar.L/2)/(obj.eEPSD.TMWSEPar.TIME_Windows_f);
            sh = (obj.eEPSD.TMWSEPar.L-1)/2; eq_pad=padarray(obj.GMi.eq,[sh 0]);
            t_pad = (0:numel(eq_pad)-1)*obj.GMi.dt;
            teq = (0:obj.GMi.L-1)*obj.GMi.dt;
            w = obj.eEPSD.w;  tw = (0:numel(w)-1)*obj.GMi.dt;
            % % prepare data for simulated EPSD
            K_smooth = 0;K = obj.GMi.L+K_smooth;K = max(K,1000);CUT_off_FREQ = 25;
            wk = linspace(0,CUT_off_FREQ*2*pi,floor(K/2)+rem(K,2));
            SGMPar.Tmodulate.type = 'spline'; SGMPar.Tmodulate.s_time = 0;
            SGMPar = FitSGMPar_Tmodulate(obj.GMi,SGMPar);
            SGMPar.Fmodulate.wk = wk; SGMPar.Fmodulate.type = Type;
            SGMPar.Fmodulate.Paramtericmodel = 'Linearmodel';
            [SGMPar,FitResult] = FitSGMPar_Fmodulate(obj.GMi,SGMPar,obj.eEPSD);
            obj.eEPSD.FitResult = FitResult;

            % % Initial plot
            myfig = figure('OuterPosition',[100,100,800,400]);
            tilefig = tiledlayout(2,4,"TileSpacing","tight","Padding","tight");
            nexttile([1,4])
            plot(t_pad-2,eq_pad)
            ylim([-max(abs(eq_pad)),max(abs(eq_pad))]);
            xlabel('$t(s)$');ylabel('$\ddot{x}_g(t)$') 
            hold on; yyaxis right
            p1 = plot(tw-2,w,'LineWidth',1.5);
            ylabel('$w_{i}(t)$');   hold on
            scatter([0,max(teq)],[1,1]*max(w)*1.05,'Marker','v','MarkerEdgeColor','k','MarkerFaceColor','g')
            TIME_Windows_t = obj.eEPSD.TMWSEPar.TIME_Windows_t;
            xlim([min(t_pad)-TIME_Windows_t/2,max(t_pad)-TIME_Windows_t/2])
            grid on; 
            
            nexttile([1,2]);p2 = plot(tw,obj.eEPSD.S(1,:));
            xlabel('$t(s)$'); ylabel('$\ddot{x}_g(t) \times w_{i}(t)$')
            grid on; box on;
            ylim(max(max(abs(obj.eEPSD.S)))*[-1,1])
            
            nexttile([1,2])
            ind_t05 = obj.GMi.ID_p(2);         % index of time at 0.05% of Arias Intensity
            ind_tf_9950 = obj.GMi.ID_p(end-1); % index of time at 99.50% of Arias Intensity
            p3 = plot(FitResult.x_opt{ind_t05}/2/pi,FitResult.y_opt{ind_t05});
            hold on
            p4 = plot(FitResult.x_opt{ind_t05}/2/pi,FitResult.phi{ind_t05});
            grid on
            xlabel('$f(Hz)$');ylabel('PSD$_i(f)$');  
            legend('empircal','fitted')
            for nn = ind_t05:ind_tf_9950
                 m1(nn) = max(FitResult.y_opt{nn});m2(nn) = max(FitResult.phi{nn});      
            end    
            ylim([0, max([m1;m2],[],'all') ])
            xlim([0,25])
            % % Generate the fame
            % Create file name variable 
%             frameNo = [1:10:length(obj.GMi.eq),numel(obj.GMi.eq)];
            frameNo  =  [1,floor(linspace(ind_t05,ind_tf_9950,200)),obj.GMi.L]
            loops = numel(frameNo);
            ii = 0;im = [];   
            for kk = frameNo
                % update the data
                % 1
                p1.XData = tw-2+(kk-1)*obj.GMi.dt;
                % 2
                p2.YData = obj.eEPSD.S(kk,:);
                if kk>=ind_t05&kk<=ind_tf_9950
                    % 3
                    p3.XData = FitResult.x_opt{kk}/2/pi;
                    p3.YData = FitResult.y_opt{kk};
                    % 4
                    p4.XData = FitResult.x_opt{kk}/2/pi;
                    p4.YData = FitResult.phi{kk};
                end
                % save the frame
                ii = ii+1;
                myframe = getframe(myfig);
                im{ii} = frame2im(myframe);
            end
             
            % save vedios
            for idx = 1:loops
                [A,map] = rgb2ind(im{idx},256);
                if idx == 1
                    imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",0.01);
                else
                    imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",0.01);
                end
            end
        end
        
        % 5. AnimateEPSD - log
        function obj = animateLogEPSD(obj,Type,filename)
            if ~isfield(obj,'eEPSD')
                GMi = obj.GMi;
                TMWSEPar.TIME_Windows_f = 4;TMWSEPar.TIME_Windows_t = 4;
                TMWSEPar.L = ceil((TMWSEPar.TIME_Windows_f/GMi.dt)/2)*2+1;
                TMWSEPar.Lt = ceil((TMWSEPar.TIME_Windows_t/GMi.dt)/2)*2+1;
                TMWSEPar.K = 5;  TMWSEPar.dt = GMi.dt;
                [PHI, PHIn, sgt_2, S,w] = TMWSE(GMi.eq,TMWSEPar);
                eEPSD.PHI = PHI;eEPSD.sgt_2 = sgt_2;eEPSD.PHIn = PHIn;
                eEPSD.wk = (0:1:(TMWSEPar.L-1)/2)*2*pi/(GMi.dt*TMWSEPar.L);
                eEPSD.S = S;eEPSD.w = w;eEPSD.L = numel(w);
                eEPSD.TMWSEPar = TMWSEPar;
                obj.eEPSD = eEPSD;
            end

            % % prepare data for real EPSD
            fs_win = (0:obj.eEPSD.TMWSEPar.L/2)/(obj.eEPSD.TMWSEPar.TIME_Windows_f);
            sh = (obj.eEPSD.TMWSEPar.L-1)/2; eq_pad=padarray(obj.GMi.eq,[sh 0]);
            t_pad = (0:numel(eq_pad)-1)*obj.GMi.dt;
            teq = (0:obj.GMi.L-1)*obj.GMi.dt;
            w = obj.eEPSD.w;  tw = (0:numel(w)-1)*obj.GMi.dt;
            % % prepare data for simulated EPSD
            K_smooth = 0;K = obj.GMi.L+K_smooth;K = max(K,1000);CUT_off_FREQ = 25;
            wk = linspace(0,CUT_off_FREQ*2*pi,floor(K/2)+rem(K,2));
            SGMPar.Tmodulate.type = 'spline'; SGMPar.Tmodulate.s_time = 0;
            SGMPar = FitSGMPar_Tmodulate(obj.GMi,SGMPar);
            SGMPar.Fmodulate.wk = wk; SGMPar.Fmodulate.type = Type;
            SGMPar.Fmodulate.Paramtericmodel = 'Linearmodel';
            [SGMPar,FitResult] = FitSGMPar_Fmodulate(obj.GMi,SGMPar,obj.eEPSD);
            obj.eEPSD.FitResult = FitResult;

            % % Initial plot
            myfig = figure('OuterPosition',[100,100,800,400]);
            tilefig = tiledlayout(2,4,"TileSpacing","tight","Padding","tight");
            nexttile([1,4])
            plot(t_pad-2,eq_pad)
            ylim([-max(abs(eq_pad)),max(abs(eq_pad))]);
            xlabel('$t(s)$');ylabel('$\ddot{x}_g(t)$') 
            hold on; yyaxis right
            p1 = plot(tw-2,w,'LineWidth',1.5);
            ylabel('$w_{i}(t)$');   hold on
            scatter([0,max(teq)],[1,1]*max(w)*1.05,'Marker','v','MarkerEdgeColor','k','MarkerFaceColor','g')
            TIME_Windows_t = obj.eEPSD.TMWSEPar.TIME_Windows_t;
            xlim([min(t_pad)-TIME_Windows_t/2,max(t_pad)-TIME_Windows_t/2])
            grid on; 
            
            nexttile([1,2]);p2 = plot(tw,obj.eEPSD.S(1,:));
            xlabel('$t(s)$'); ylabel('$\ddot{x}_g(t) \times w_{i}(t)$')
            grid on; box on;
            ylim(max(max(abs(obj.eEPSD.S)))*[-1,1])
            
            nexttile([1,2])
            ind_t05 = obj.GMi.ID_p(2);         % index of time at 0.05% of Arias Intensity
            ind_tf_9950 = obj.GMi.ID_p(end-1); % index of time at 99.50% of Arias Intensity
            p3 = plot(FitResult.x_opt{ind_t05}/2/pi,log10(FitResult.y_opt{ind_t05}));
            hold on
            p4 = plot(FitResult.x_opt{ind_t05}/2/pi,log10(FitResult.phi{ind_t05}));
            grid on
            xlabel('$f(Hz)$');ylabel('PSD$_i(f)$');  
            legend('empircal','fitted')
            for nn = ind_t05:ind_tf_9950
                 m1(nn) = max(FitResult.y_opt{nn});m2(nn) = max(FitResult.phi{nn});      
            end    
            ylim([max([log10(m1);log10(m2)],[],'all')-3,max([log10(m1);log10(m2)],[],'all')]);
            xlim([0,25])
            % % Generate the fame
            % Create file name variable 
%             frameNo = [1:10:length(obj.GMi.eq),numel(obj.GMi.eq)];
            frameNo  =  [1,floor(linspace(ind_t05,ind_tf_9950,200)),obj.GMi.L];
            loops = numel(frameNo);
            ii = 0;im = [];   
            for kk = frameNo
                % update the data
                % 1
                p1.XData = tw-2+(kk-1)*obj.GMi.dt;
                % 2
                p2.YData = obj.eEPSD.S(kk,:);
                if kk>=ind_t05&kk<=ind_tf_9950
                    % 3
                    p3.XData = FitResult.x_opt{kk}/2/pi;
                    p3.YData = log10(FitResult.y_opt{kk});
                    % 4
                    p4.XData = FitResult.x_opt{kk}/2/pi;
                    p4.YData = log10(FitResult.phi{kk});
                end
                % save the frame
                ii = ii+1;
                myframe = getframe(myfig);
                im{ii} = frame2im(myframe);
            end
             
            % save vedios
            for idx = 1:loops
                [A,map] = rgb2ind(im{idx},256);
                if idx == 1
                    imwrite(A,map,filename,"gif","LoopCount",Inf,"DelayTime",0.01);
                else
                    imwrite(A,map,filename,"gif","WriteMode","append","DelayTime",0.01);
                end
            end
        end
        

        % 6. Plot constant-ductlity spectrum
        function obj = PlotNLSa(obj,mu_Sa,zeta,T_t)
            % Input: mu_Sa, target ductlity constant
            g = 9.816;

            [Ay,Sa,E_Y,Neq] = ...
                             Compute_NonL_Sa(obj.GMi,g,T_t,mu_Sa,zeta);
            
            obj.NLSA.Ay = Ay; obj.NLSA.Sa = Sa; obj.NLSA.T_t = T_t;
        
%             figure(Position=[100,100,500,400])
            for nn = 1:numel(mu_Sa)
                lstr{nn+1} = ['$\mu=$',num2str(mu_Sa(nn))];
            end
            lstr{1} = ['$\mu=$',num2str(1)];
            plot(T_t,Sa)
            hold on
            plot(T_t,Ay)
            xlabel('T [s]')
            ylabel('$A_y$')
            xlim([0.05,max(T_t)])
            grid on
            set(gca,'XScale','log')
            % set(gca,'YScale','log')
            legend(lstr,'Location','best')
            
        end

        % % 7. Compute a n-F Boucwen Structure
        function obj = plotBWmodel(StruModelPar)
            if  isnumeric(StruModelPar)
                    StruModelPar.k = 3.0*1e8;     % N/m
                    StruModelPar.m = 5*1e6;       % kg
                    StruModelPar.kexi = 0.05;     % Rayleigh damping
                    DoF = numel(StruModelPar.m);
                    StruModelPar.DoF = DoF;
                    StruModelPar.Alpha = 0.1;StruModelPar.Uy = 0.1;
                    StruModelPar.N = 5;StruModelPar.A = 1;
                    StruModelPar.Beta = 1./(2*StruModelPar.Uy.^StruModelPar.N);  % The symbol differs from the paper 'AL-GP'
                    StruModelPar.Gamma = StruModelPar.Beta;
            end
            obj.BWResult = StoryShear_nF_BoucWen(StruModelPar,obj.GMi);     

% 
%             figure
%             plot(obj.BWResult.dy,obj.BWResult.StoryShearF)
% 
%             figure(Position=[100,100,1000,300])
%             plot(obj.BWResult.t,obj.BWResult.Disp)
%             grid on

        end
        
        % % 7. Plot simulated Specturm (mean and cofident bound)
        function obj = plotSimSa(obj,SGMPar,T_t)
            if isempty(T_t)
                T_t = 0.01:0.01:1.5;T_t = [T_t 1.6:0.1:10.1];
            end
            Option.T_t =T_t;
            Option.NofSGM = 100;
            Option.ComputeSa = 'on';
            % [simGMi,t_q_pad] = SGMgeneration(GMi,SGMPar,Option);
            [simGMi,~]  = SGMgeneration(SGMPar,obj.GMi.dt,obj.GMi.teq(end),Option);
            simSa=[];
            for nn = 1:Option.NofSGM
                simSa(nn,:) = simGMi(nn).Sa;
            end
            ID = randperm(Option.NofSGM,3);
            p1 = plot(T_t,log(simSa(ID,:)),'--','Color',[6,6,6]/10,'LineWidth',1.5);
            hold on
            mu_log_Sa = mean(log(simSa));
            std_log_Sa = std(log(simSa));
            p2 = plot(T_t,mu_log_Sa,'--','Color','r','LineWidth',2);
            hold on
            g = 9.816;

            p3 = plot(obj.SA.T_t,log(obj.SA.Sa),'b-','LineWidth',1.5);
            hold on
            p4 = plot(T_t,[mu_log_Sa+std_log_Sa;mu_log_Sa-std_log_Sa],'g-.','LineWidth',1);
            hold on
            p5 = plot(T_t,[mu_log_Sa+2*std_log_Sa;mu_log_Sa-2*std_log_Sa],'m-.','LineWidth',1);
            % set(gca,'Yscale','log')
            % set(gca,'Xscale','log')
%             legend([p2(1),p4(1),p5(1),p3(1),p1(1)],'simulation mean: $\mu$','$\mu\pm{1\cdot\sigma}$',...
%                 '$\mu\pm{2\cdot\sigma}$','real','trajectories of simulation','location','best')
            grid on
            set(gca,'Xscale','log')
            xlabel('$T$ [s]')
            ylabel('$\log(S_a)$')
%             xlim([0.05,10.1])
%             Filter_cut_t = 
%             title(['$f_c=$',num2str(Filter_cut_t/2/pi),'Hz'])

        end

        % % 7-1. Plot simulated Specturm (mean and cofident bound)
        function obj = plotLinearSimSa(obj,SGMPar,T_t)
            if isempty(T_t)
                T_t = 0.01:0.01:1.5;T_t = [T_t 1.6:0.1:10.1];
            end
            Option.T_t =T_t;
            Option.NofSGM = 100;
            Option.ComputeSa = 'on';
            % [simGMi,t_q_pad] = SGMgeneration(GMi,SGMPar,Option);
            [simGMi,~]  = SGMgeneration(SGMPar,obj.GMi.dt,obj.GMi.teq(end),Option);
            simSa=[];
            for nn = 1:Option.NofSGM
                simSa(nn,:) = simGMi(nn).Sa;
            end
            ID = randperm(Option.NofSGM,10);
            p1 = plot(T_t,(simSa(ID,:)),':','Color',[9,9,9]/10,'LineWidth',1.5);
            hold on
            mu_Sa = mean((simSa));
            std_Sa = std((simSa));
%             q16 = quantile(simSa,16/100);
%             q84 = quantile(simSa,84/100);
%             q84 = quantile(simSa,2.5/100);
%             q84 = quantile(simSa,97.5/100);
            p2 = plot(T_t,mu_Sa,'--','Color','r','LineWidth',2);
            hold on
            g = 9.816;

            if ~isfield(obj.SA,'Sa')
                GMi = obj.GMi;g = 9.816;
                zeta = 0.05;
                [Sd,Sv,Sa]  = ResponseSpectraNM(GMi.eq,GMi.dt,g,T_t,zeta,0,0,1/2,1/6);
                obj.SA.T_t = T_t; obj.SA.zeta = zeta;
                obj.SA.Sa = Sd; obj.SA.Sv = Sv; obj.SA.Sa = Sa;
            end
            p3 = plot(obj.SA.T_t,(obj.SA.Sa),'b-','LineWidth',1.5);
            hold on
            p4 = plot(T_t,[quantile(simSa,84/100);quantile(simSa,16/100)],'g-.','LineWidth',1);
            hold on
            p5 = plot(T_t,[quantile(simSa,97.5/100);quantile(simSa,2.5/100)],'m-.','LineWidth',1);
            % set(gca,'Yscale','log')
%             set(gca,'Xscale','log')
%             legend([p2(1),p4(1),p5(1),p3(1),p1(1)],'simulation mean: $\mu$','$\mu\pm{1\cdot\sigma}$',...
%                 '$\mu\pm{2\cdot\sigma}$','real','trajectories of simulation','location','best')
            grid on
            set(gca,'Xscale','log')
            xlabel('$T$ [s]')
            ylabel('$S_a(T)$')
            xlim([0.05,10.1])
%             Filter_cut_t = 
%             title(['$f_c=$',num2str(Filter_cut_t/2/pi),'Hz'])

        end

        % % 8. Comput the time reaching given p% on the AI Q-Q plot
        function t_AI_p = GetSplineNodes(obj,AI_p)
            % Compute the time instant 'GM.t_AI_p' that reaches p% of the total Arial Intensity
            AI = cumtrapz(obj.GMi.teq,obj.GMi.eq.^2);
            AI_n = AI/AI(end)*100;
            [~,ID_repeat] = unique(AI_n);      
            t_AI_p = interp1(AI_n(ID_repeat),obj.GMi.teq(ID_repeat),AI_p);            
%             for nn = 1:numel(AI_p)
%                 ID_p(nn) = find(AI_n>=AI_p(nn),1,'first');
%             end     
        end

         % % 8. Comput the time reaching given p% on the AI Q-Q plot
        function [t_AI_p,t_q_pad,qmod] = PlotQQAI(obj,SGMPar)
            % Compute the time instant 'GM.t_AI_p' that reaches p% of the total Arial Intensity
            AI_p = SGMPar.Tmodulate.Par.p;
            AI = cumtrapz(obj.GMi.teq,obj.GMi.eq.^2);
            AI_n = AI/AI(end)*100;
            [~,ID_repeat] = unique(AI_n);      
            t_AI_p = interp1(AI_n(ID_repeat),obj.GMi.teq(ID_repeat),AI_p);                      
            plot(obj.GMi.teq,AI_n*AI(end)/100,'-b','LineWidth',2)
            hold on
            plot(t_AI_p,AI_p*AI(end)/100,'o','MarkerEdgeColor','k',...
                'MarkerFaceColor','g','LineWidth',0.5)
            hold on
            [t_q_pad,qmod,padn] = TModulate_Estimator(SGMPar.Tmodulate.Par,...
                                obj.GMi.dt,obj.GMi.L);
            plot(t_q_pad,cumtrapz(t_q_pad,qmod.^2),'r--','LineWidth',1.5)
            ylabel('$I_a(t)$');xlabel('$t$');
        end

        % % 9
        function E_lowfeq_fn = PlotLowFre(obj,fn)
            GMi = obj.GMi;
            fs = 1/GMi.dt;
            Y = fft(GMi.eq);f_real = fs*(0:(GMi.L/2))/GMi.L;
            P2 = abs(Y/GMi.L);P1_real = P2(1:GMi.L/2+1);
            P1_real(2:end-1) = 2*P1_real(2:end-1);
            obj.FFT.P1 = P1_real;
            obj.FFT.f = f_real;

            E_f = cumtrapz(f_real,P1_real.^2);
            E_f_n = E_f/E_f(end)*100;
            E_lowfeq_fn = interp1(f_real,E_f_n,fn);

        end

    end



    methods (Static)

        % 1. Write the SGMPar into matrix 
        function XED = getSGMPar(SGMPar_aug)
            NofGM = size(SGMPar_aug,2);
            for ir = 1:NofGM
                SGMPar = SGMPar_aug{ir};
                Wf(ir,1) = SGMPar.HighPassFilter.Wf;
                AI(ir,1) = SGMPar.Tmodulate.Par.AI(end);
                t_AI_p(ir,:) = SGMPar.Tmodulate.Par.t_AI_p(1:end);
                wg_mid(ir,1)  = SGMPar.Fmodulate.Par.wg_tmid;
                wg_rate(ir,1) = SGMPar.Fmodulate.Par.wg_rate;
                zeta_mid(ir,1)  = SGMPar.Fmodulate.Par.zeta_tmid;
                zeta_rate(ir,1) = SGMPar.Fmodulate.Par.zeta_rate;
            end
            AI_log = log(AI);
            tmid = t_AI_p(:,4);
            D5_95 = t_AI_p(:,6)-t_AI_p(:,2);
            % Assemble Pars
            XED = [AI_log,wg_mid,wg_rate,zeta_mid,zeta_rate,tmid,D5_95,Wf/2/pi];
        end
      
    end

end
