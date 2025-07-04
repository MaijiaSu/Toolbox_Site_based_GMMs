function simGMs = SiteSGMGenerator(SiteGMM,ScenVar,NofsimGMs,dt)
% Default Paramter
% h = 6;
% Ztor = 0;
% f_fltz = Ztor;
% f_fltz(Ztor<1) = 1;
% f0  = 0; % srike-slip

%%
beta = SiteGMM.GMPEs.beta;
Sigma2_T = SiteGMM.GMPEs.Sigma2_T;
err_corrMat = SiteGMM.GMPEs.err_corrMat;
x0 = SiteGMM.GMPEs.basisfuns(ScenVar);


% 1.  GMPE mean
Z_mean =  beta*x0(:);

% 2. GMPE residual
Ndim = size(beta,1);
CoV_ZZ = diag(sqrt(Sigma2_T))*err_corrMat*diag(sqrt(Sigma2_T));
Z_error =  mvnrnd(zeros(Ndim,1),CoV_ZZ,NofsimGMs);

% 3. GMPE
GMPE_pre = repmat(Z_mean,1,NofsimGMs)' + Z_error;

% 4. Transfrom back orignal space
for nn = 1:SiteGMM.Ndim
    CDF_XED_pre(:,nn) = normcdf(GMPE_pre(:,nn),0,1);
end
X_GMPE = uq_all_invcdf(CDF_XED_pre, SiteGMM.Pars_JointPDF.UqLabInput.Marginals);


%% Bound the GMPE prediction
% Loop through each dimension and clamp values to the specified bounds
DisBounds = SiteGMM.Pars_JointPDF.DisBounds;
for nn = 1:Ndim
    min_val = DisBounds(1,nn);  % Minimum bound for this dimension
    max_val = DisBounds(2,nn);  % Maximum bound for this dimension
    % Clamp values to the min and max bounds for this dimension
    id1 = X_GMPE(:,nn)<min_val;
    X_GMPE(id1, nn) = min_val;
    id2 = X_GMPE(:,nn)>max_val;
    X_GMPE(id2, nn) = max_val;
end

%% Generate GMs using predicted GM parmeters
NofWN = 1;
NofComp = SiteGMM.NofComp;
ndim =  SiteGMM.Ndim/NofComp;
for nn = 1:NofComp
    simGMpars{nn} = X_GMPE(:,(nn-1)*ndim+1:nn*ndim);
end
if NofComp >= 1
    simGMs.H1 = SGMGenerator_simPar(SiteGMM.FitOption,simGMpars{1},dt);
elseif  NofComp >= 2
    simGMs.H2 = SGMGenerator_simPar(SiteGMM.FitOption,simGMpars{2},dt);
elseif  NofComp >= 3
    simGMs.V = SGMGenerator_simPar(SiteGMM.FitOption,simGMpars{3},dt);
end
simGMs.dt = dt;


% simGM_GMPE_marjor = SGMGenerator_simPar(SiteGMM.FitOption,GMpar_marjor,dt);
% simGM_GMPE_inter  = SGMGenerator_simPar(SiteGMM.FitOption,GMpar_inter,dt);
% simGM_GMPE_V = SGMGenerator_simPar(SiteGMM.FitOption,GMpar_V,dt);
% % %% Compute RotD50 Sa
% % Tt = logspace(log10(5e-2),log10(10),101);
% % zeta = 0.05;
% % RotD50_Sa=[];
% % dt = 0.02;
% % for ir = 1:NofGMPEsample
% %     ag1 = simGM_GMPE_marjor{ir}.eq;
% %     ag2 = simGM_GMPE_inter{ir}.eq;
% %     RotD50_Sa(ir,:) = ComputeRotD50_Sa(ag1,ag2,Tt,zeta,dt);
% % end
% %
% % %% V/H ratio of simulated motions
% % g = 9.816;
% % for nn = 1:NofGMPEsample
% %     [~,~,Sa_V(nn,:)] =  ResponseSpectraNM(simGM_GMPE_V{ir}.eq,dt,g,Tt,zeta,0,0,1/2,1/6);
% % end
% % sim_VH = Sa_V./RotD50_Sa;
% %
% % %% Sa of NGA-GMPE
% % z1 = 999;
% %
% % % gmpe_bssa_2014(M, T, Rjb, Fault_Type, region, z1, Vs30
% % Fault_Type=1;
% % % Fault_Type    = 0 for unspecified fault
% % %               = 1 for strike-slip fault
% % %               = 2 for normal fault
% % %               = 3 for reverse fault
% % region=0;
% % % region        = 0 for global (incl. Taiwan)
% % %               = 1 for California
% % %               = 2 for Japan
% % %               = 3 for China or Turkey
% % %               = 4 for Italy
% % T0 = logspace(log10(5e-2),log10(10),12);
% % % GMPE#1
% % Rjb = r0;Rx=0;Ry0 = 0;Ztor=0;delta=45; lambda = 0;fas=1; HW=0;W = 999;Z10 = 999;FVS30=1;
% % [mu_GMPE_1, std_GMPE_1] = ASK_2014_nga(m0, Tt, r0, Rjb, Rx, Ry0, Ztor, delta, lambda, fas, HW, W, Z10, vs30, FVS30, region);
% % [mu01, std01] = ASK_2014_nga(m0, T0, r0, Rjb, Rx, Ry0, Ztor, delta, lambda, fas, HW, W, Z10, vs30, FVS30, region);
% % % GMPE#2
% % [mu_GMPE_2, std_GMPE_2] = gmpe_bssa_2014(m0, Tt, r0, Fault_Type, region, z1, vs30);
% % [mu02, std02] = gmpe_bssa_2014(m0, T0, r0, Fault_Type, region, z1, vs30);
% % % GMPE#3
% % Rjb = r0;Rx=0;Ztor=0;delta=45; lambda = 0;fas=1; HW=0;W = 999;Ztor=999;Zbot=nan;Fhw=0;Z25=exp(7.089 - 1.144*log(vs30));Zhyp=9;region=0;
% % [mu_GMPE_3, std_GMPE_3] = CB_2014_nga(m0, Tt, r0, Rjb, Rx, W, Ztor, Zbot, delta, lambda, Fhw, vs30, Z25, Zhyp, region);
% % [mu03, std03] = CB_2014_nga(m0, T0, r0, Rjb, Rx, W, Ztor, Zbot, delta, lambda, Fhw, vs30, Z25, Zhyp, region);
% % % GMPE#4
% % Rjb = r0;
% % [mu_GMPE_4, std_GMPE_4] = CY_2014_nga(m0, Tt, r0, Rjb, Rx, Ztor, delta, lambda, Z10, vs30,Fhw, FVS30, region);
% % [mu04, std04] = CY_2014_nga(m0, T0, r0, Rjb, Rx, Ztor, delta, lambda, Z10, vs30,Fhw, FVS30, region);
% %
% % %% Compare Sa #1: horizontal
% % close all
% % selfGrootDefault(0)
% %
% % figure('Position',[100,100,1000,450],'Visible','on')
% % tt = tiledlayout(2,2,"TileSpacing","tight","Padding","tight")
% % nexttile(1,[2,1])
% % p1 = plot(T0,mu01,'v','Color',[4,4,4]/10);
% % hold on
% % p2 = plot(T0,mu02,'s','Color',[4,4,4]/10);
% % hold on
% % p3 = plot(T0,mu03,'^','Color',[4,4,4]/10);
% % hold on
% % p4 = plot(T0,mu04,'o','Color',[4,4,4]/10);
% % hold on
% % ALLGMPE_mean = exp(mean([log(mu_GMPE_1);log(mu_GMPE_2);log(mu_GMPE_3);log(mu_GMPE_4)]));
% % p5 = plot(Tt,ALLGMPE_mean,'-','Color',[4,4,4]/10);
% % hold on
% % ALLGMPE_std = sqrt(mean([std01.^2;std02.^2;std03.^2;std04.^2]));
% % ALLGMPE_mean = exp(mean([log(mu01);log(mu02);log(mu03);log(mu04)]));
% % yneg = exp(log(ALLGMPE_mean) - 1*ALLGMPE_std);
% % ypos = exp(log(ALLGMPE_mean) + 1*ALLGMPE_std);
% % p6 = errorbar(T0,ALLGMPE_mean,yneg,ypos);
% % p6.Color = [4,4,4]/10;
% % set(gca,'XScale','log','Yscale','log')
% % xlim([0.04,12])
% %
% % % Plot Sa of simulated motions
% % logmu_sim = mean(log(RotD50_Sa));
% % logstd_sim = std(log(RotD50_Sa));
% % p7=plot(Tt,exp(logmu_sim),'--r');
% % hold on
% % plot(Tt,exp(logmu_sim+1*logstd_sim),'--r')
% % hold on
% % plot(Tt,exp(logmu_sim-1*logstd_sim),'--r')
% % set(gca,'XScale','log','Yscale','log')
% % legend([p1,p2,p3,p4,p6,p7],'ASK14: median','BSSA14: median','CB14: median',...
% %     'CY14: medan','NGA west2 GMPEs: $\mu$, $\pm \sigma$',...
% %     '200 site-based motions: $\mu$, $\mu \pm \sigma$',...
% %     'Location','southwest')
% %
% % grid on
% % title(['$M=',num2str(m0),',','R=',num2str(r0),'km',',','Vs30=',num2str(vs30),'$'])
% % xlabel('$T$ [s]')
% % ylabel('RotD50 $Sa(T)$')
% %
% %
% % %% Compare Sa #2: Vertical
% % % GMPE
% % FRV=0;FNM=1;
% % [mu2_VH, sig_logVH, ~, ~] = gmpeVH_ga_2011(m0, r0, vs30, FRV, FNM, Tt);
% %
% % nexttile(2,[1,1])
% % p1 = plot(Tt,mu2_VH,'Color',[4,4,4]/10)
% % hold on
% % hold on
% % set(gca,'XScale','log','Yscale','linear')
% % hold on
% % median_sim_VH = median(sim_VH);
% % plot(Tt,median_sim_VH,'r')
% % ylim([0,1.5])
% % grid on
% % title(['$M=',num2str(m0),',','R=',num2str(r0),'km,','Vs30=',num2str(vs30),'$'])
% % xlabel('$T$ [s]')
% % ylabel('median $V/H$')
% %
% % nexttile(4,[1,1])
% % p1 = plot(Tt,sig_logVH,'Color',[4,4,4]/10)
% % hold on
% % std_logVH_sim = std(log(sim_VH));
% % plot(Tt,std_logVH_sim,'r')
% % ylabel('$\sigma[\log(V/H)]$')
% % grid on
% % xlabel('$T$ [s]')
% % set(gca,'XScale','log','Yscale','linear')
% % legend('GMPE: GA2011','Simulated motions',...
% %     'Location','best')
% % ylim([0,1])