function simGM_GMPE = SGMGenerator_simPar(FitOption,X_GMPE,dt)
    NofGM = size(X_GMPE,1);
    ndim =  size(X_GMPE,2);
    if strcmp(FitOption.Tmodulate.type,'spline')
        nPar_T = 7;
        SGMPar.Tmodulate.Par.p = [0,5,33,45,75,95,100];
        SGMPar.Tmodulate.Par.t_AI_p = [];
    elseif strcmp(SGMPar_ir.Tmodulate.type,'Gamma')
        nPar_T = 3;
        SGMPar_ir.Tmodulate.Par.D5_95 = [];
        SGMPar_ir.Tmodulate.Par.tmid = [];
    end
    SGMPar.GenMethod = FitOption.GenMethod;
    SGMPar.Tmodulate.type = FitOption.Tmodulate.type;
    SGMPar.Fmodulate = FitOption.Fmodulate;
    SGMPar.Fmodulate.opt_par = [];
    SGMPar.HighPassFilter.Type = 'Oscillator';
    SGMPar.HighPassFilter.Wf = [];

    trendmodel = SGMPar.Fmodulate.Paramtericmodel;
    if ischar(trendmodel)
        if strcmp(trendmodel,'ReducedLinearmodel')
            trendmodel = [2,1];
        end
    end
    [~,FilterNofPar] = SetFilterModel(SGMPar.Fmodulate.type);
    nofTrendPars = zeros(1,FilterNofPar);
    nofTrendPars(find(trendmodel==1)) = 1;
    nofTrendPars(find(trendmodel==2)) = 2;
    nofTrendPars(find(trendmodel==3)) = 5;
    nPar_F = ndim-nPar_T-1;

    % Prepare the data
    log_AI = X_GMPE(:,1);
    Pars_T = X_GMPE(:,2:nPar_T);
    Wf = X_GMPE(:,end)*2*pi;

    for ir = 1:NofGM
        tempPars_F = X_GMPE(ir,nPar_T+1:nPar_T+nPar_F);
        temp_optF = zeros(FilterNofPar-1,max(nofTrendPars));
        for n_trend = 1:FilterNofPar-1
             temp_optF(1:nofTrendPars(n_trend),n_trend) = tempPars_F(1:nofTrendPars(n_trend))';
             tempPars_F(1:nofTrendPars(n_trend)) = [];
        end
        Pars_F{ir} = temp_optF;
    end

%% 2. Simulate GM given the GMPar PDF model
tstart = tic;

for ir = 1:NofGM
    
    % Time modulating function
    if strcmp(SGMPar.Tmodulate.type,'spline')
        AI_MCS = exp(log_AI(ir,1));
        SGMPar.Tmodulate.Par.AI = AI_MCS;

        SGMPar.Tmodulate.Par.t_AI_p = [];
        Dt_MCS = Pars_T(ir,:);
        t_AI_p_MCS = [0,cumsum(Dt_MCS)];
        SGMPar.Tmodulate.Par.t_AI_p = t_AI_p_MCS;
    elseif strcmp(SGMPar.Tmodulate.type,'Gamma')
        AI_MCS = exp(log_AI(ir,1));
        SGMPar.Tmodulate.Par.AI = AI_MCS;

        SGMPar.Tmodulate.Par.tmid = Pars_T(ir,1);
        SGMPar.Tmodulate.Par.D5_95 = Pars_T(ir,2);

    end
  
    % Frequency-modulating function
     SGMPar.Fmodulate.opt_par = Pars_F{ir};

     % High-pass filter
    SGMPar.HighPassFilter.Wf = Wf(ir);

    
%%
    % % Start GM generation
    Option.NofSGM = 1; % number of white noise samples
    Option.ComputeSa = 'off';
    [simGM_GMPE{ir}] = SGMgeneration(SGMPar,dt,[],Option);
    tend = toc(tstart);
%     disp(['GMPE-based simulator:',num2str(ir),'/',num2str(size(X_GMPE,1)),':','Elapsed time = ', ...
%         num2str((tend)/60,'%.1f'),'mins']);
end
return

