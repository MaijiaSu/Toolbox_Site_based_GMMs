function Filter_cut_t = FitSGMPar_HighPassFilter(SGMPar,GMi)
%% Compute the conner frequency of the filter
% Input: i. the full-nonstationary EPSD
%        ii. the discrete size in frequency domain.i.e., wk=k*dw
%        ii. the discrete sequence in time dimain, i.e., t_eq
% Output: the estimated cornner frequency

method = 4;
%% ------------------------------------------------------------------------
%       Method#1:based on energy fitting (from Broccardo and Dabaghi)
%% ------------------------------------------------------------------------
if method == 1
    
    % 0. Compare data
    dw = SGMPar.Fmodulate.dw;
    wk = SGMPar.Fmodulate.wk; 
     [t_q_pad,qmod,padn]  = TModulate_Estimator(SGMPar.Tmodulate.Par,GMi.dt,GMi.L);
    % Compute the unit-variance EPSD
    EPSD_unit_var =  FModulateComputeEPSD(SGMPar.Tmodulate,SGMPar.Fmodulate,GMi.L,GMi.dt,wk);

    % unit-variance EPSD with padding
    EPSD_unit_var_pad = [repmat(EPSD_unit_var(:,1),1,padn),...
        EPSD_unit_var,...
        repmat(EPSD_unit_var(:,end),1,padn)];
    % Full-nonstationary EPSD
    EPSD_q = EPSD_unit_var_pad.*repmat(qmod.^2,numel(wk),1);
    
    % 1. Initialization
    %     w=PHI.dw:PHI.dw:PHI.K*PHI.dw;
    %     PHI.psih=0.9;
    v_eq_ref = cumsum(GMi.eq)*GMi.dt;
    u_eq_ref = cumsum(cumsum(GMi.eq))*GMi.dt^2;
    en_v = cumsum(v_eq_ref.^2)*GMi.dt;
    en_u = cumsum(u_eq_ref.^2)*GMi.dt;
    Wf0 = dw*4;
    Phi_curr= Wf0;

    M_2 = ComputeM_2(Wf0,t_q_pad,EPSD_q,wk);

    % figure('Name','new')
    % [X1,X2] = meshgrid(t_q_pad,wk);
    % contour(X1,X2,log(M_2),50)

    dts = GMi.dt;
    M=sqrt(M_2);
    M_0 = M;%
    M_2_0 = M_2;
    M_d=diff([zeros(size(M,1),1), M],1,2)/(dts);       % Derivative of M
    M_dd=diff([zeros(size(M,1),1), M_d],1,2)/(dts);    % Derivative of M
    %WW=repmat(w,[size(M,2) 1])';
    WW_2=repmat(wk.^2,[size(M,2) 1])';
    M_2_d = M_d.^2+WW_2.*M.^2;
    M2_dd = (M_dd-WW_2.*M).^2+4*WW_2.*M_d.^2;  % Theoretical EVPSD

    v_q_2 = cumsum(2*sum(M_2_d)*dw)*dts;
    u_q_2 = cumsum(2*sum(M_2)*dw)*dts;

    % 2. Compute fc iteratively
    fun_curr = abs(en_u(end)-u_q_2(end));
    ii =0;
    while abs(fun_curr)/u_q_2(end) > 0.01
        %     if fun_curr/u_q_2(end) < 0.01
        %         break;
        %     end
        ii = ii+1;
        % Gradient
        Wf = Phi_curr+dw/5;
        %     [~,~,M_2] = EPSD_phi_newV2(qmod,PHI);

        M_2 = ComputeM_2(Wf,t_q_pad,EPSD_q,wk);

        u_q_2 = cumsum(2*sum(M_2)*dw)*dts;
        fun_curr_eps = abs(en_u(end)-u_q_2(end));
        Grad = (fun_curr_eps-fun_curr)/(dw/5);
        %
        Phi_new = Phi_curr-fun_curr/Grad;
        %
        Phi_curr = Phi_new;
        Wf = Phi_curr;

        %     [~,EPSD_q,M_2] = EPSD_phi_newV2(qmod,PHI);
        M_2 = ComputeM_2(Wf,t_q_pad,EPSD_q,wk);
        u_q_2 = cumsum(2*sum(M_2)*dw)*dts;
        fun_curr = abs(en_u(end)-u_q_2(end));
    end


    %
    Filter_cut_t = Wf;
    if Filter_cut_t < 0
        Filter_cut_t = Wf0 ;
        disp('WARNING the filter has a negative value!')
        M_2 = M_2_0;
        % M = M_0;
    end

end
%% ------------------------------------------------------------------------
%             Method#4:based on Spectrum compatiable -v4 （mim -mean(k)）
%% ------------------------------------------------------------------------

if method == 4

    Option.NofSGM = 100;
    Option.ComputeSa = 'on';
    Option.GenMethod = 'spectral';               % spectral or temporal 
    Option.Sa_Zeta = 0.05;                       % Sa damper ratio 
    
   % --- Preprocess_dataset
       CUT_off_FREQ = 25;  % so that dt=0.02s
    [acc,dt] = MyGMdownsampling(CUT_off_FREQ,GMi.dt,GMi.eq);
    cut_p_lower = 0.1; % unit：%
    cut_p_upper = 99.9; % unit：%
    acc = MyGMcut_quasizeros(cut_p_lower,cut_p_upper,dt,acc);
    GMi.eq = acc;
    GMi.dt = dt;
    GMi.L = numel(acc);
    GMi.teq = (0:GMi.L-1)*dt;
    SGMPar = FitSGMPar_Tmodulate(GMi,SGMPar);
    % ----------------------------------

    % Sa-compatiable region T>T_cut = 1s
    T_cut = 1;

    % Sa for real Ground motion
%     T_t = 0.01:0.01:1.5;
%     T_t = [T_t 1.6:0.1:10.1];
%     T_t = linspace(T_cut,10,50);
    T_t = logspace(log10(T_cut),log10(10),30);
    g = 9.816;
    [~,~,Sa]  = ResponseSpectraNM(GMi.eq,GMi.dt,g,T_t,0.05,0,0,1/2,1/6);
    GMi.Sa = Sa;
    
    % loop#1
    loop1_n = 11;
    WWF1 = linspace(0,2,loop1_n)*2*pi;
    for n_wf = 1:numel(WWF1)
        % simulate SGM and compute the Sa
        SGMPar.HighPassFilter.Wf = WWF1(n_wf);
        Option.T_t = T_t;
%         Option.GenMethod = FitOption.GenMethod; 
        
%         [simGMi,~] = SGMgeneration(GMi,SGMPar,Option);
        [simGMi,~]  = SGMgeneration(SGMPar,GMi.dt,GMi.teq(end),Option);
        for nn = 1:Option.NofSGM
            simSa(nn,:) = simGMi(nn).Sa;
        end

        % compute the confident bound
        mu_Sa = mean(log(simSa));
        std_Sa = std(log(simSa));

        % compute the bias
        b_T = log(GMi.Sa)-mu_Sa;

        % compute the relative bias
%         k = sign(b_T).*((b_T)./std_log_Sa).^2;
          k = ((b_T)./std_Sa);

        % compute the mean of k 
%         id_cut = find(T_t>=T_cut,1,'first');
        k1_mean(n_wf) = abs(trapz(T_t,k));
    end

    % loop#2
    [~,id0] = min(k1_mean);
    if id0==1
        id1 = 1;
        id2 = 3;
        loop2_n = 11;
    elseif id0==loop1_n
        id1 = loop1_n-2;
        id2 = loop1_n;
        loop2_n = 11;
    else
        id1 = id0-1;
        id2 = id0+1;
        loop2_n = 21;
    end   
    wf1 = WWF1(id1);
    wf2 = WWF1(id2);
    WWF2 = linspace(wf1,wf2,loop2_n);
    for n_wf = 1:numel(WWF2)
        % simulate SGM and compute the Sa
        SGMPar.HighPassFilter.Wf = WWF2(n_wf);

        [simGMi,~] = SGMgeneration(SGMPar,GMi.dt,GMi.teq(end),Option);
        for nn = 1:Option.NofSGM
            simSa(nn,:) = simGMi(nn).Sa;
        end

        % compute the confident bound
        mu_Sa = mean(log(simSa));
        std_Sa = std(log(simSa));

        % compute the bias
        b_T = log(GMi.Sa)-mu_Sa;

        % compute the k
%         k = sign(b_T).*((b_T)./std_log_Sa).^2;
          k = ((b_T)./std_Sa);
        
%         id_cut = find(T_t>=T_cut,1,'first');
          k2_mean(n_wf) =  abs(trapz(T_t,k));
    end
   
    % best-estimation
    [~,id] = min(k2_mean);
    Filter_cut_t = WWF2(id);
end  
%% ------------------------------------------------------------------------
%             Method#5:based on Spectrum compatiable -v4 （mim -mean(k)）
%% ------------------------------------------------------------------------

if method == 5

    Option.NofSGM = 100;
    Option.ComputeSa = 'on';
    Option.GenMethod = 'spectral';               % spectral or temporal 
    Option.Sa_Zeta = 0.05;                       % Sa damper ratio 
    
    % Sa-compatiable region T>T_cut = 1s
    T_cut = 5;

    % Sa for real Ground motion
%     T_t = 0.01:0.01:1.5;
%     T_t = [T_t 1.6:0.1:10.1];
%     T_t = linspace(T_cut,10,50);
    T_t = logspace(log10(T_cut),log10(30),30);
    g = 9.816;
    [~,~,Sa]  = ResponseSpectraNM(GMi.eq,GMi.dt,g,T_t,0.05,0,0,1/2,1/6);
    GMi.Sa = Sa;
    
    % loop#1
    loop1_n = 11;
    WWF1 = linspace(0,2,loop1_n)*2*pi;
    for n_wf = 1:numel(WWF1)
        % simulate SGM and compute the Sa
        SGMPar.HighPassFilter.Wf = WWF1(n_wf);
        Option.T_t = T_t;
%         Option.GenMethod = FitOption.GenMethod; 
        
%         [simGMi,~] = SGMgeneration(GMi,SGMPar,Option);
        [simGMi,~]  = SGMgeneration(SGMPar,GMi.dt,GMi.teq(end),Option);
        for nn = 1:Option.NofSGM
            simSa(nn,:) = simGMi(nn).Sa;
        end

        % compute the confident bound
        mu_Sa = mean(log(simSa));
        std_Sa = std(log(simSa));

        % compute the bias
        b_T = log(GMi.Sa)-mu_Sa;

        % compute the relative bias
%         k = sign(b_T).*((b_T)./std_log_Sa).^2;
          k = ((b_T)./std_Sa);

        % compute the mean of k 
%         id_cut = find(T_t>=T_cut,1,'first');
        k1_mean(n_wf) = abs(trapz(T_t,k));
    end

    % loop#2
    [~,id0] = min(k1_mean);
    if id0==1
        id1 = 1;
        id2 = 3;
        loop2_n = 11;
    elseif id0==loop1_n
        id1 = loop1_n-2;
        id2 = loop1_n;
        loop2_n = 11;
    else
        id1 = id0-1;
        id2 = id0+1;
        loop2_n = 21;
    end   
    wf1 = WWF1(id1);
    wf2 = WWF1(id2);
    WWF2 = linspace(wf1,wf2,loop2_n);
    for n_wf = 1:numel(WWF2)
        % simulate SGM and compute the Sa
        SGMPar.HighPassFilter.Wf = WWF2(n_wf);

        [simGMi,~] = SGMgeneration(SGMPar,GMi.dt,GMi.teq(end),Option);
        for nn = 1:Option.NofSGM
            simSa(nn,:) = simGMi(nn).Sa;
        end

        % compute the confident bound
        mu_Sa = mean(log(simSa));
        std_Sa = std(log(simSa));

        % compute the bias
        b_T = log(GMi.Sa)-mu_Sa;

        % compute the k
%         k = sign(b_T).*((b_T)./std_log_Sa).^2;
          k = ((b_T)./std_Sa);
        
%         id_cut = find(T_t>=T_cut,1,'first');
          k2_mean(n_wf) =  abs(trapz(T_t,k));
    end
   
    % best-estimation
    [~,id] = min(k2_mean);
    Filter_cut_t = WWF2(id);
end

%% ------------------------------------------------------------------------
%             Method#4:based on Spectrum compatiable -v4 （linear scale）
%% ------------------------------------------------------------------------

if method == 41

    Option.NofSGM = 100;
    Option.ComputeSa = 'on';
    Option.GenMethod = 'spectral';               % spectral or temporal 
    Option.Sa_Zeta = 0.05;                       % Sa damper ratio 
    
    % Sa-compatiable region T>T_cut = 1s
    T_cut = 1;

    % Sa for real Ground motion
%     T_t = 0.01:0.01:1.5;
%     T_t = [T_t 1.6:0.1:10.1];
%     T_t = linspace(T_cut,10,50);log
    T_t = logspace(T_cut,10,101);
    g = 9.816;
    [~,~,Sa]  = ResponseSpectraNM(GMi.eq,GMi.dt,g,T_t,0.05,0,0,1/2,1/6);
    GMi.Sa = Sa;
    
    % loop#1
    loop1_n = 11;
    WWF1 = linspace(0,4,loop1_n)*2*pi;
    for n_wf = 1:numel(WWF1)
        % simulate SGM and compute the Sa
        SGMPar.HighPassFilter.Wf = WWF1(n_wf);
        Option.T_t = T_t;
%         Option.GenMethod = FitOption.GenMethod; 
        
%         [simGMi,~] = SGMgeneration(GMi,SGMPar,Option);
        [simGMi,~]  = SGMgeneration(SGMPar,GMi.dt,GMi.teq(end),Option);
        for nn = 1:Option.NofSGM
            simSa(nn,:) = simGMi(nn).Sa;
        end

        % compute the confident bound
        mu_Sa = mean((simSa));
        std_Sa = std((simSa));

        % compute the bias
        b_T = (GMi.Sa)-mu_Sa;

        % compute the relative bias
%         k = sign(b_T).*((b_T)./std_log_Sa).^2;
          k = ((b_T)./std_Sa);

        % compute the mean of k 
%         id_cut = find(T_t>=T_cut,1,'first');
        k1_mean(n_wf) = abs(trapz(T_t,k));
    end

    % loop#2
    [~,id0] = min(k1_mean);
    if id0==1
        id1 = 1;
        id2 = 3;
        loop2_n = 11;
    elseif id0==loop1_n
        id1 = loop1_n-2;
        id2 = loop1_n;
        loop2_n = 11;
    else
        id1 = id0-1;
        id2 = id0+1;
        loop2_n = 21;
    end   
    wf1 = WWF1(id1);
    wf2 = WWF1(id2);
    WWF2 = linspace(wf1,wf2,loop2_n);

    for n_wf = 1:numel(WWF2)
        % simulate SGM and compute the Sa
        SGMPar.HighPassFilter.Wf = WWF2(n_wf);

        [simGMi,~] = SGMgeneration(SGMPar,GMi.dt,GMi.teq(end),Option);
        for nn = 1:Option.NofSGM
            simSa(nn,:) = simGMi(nn).Sa;
        end

        % compute the confident bound
        mu_Sa = mean((simSa));
        std_Sa = std((simSa));

        % compute the bias
        b_T = (GMi.Sa)-mu_Sa;

        % compute the k
%         k = sign(b_T).*((b_T)./std_log_Sa).^2;
          k = ((b_T)./std_Sa);
        
%         id_cut = find(T_t>=T_cut,1,'first');
          k2_mean(n_wf) =  abs(trapz(T_t,k));
    end
   
    % best-estimation
    [~,id] = min(k2_mean);
    Filter_cut_t = WWF2(id);
end  

%% Figure
%     plot_fit_optfc(WWF1,WWF2,k1_mean,k2_mean,Filter_cut_t,SGMPar,GMi)

  

end