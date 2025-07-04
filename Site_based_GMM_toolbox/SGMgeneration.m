function [simGM,t_q_pad,K_SGM] = SGMgeneration(SGMPar,dt,tf,Option)

%% Defacult option
DefaultPar = inputParser; 
DefaultPar.KeepUnmatched = 1;
DefaultPar.addParameter('ComputeSa','off');
DefaultPar.addParameter('TextShow','off');
DefaultPar.addOptional('Sa_Zeta',0.05,@isnumeric);
DefaultPar.addOptional('NofSGM',1,@isnumeric);
% DefaultPar.addOptional('GenMethod','spectral'); % spectral or temporal
DefaultPar.addOptional('RVs',1,@isnumeric); % Methods of generating random seed
T_t = 0.01:0.01:1.5; T_t = [T_t 1.6:0.1:10.1];
DefaultPar.addOptional('T_t',T_t,@isnumeric);
parse(DefaultPar,Option);
Option = DefaultPar.Results;

if ~isfield(SGMPar,'GenMethod')
    SGMPar.GenMethod = 'spectral'; %temporal
end

if ~isfield(SGMPar.HighPassFilter,'Type')
    SGMPar.HighPassFilter.Type = 'Oscillator';
end

if ~isfield(SGMPar.Tmodulate,'s_time')
    if strcmp(SGMPar.Tmodulate.type,'spline')
        SGMPar.Tmodulate.s_time = 0;
    else
        SGMPar.Tmodulate.s_time = 3;
    end
end


%% Discretization of the simulatede GM

% Default value of tf only if tf =[];
if isempty(tf)
    if strcmp(SGMPar.Tmodulate.type,'Gamma')
        tmid = SGMPar.Tmodulate.Par.tmid;
        D5_95 = SGMPar.Tmodulate.Par.D5_95;
        if 2*D5_95>tmid
            tf = 2*D5_95;
        else
            tf = D5_95+tmid;
        end
    elseif strcmp(SGMPar.Tmodulate.type,'spline')
        tf = SGMPar.Tmodulate.Par.t_AI_p(end);
        t_AI_p = SGMPar.Tmodulate.Par.t_AI_p;
D5_95= t_AI_p(6)- t_AI_p(2);
    end
end

if isnumeric(tf) == 1
    L = ceil(tf/dt)+1;
else ischar(tf) == 1
    L = str2num(tf);  
end
s_time = SGMPar.Tmodulate.s_time;
K_smooth = ceil(s_time/dt);
K = L+K_smooth;
CUT_off_FREQ = 25; 
wk  = linspace(0,CUT_off_FREQ*2*pi,floor(K/2));
% K = 2*numel(wk);
dw = wk(2)-wk(1);
% teq = (0:L-1)*dt;

wk  = linspace(0,CUT_off_FREQ*2*pi,floor(K/2)+1);
wk(1) = [];

%% Check the Lower bound of fc
fc = SGMPar.HighPassFilter.Wf/2/pi;
if strcmp(SGMPar.Tmodulate.type,'Gamma')
    D5_95 = SGMPar.Tmodulate.Par.D5_95;
elseif strcmp(SGMPar.Tmodulate.type,'spline')
    t_AI_p = SGMPar.Tmodulate.Par.t_AI_p;
    D5_95= t_AI_p(6)- t_AI_p(2);
end
fc_min = 1/D5_95;
fc_min = 0;
fc = max([fc,fc_min]);
SGMPar.HighPassFilter.Wf = fc*2*pi;

%% Prepare the parameters
% % 1 Time-modulating function
% SGMPar.Tmodulate.Par
% qmod = SGMPar.Tmodulate.fun(SGMPar.Tmodulate.Par,GM.teq);
if ~isfield(SGMPar.Tmodulate.Par,'type')
    SGMPar.Tmodulate.Par.type = SGMPar.Tmodulate.type;
end
[t_q_pad,qmod,padn] = TModulate_Estimator(SGMPar.Tmodulate.Par,dt,L);

% % 2 Spectral-modulating function (Only required for Temporal simulation)
% Par = Fmodulate_Estimator(Paramtericmodel,BN,Par,GM,Tmodulate)
if strcmp(SGMPar.GenMethod,'temporal')  
    [~,fitted_opt_par] =  FModulateComputeEPSD(SGMPar.Tmodulate,SGMPar.Fmodulate,L,dt,wk);
    wg_ti = fitted_opt_par(:,1);
    zeta_ti = fitted_opt_par(:,2);
    % Bound the parameters
    wg_ti(wg_ti<=25/500*2*pi) = 25/500*2*pi;
    zeta_ti(zeta_ti<=0.02)=0.02;
    zeta_ti(zeta_ti>=0.98)=0.98;
    % Pading
    wg_ti_pad = [ones(padn,1)*wg_ti(1);wg_ti;ones(padn,1)*wg_ti(end)];
    zeta_ti_pad = [ones(padn,1)*zeta_ti(1);zeta_ti;ones(padn,1)*zeta_ti(end)];
end

% % 3. Conner Frequency
Wf = SGMPar.HighPassFilter.Wf;
fc = Wf/2/pi;

%% 4. Compute the full-nonstationary Matrix K_SGM
% Remark: Both the temporal&spectral representaion methods to generate SGM
% can be fomulated as the form: Xg = K_SGM * W, where W is a discrte white
% noise vector
% Option.GenMethod = 'temporal';
% -------------------------------------------------------------------------
% % 4-1: temporal representaion method
if strcmp(SGMPar.GenMethod,'temporal')
    % Pseudo Acceleration IRF
    h = @(wf,kexif,dt) wf/sqrt(1-kexif^2)*exp(-kexif*wf*dt).*...
        sin(wf*sqrt(1-kexif)*dt);
    H = zeros(K,K);
    for nn = 1:K
        t_temp = dt*(1:K-nn+1);
        H(nn:K,nn) = h(wg_ti_pad(nn),zeta_ti_pad(nn),t_temp);
    end
    % Normalization to ensure variance == 1
    tempNorm = sqrt(sum(H.^2,2));
    H2 = H./repmat(tempNorm,1,K);
    % Multplication with the time-modulating function
    K_SGM_b = H2.*repmat(qmod',1,K);

    % Pass K_SGM_b through high-pass filter
    if Wf>0
        if strcmp(SGMPar.HighPassFilter.Type,'Oscillator')
            hf = t_q_pad.*exp(-Wf.*t_q_pad);
            Hf = repmat(hf',1,K);
            K_SGM_D  =  fastConv(K_SGM_b,Hf,1);
            K_SGM_D  =  real(K_SGM_D(1:K,:))*dt;
            K_SGM_V  =  diff([zeros(1,K);K_SGM_D],1,1)/dt;
            K_SGM  =  diff([zeros(1,K);K_SGM_V],1,1)/dt;
        elseif strcmp(SGMPar.HighPassFilter.Type,'Butterworth')
            Fs = 1/dt;  % Nyquist frequency
            cutoff = fc/(Fs/2); % cut-off ratio
            if ~isfield(SGMPar.HighPassFilter,'order')
                order = 2;
            else
                order = SGMPar.HighPassFilter.order;
            end
            [b,a] = butter(order,cutoff,'high');
            K_SGM = filter(b,a,K_SGM_b,[],1);
        end
    else %i.e, if Wf<=0, high-pass filter is not used
        K_SGM = K_SGM_b;
    end
    
    % Compute the energy-correction factor
    q_b = sum(K_SGM_b.^2,2);
    q_a = sum(K_SGM.^2,2);
    corrE = trapz(t_q_pad,q_b)/trapz(t_q_pad,q_a);
    K_SGM = sqrt(corrE)*K_SGM;

% -------------------------------------------------------------------------
% % 4-2: spectral representaion method
elseif strcmp(SGMPar.GenMethod,'spectral')
    % Compute the unit-variance EPSD
    EPSD_unit_var =  FModulateComputeEPSD(SGMPar.Tmodulate,SGMPar.Fmodulate,L,dt,wk);

    % unit-variance EPSD with padding
    EPSD_unit_var_pad = [repmat(EPSD_unit_var(:,1),1,padn),...
        EPSD_unit_var,...
        repmat(EPSD_unit_var(:,end),1,padn)];

%     temp = sum(sqrt(EPSD_unit_var_pad*dw),1);
%     temp = sum(EPSD_unit_var_pad,1);
% temp = sum(sqrt(EPSD_unit_var_pad),1);

    % Full-nonstationary EPSD
    EPSD_q = EPSD_unit_var_pad.*repmat(qmod.^2,numel(wk),1);
    % Final Generation Martix K_SGM_b (before applying high-pass filter)
    corrE = 1;
    s1 = cos(t_q_pad'*wk).*sqrt(2*corrE*abs(EPSD_q)'*dw);
    s2 = sin(t_q_pad'*wk).*sqrt(2*corrE*abs(EPSD_q)'*dw);
    K_SGM_b = [s1,s2];

    % High-pass filter
    if  Wf>0
        if strcmp(SGMPar.HighPassFilter.Type,'Oscillator')  % critically damped oscillator
            % Method I: Priestley theory (discard)
            % Energy coorection factor
%             [corrE,~,~,M2_dd_n] = ComputeCorrE(Wf,EPSD_q,dw,wk,t_q_pad,dt);
%             % Final matrix K_SGM
%             s1 = cos(t_q_pad'*wk).*sqrt(2*corrE*abs(M2_dd_n)'*dw);
%             s2 = sin(t_q_pad'*wk).*sqrt(2*corrE*abs(M2_dd_n)'*dw);
%             K_SGM = [s1,s2];
            % Method II:Discretization+Convolution
            hf     =   t_q_pad.*exp(-Wf.*t_q_pad);
            Apf    =   K_SGM_b;
            Bpf    =   repmat(hf(:),[1,numel(wk)*2]).*dt;
            Cpf    =   fastConv(Apf,Bpf,1);
            K_SGM_D  =  Cpf(1:numel(t_q_pad),:);
            K_SGM_V  =  diff([zeros(1,numel(wk)*2);K_SGM_D],1,1)/dt;
            K_SGM = diff([zeros(1,numel(wk)*2);K_SGM_V],1,1)/dt;
        elseif strcmp(SGMPar.HighPassFilter.Type,'Butterworth')
            Fs = 1/dt;  % Nyquist frequency
            cutoff = fc/(Fs/2); % cut-off ratio
            if ~isfield(SGMPar.HighPassFilter,'order')
                order = 2;
            else
                order = SGMPar.HighPassFilter.order;
            end
            [b,a] = butter(order,cutoff,'high');
            K_SGM = filter(b,a,K_SGM_b,[],1);
        end
    else %i.e, if Wf<=0, high-pass filter is non-used
        K_SGM = K_SGM_b;
    end

    % Enery correction
    q_b = sum(K_SGM_b.^2,2);
    q_a = sum(K_SGM.^2,2);
    corrE = trapz(t_q_pad,q_b)/trapz(t_q_pad,q_a);
%     corrE = 1;
    K_SGM = K_SGM*sqrt(corrE);
end

% save data
% corrE = 1;
% SGMPar.HighPassFilter.corrE = corrE;

%%  5. Generate synthetic ground motion and Compute Sa
simGM = [];
for ir = 1:Option.NofSGM
    % 5.1 Generate seismic acceleration
    if strcmp(Option.TextShow,'on')
        disp(['Synthetic Ground motion Simulating:',num2str(ir),'/',num2str(Option.NofSGM)]);
    end
    
    if Option.RVs == 1  % Random seed#1
        U2 = normrnd(0,1,[size(K_SGM,2),1]);
    elseif Option.RVs == 2  % Random seed#2
%         cos(a+b) = cos(a)cos(b) - sin(a)sin(b)
        phi =  rand(size(wk,2),1)*2*pi-pi;
        U1cos = sqrt(2)*cos(wk(:).*phi); U1sin = -sqrt(2)*sin(wk(:).*phi);
        U2 = [U1cos;U1sin];
    end
    a_g = K_SGM*U2;
    a_g = a_g(:)';
%     v_g = cumtrapz(t_q_pad,a_g);
%     u_g = cumtrapz(t_q_pad,v_g);

    % 5.2 Compute the Spectrum accleration
    if strcmp(Option.ComputeSa,'on')
        g = 9.816;
        zeta = Option.Sa_Zeta;
        [~,~,Sa,maxT]  = ResponseSpectraNM(a_g,dt,g,Option.T_t,zeta,0,0,1/2,1/6);  

        simGM(ir).Sa = Sa;
%         simGM(ir).maxT = maxT;
    end
    
    %  5.3. Storage
    %     simGM(ir).u_g = u_g;
    %     simGM(ir).v_g = v_g;
    simGM(ir).eq = a_g;
    simGM(ir).L = numel(a_g);
    simGM(ir).corrE = corrE;
    simGM(ir).fc = fc;
    simGM(ir).dt = dt;
end