function SGMPar = SGMParFitting(GMi,FitOption)
% %
if nargin==1
    FitOption = struct();
    FitOption.Tmodulate = struct();
    FitOption.Fmodulate = struct();
    FitOption.HighPassFilter = struct();
end
if ~isfield(FitOption,'GenMethod')
    FitOption.GenMethod = 'spectral';
end


%% ------------------------------------------------------------------------
%                            Input Paramets
%% ------------------------------------------------------------------------
%% Input Paramets
% 1. Seimic wave
% Input GMi as a struct variables, including following fileds
GMi.dt = GMi.dt;
GMi.eq = GMi.eq(:);
GMi.L = GMi.L;

% 2. Pars setting for TMWSE
TMWSEPar.TIME_Windows_f = 4;                                   % Default time windows for short-time FFT
TMWSEPar.TIME_Windows_t = 4;                                   % Default time windows for smoothing q(t)
TMWSEPar.L = ceil((TMWSEPar.TIME_Windows_f/GMi.dt)/2)*2+1;     % Number of discrete point for the windows choose an odd number
TMWSEPar.Lt = ceil((TMWSEPar.TIME_Windows_t/GMi.dt)/2)*2+1;    % Number of discrete point for time avaraging
TMWSEPar.K = 5;                                                % Default number of leakage-resistant orthogonal windows
TMWSEPar.dt = GMi.dt;                                          % Integration time step

% 3. Pars setting for Ground Motion Model
SGMPar = struct();

% 3.0 Generation method
SGMPar.GenMethod = FitOption.GenMethod;

% 3.1 Time-modulating function
DefaultPar = inputParser;
DefaultPar.KeepUnmatched = 1;
DefaultPar.addParameter('type','spline');       % 'spline','Gamma'
DefaultPar.addOptional('s_time',0,@isnumeric);  % length of zeros-padding, unit:[s]
parse(DefaultPar,FitOption.Tmodulate);
SGMPar.Tmodulate = DefaultPar.Results;

% 3.2 Time-frequency-modulating function phi(t_i,w|theta)
% 3.2.1 Discretization (for estimating Pars of spectral nonstaionary)
% scheme#1 (Mayssa, Marco)
% CUT_off_FREQ = 25;
% dw = 2*pi/GM.teq(end);
% K = ceil(CUT_off_FREQ*2*pi/dw);
% K = max(K,1000);
% wk = linspace(0,CUT_off_FREQ*2*pi,K);
% dw = wk(2)-wk(1);
% scheme#2
K_smooth = ceil(SGMPar.Tmodulate.s_time/GMi.dt);
K = GMi.L+K_smooth;
K = max(K,1000);
CUT_off_FREQ = 25;
wk = linspace(0,CUT_off_FREQ*2*pi,floor(K/2)+rem(K,2));
dw = wk(2)-wk(1);

% 3.2.2 Methods for capturing the characteristic of spectral non-stationaritly
DefaultPar = inputParser;
DefaultPar.KeepUnmatched = 1;
DefaultPar.addParameter('Paramtericmodel','Linearmodel'); % 'Constantmodel' ,'Linearmodel','Polymodel';
DefaultPar.addParameter('type','IIorderFilter'); %
parse(DefaultPar,FitOption.Fmodulate);
SGMPar.Fmodulate = DefaultPar.Results;
SGMPar.Fmodulate.dw = dw;
SGMPar.Fmodulate.wk = wk;

% 4. Setting for high-pass filter
DefaultPar = inputParser;
DefaultPar.KeepUnmatched = 1;
DefaultPar.addParameter('Type','Oscillator'); % 'Oscillator' and 'Butterworth'
parse(DefaultPar,FitOption.HighPassFilter);
SGMPar.HighPassFilter = DefaultPar.Results;

%% ------------------------------------------------------------------------
%                    Check the Parameter setting
%% ------------------------------------------------------------------------
if strcmp(SGMPar.GenMethod,'temporal')&&...
        ~strcmp(SGMPar.Fmodulate.type,'IIorderFilter')
    error('Unmatched filter model for temporal representation methods')
end

%% ------------------------------------------------------------------------
%                               Output
%% ------------------------------------------------------------------------
% 1. Estimated EPSD
% 2.   Pars of Synthetic Ground Motion (SGM)
% 2.1. Pars of temporal modulating function
% 2.2. Pars of spectral modulating function
% 2.3. Pars of high-pass filter, i.e., conner frequency

%% ------------------------------------------------------------------------

%                               Start fitting

%% ------------------------------------------------------------------------

%% ------------------------------------------------------------------------
%                         Step.0 Initialization
%% ------------------------------------------------------------------------
% Data Preparation
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


%% ------------------------------------------------------------------------
%                 Step.1 Compute the empirical EPSD
%% ------------------------------------------------------------------------
[PHI, PHIn, sgt_2, S,w] = TMWSE(GMi.eq,TMWSEPar);
eEPSD.PHI = PHI;
eEPSD.sgt_2 = sgt_2;
eEPSD.PHIn = PHIn;
eEPSD.wk = (0:1:(TMWSEPar.L-1)/2)*2*pi/(GMi.dt*TMWSEPar.L);
eEPSD.S = S;
eEPSD.Hannwin = w;
eEPSD.L = numel(w);


%% ------------------------------------------------------------------------
%                  Step.2 Estimate the Pars of SGM
%% ------------------------------------------------------------------------
%% 2.1. Compute the Pars of the Temporal nonstationary Mdeol
% input: GM,  model type
% output:
% 1. Estimated Pars of the given model ('spline','Gamma',...)
% 2. and return the Time modulating function
disp('Perform temporal nonstationary fitting')
SGMPar = FitSGMPar_Tmodulate(GMi,SGMPar);
% Remark:
% i.time-modulating functions can be smoothed by a hann-window
% ii. [t_q_pad,qmod,padn] = TModulate_Estimator(SGMPar.Tmodulate.Par,teq),
% which returns a discrte time-modulating function

%% 2.2. Compute the pars of the Spectral nonstationary Model
% input: GM, eEPSD, model_type
% output:
% 1. Estimated Pars of the given spectral-nonstationary model (functional types of EPSD and parameteric model for estimated pars)
disp('Perform spectral nonstationary fitting')
SGMPar = FitSGMPar_Fmodulate(GMi,SGMPar,eEPSD);

%% 2.3. Compute the conner frequency
disp('Perform conner frequency fitting')
SGMPar.HighPassFilter.Wf = FitSGMPar_HighPassFilter(SGMPar,GMi);

% (eEPSD,SGMPar,EPSD_q,dw,wk,t_q_pad,GM)
% ,EPSD_q,dw,wk
end