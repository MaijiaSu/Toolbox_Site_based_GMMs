clc,clear,close all
addpath('Site_based_GMM_toolbox')
set(0,'DefaultFigureColor',[1 1 1])
set(groot, 'defaultAxesTickLabelInterpreter','latex');
set(groot, 'defaultLegendInterpreter','latex');
set(groot, 'defaultTextInterpreter','latex');
set(0,'defaulttextinterpreter','latex','DefaultAxesFontSize',12)
load GMDataset.mat

%% Fitting Option of Stochastic Ground Motion Models
% I. Simulation Methods
FitOption.GenMethod = 'spectral';                     % 'spectral' or 'temporal'
% II. Time-modulating functions
FitOption.Tmodulate.type = 'spline';                  % 'spline', 'Gamma'
% III. Frequency filters
FitOption.Fmodulate.type = 'IIorderFilter';           % 'IIorderFilter','KanaiTajimi','ConvexII-2','CascadeII-2','CloughPenzien'
% IV. Trend functions of slected filters
FitOption.Fmodulate.Paramtericmodel = [2,1];          %  type of trend functions, a vector input with size consitent with free paramter of the filter type
                                                      %  = 1, costant;
                                                      %  = 2, linear model;
                                                      %  = 3, smoothed polyline    
% V. High-pass filters                                                   
FitOption.HighPassFilter.Type = 'Oscillator';         % 'Oscillator' or 'Butterworth'


%% Step 1. Fit the GMM Paramters
SiteGMM = Fit_SGMPars(FitOption,GM);

%% Step 2. Fit the Ground Motion Prediction Equations
% load GMDataset.mat
% load SiteGMM_607NGArecords
M = [MetaData.M];
R_RUP = [MetaData.Clstd];
VS30 = [MetaData.VS30];
EQID = [MetaData.EQID];
SW_Rake = [MetaData.Rake];
id1 = (logical(SW_Rake>=60&SW_Rake<=120)); % Reverse
id2 = (logical((SW_Rake>=-30&SW_Rake<=30)+(SW_Rake>=150&SW_Rake<=270))); % Strike-slip
Ftype = ones(1,SiteGMM.NofGM);Ftype(id2) = 0;
Ztor = [MetaData.Ztor];
ScenVar = [M',R_RUP',VS30',Ftype',Ztor'];

initial_explainTerm = ...
[1	1	1	1	1	1	1
1	1	0	1	1	0	1
1	1	0	1	1	0	1
1	1	0	1	1	0	1
1	1	0	1	1	0	1
1	1	0	1	1	0	1
1	1	0	1	1	0	1
1	1	0	0	1	0	1
1	1	0	0	1	0	0
1	1	0	0	1	0	0
1	1	0	1	1	0	1];

basisfuns = @(ScenVar) GMPE_BasisFuns(ScenVar);
SiteGMM = Fit_GMPEs(SiteGMM,basisfuns,ScenVar,EQID,initial_explainTerm);

%% Save data
save FitResults SiteGMM

