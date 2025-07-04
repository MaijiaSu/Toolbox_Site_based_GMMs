clc,clear,close all
load FitResults

%% #1: Generate ground motion using real GM parameters
dt = 0.02
simGMs = SiteSGMGenerator(1,SiteGMM,[],[],dt);


%% #2: Generate ground motion similar to a dataset
dt = 0.02;
NofsimGMs = 100;
simGMs = SiteSGMGenerator(2,SiteGMM,[],NofsimGMs,dt);
% given a dataset
% simGMs = SiteSGMGenerator(SiteGMM,ScenVar,NofsimGMs,dt)

%% #3: Generate ground motion for specifc earthquake scenarios
% Given a scenarious
vs30 = 760;m0 = 6.5;r0 = 10;
Ftype = 0;ztor = 0;
ScenVar = [m0,r0,vs30,Ftype,ztor];
dt = 0.02;
NofsimGMs = 100;
simGMs = SiteSGMGenerator(3,SiteGMM,ScenVar,NofsimGMs);


