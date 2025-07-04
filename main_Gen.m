%% Generate ground motion for specifc earthquake scenarios
vs30 = 760;m0 = 6.5;r0 = 10;
Ftype = 0;ztor = 0;
ScenVar = [m0,r0,vs30,Ftype,Ztor];
dt = 0.02;
NofsimGMs = 100;
simGMs = SiteSGMGenerator(SiteGMM,ScenVar,NofsimGMs,dt);