function [ALLGMPE_mean,ALLGMPE_std,mu01,mu02,mu03,mu04,std01,std02,std03,std04] = ...
                                           NGA_GMPE(vs30,m0,r0,Tt)

% Tt = logspace(log10(5e-2),log10(10),101);
% vs30 
% m0 
% r0

z1 = 999;
% gmpe_bssa_2014(M, T, Rjb, Fault_Type, region, z1, Vs30
Fault_Type=1;
% Fault_Type    = 0 for unspecified fault
%               = 1 for strike-slip fault
%               = 2 for normal fault
%               = 3 for reverse fault
region=0;
% region        = 0 for global (incl. Taiwan)
%               = 1 for California
%               = 2 for Japan
%               = 3 for China or Turkey
%               = 4 for Italy

T0 = Tt;

% GMPE#1
Rjb = r0;Rx=0;Ry0 = 0;Ztor=0;delta=45; lambda = 0;fas=1; HW=0;W = 999;Z10 = 999;FVS30=1;
% [mu_GMPE_1, std_GMPE_1] = ASK_2014_nga(m0, Tt, r0, Rjb, Rx, Ry0, Ztor, delta, lambda, fas, HW, W, Z10, vs30, FVS30, region);
[mu01, std01] = ASK_2014_nga(m0, T0, r0, Rjb, Rx, Ry0, Ztor, delta, lambda, fas, HW, W, Z10, vs30, FVS30, region);
% GMPE#2
% [mu_GMPE_2, std_GMPE_2] = gmpe_bssa_2014(m0, Tt, r0, Fault_Type, region, z1, vs30);
[mu02, std02] = gmpe_bssa_2014(m0, T0, r0, Fault_Type, region, z1, vs30);
% GMPE#3
Rjb = r0;Rx=0;Ztor=0;delta=45; lambda = 0;fas=1; HW=0;W = 999;Ztor=999;Zbot=nan;Fhw=0;Z25=exp(7.089 - 1.144*log(vs30));Zhyp=9;region=0;
% [mu_GMPE_3, std_GMPE_3] = CB_2014_nga(m0, Tt, r0, Rjb, Rx, W, Ztor, Zbot, delta, lambda, Fhw, vs30, Z25, Zhyp, region);
[mu03, std03] = CB_2014_nga(m0, T0, r0, Rjb, Rx, W, Ztor, Zbot, delta, lambda, Fhw, vs30, Z25, Zhyp, region);
% GMPE#4
Rjb = r0;
% [mu_GMPE_4, std_GMPE_4] = CY_2014_nga(m0, Tt, r0, Rjb, Rx, Ztor, delta, lambda, Z10, vs30,Fhw, FVS30, region);
[mu04, std04] = CY_2014_nga(m0, T0, r0, Rjb, Rx, Ztor, delta, lambda, Z10, vs30,Fhw, FVS30, region);


% Version#1
% ALLGMPE_mean = mean([mu_GMPE_1;mu_GMPE_2;mu_GMPE_3;mu_GMPE_4]);
% ALLGMPE_std = sqrt(mean([std01.^2;std02.^2;std03.^2;std04.^2]));

% Version#2
% ALLGMPE_mean = exp(mean([log(mu_GMPE_1);log(mu_GMPE_2);log(mu_GMPE_3);log(mu_GMPE_4)]));
% ALLGMPE_std = sqrt(mean([std01.^2;std02.^2;std03.^2;std04.^2]));
% ALLGMPE_mean = exp(mean([log(mu01);log(mu02);log(mu03);log(mu04)]));
            

% Dabaghi 2018
ALLGMPE_mean=exp(1/4*(log(mu01))+1/4*(log(mu02))+1/4*(log(mu03))+1/4*(log(mu04))); %unEQUAL WEIGHTS
ALLGMPE_1std=exp(1/4*(log(mu01)+std01)+1/4*(log(mu02)+std02)+1/4*(log(mu03)+std03)+2/4*(log(mu04)+std04)); %unEQUAL WEIGHTS
ALLGMPE_std=log(ALLGMPE_1std)-log(ALLGMPE_mean);
