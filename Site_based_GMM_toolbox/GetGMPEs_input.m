function [XX,Z_corr,explainTerm] = GetGMPEs_input(MetaData)

EQID
ScenVar = [M,Ftype,Ztor,R_RUP,VS30]

%% Explanatory variables of GMPE
NofGM = size(MetaData,2);
M = [MetaData.M];
R_RUP = [MetaData.Clstd];
VS30 = [MetaData.VS30];
EQID = [MetaData.EQID];
SW_Rake = [MetaData.Rake];
id1 = (logical(SW_Rake>=60&SW_Rake<=120)); % Reverse
id2 = (logical((SW_Rake>=-30&SW_Rake<=30)+(SW_Rake>=150&SW_Rake<=270))); % Strike-slip
Ftype = ones(1,NofGM);
Ftype(id2) = 0;
Ztor = [MetaData.Ztor];
f_fltz = Ztor;
f_fltz(Ztor<1) = 1;
h = 6; % km, focal depth
VS30_hat = min([VS30;1100*ones(1,NofGM)]);

% inter-event
EQID = [MetaData.EQID];
uniqueEQID = unique(EQID);
NofEQ = numel(uniqueEQID);
% construct correlation matrix
R_InterCorr = zeros(NofGM);
for nn = 1:NofEQ
    eqID = uniqueEQID(nn);
    id = find(EQID==eqID);
    for ii = 1:numel(id)
        R_InterCorr(id(ii),setdiff(id,id(ii))) = 1;
    end
end

Z_corr = zeros(NofGM,NofEQ);
for nn = 1:NofEQ
    Z_corr(EQID==uniqueEQID(nn),nn) = 1;
end

% Fault type
SW_Rake = [MetaData.Rake];
% Reverse
id1 = (logical(SW_Rake>=60&SW_Rake<120));
% Strike-slip
id2 = (logical((SW_Rake>=-30&SW_Rake<30)+(SW_Rake>=150&SW_Rake<270)));
% reverse‐oblique
id3 = (logical((SW_Rake>=120&SW_Rake<150)+(SW_Rake>=30&SW_Rake<60)));
[sum(id1),sum(id2),sum(id3)];

EQID = [MetaData.EQID];
uniqueEQID = unique(EQID);
NofEQ = numel(uniqueEQID);


%%
X0 = ones(1,NofGM);
X1 = M;
X2 = (M-6.5).*(M>6.5);
X3 = Ftype.*f_fltz;
X4 = log(sqrt(R_RUP.^2+h^2));
X5 = M.*log(sqrt(R_RUP.^2+h^2));
X6 = log(VS30_hat);
XX = [X0;X1;X2;X3;X4;X5;X6]';



%%
explainTerm_H = ...
    [1	1	1	1	1	1	1
    1	1	0	0	1	0	0
    1	1	0	1	1	0	0
    1	1	0	1	0	0	1
    1	1	0	0	1	0	1
    1	1	0	0	1	0	1
    1	1	0	1	0	0	1
    1	1	0	0	0	0	1
    1	0	0	0	1	0	0
    1	1	0	0	1	0	0
    1	1	0	0	1	0	1];

explainTerm_V = ...
    [1	1	1	1	1	1	1
    1	1	0	1	1	0	0
    1	1	0	0	1	0	0
    1	1	0	0	1	0	1
    1	1	0	1	1	0	1
    1	1	0	1	1	0	1
    1	1	0	1	0	0	1
    1	1	0	0	1	0	0
    1	1	0	0	0	0	0
    1	1	0	0	1	0	0
    1	1	0	0	1	0	1];

explainTerm = [explainTerm_H;explainTerm_H;explainTerm_V];
end