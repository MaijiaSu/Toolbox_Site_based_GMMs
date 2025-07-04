function XX = GMPE_BasisFuns(ScenVar,EQID)
% ScenVar = [M,R_RUP,VS30,Ftype,Ztor]
M = ScenVar(:,1);
R_RUP  = ScenVar(:,2);
VS30  = ScenVar(:,3);
Ftype  = ScenVar(:,4);
Ztor  = ScenVar(:,5);
NofGM = size(ScenVar,1);

%% Explanatory variables of GMPE
f_fltz = Ztor;
f_fltz(Ztor<1) = 1;
h = 6; % km, focal depth
VS30_hat = min([VS30(:),1100*ones(NofGM,1)],[],2);

%%
X0 = ones(NofGM,1);
X1 = M;
X2 = (M-6.5).*(M>6.5);
X3 = Ftype.*f_fltz;
X4 = log(sqrt(R_RUP.^2+h^2));
X5 = M.*log(sqrt(R_RUP.^2+h^2));
X6 = log(VS30_hat);
XX = [X0,X1,X2,X3,X4,X5,X6];

