function DisBounds = GetSGMPars_support(SiteGMM)


FitOption = SiteGMM.FitOption;
ndim = SiteGMM.Ndim/SiteGMM.NofComp;
NofComp = SiteGMM.NofComp;


DisBounds = zeros(2,ndim);

% log_AI
DisBounds(:,1) = [-inf;inf];

% Pars of Time-modulating functions
if strcmp(FitOption.Tmodulate.type,'spline')
    nPar_T = 7;
    DisBounds(:,2:7) = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1;
        20, 15,  10,  20, 40, 40];
elseif strcmp(FitOption.Tmodulate.type,'Gamma')
     nPar_T = 3;
    DisBounds(:,2:3) = [0.5, 5;
        40, 45];
end


% Pars of Frequency-modulating functions
[~,FilterNofPar,par_type] = ...
                    SetFilterModel(FitOption.Fmodulate.type);
Paramtericmodel = FitOption.Fmodulate.Paramtericmodel;
FrePar_support = [];
for nn = 1:FilterNofPar-1
    switch Paramtericmodel(nn)
        case 1 % constant model
            switch par_type(nn)
                case 1 % wg
                    tempbound = [0.1*2*pi;inf];
                case 2 % zeta_g
                    tempbound = [0.02;1];
                case 3 % alpha
                    tempbound = [0;1];
            end
        case 2 % linear model
            switch par_type(nn)
                case 1 % wg
                    tempbound = [0.1*2*pi,-inf;
                                inf,inf];
                case 2 % zeta_g
                    tempbound = [0.02,-inf;
                                  1,inf];
                case 3 % alpha
                    tempbound = [0,-inf;
                                  1,inf];
            end
        case 3 % nonparametic model
            switch par_type(nn)
                case 1 % wg
                    tempbound = [0.1*2*pi,0.1*2*pi,0.1*2*pi,0.1*2*pi,0.1*2*pi;
                                inf,inf,inf,inf,inf];
                case 2 % zeta_g
                    tempbound = [0.02,0.02,0.02,0.02,0.02;
                                  1,1,1,1,1];
                case 3 % alpha
                     tempbound = [0,0,0,0,0;
                                  1,1,1,1,1];
            end
    end
    FrePar_support = [FrePar_support,tempbound];
end
DisBounds(:,nPar_T+1:end-1) = FrePar_support;

% Pars of high-pass filter
DisBounds(:,end) = [0.01,2];

% repeat the supports for the NofComp compnents
DisBounds = repmat(DisBounds,1,NofComp);
