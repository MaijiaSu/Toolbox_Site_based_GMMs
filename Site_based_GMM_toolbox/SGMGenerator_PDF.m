function [simGM,simGMPar] = SGMGenerator_PDF(SGMPar,GMPar_PDF,NofGM,NofMCS,Tt)
% Input: 
% 1. GM model setttings, SGMPar
% 2. PDF settings of GM Pars,GMPar_PDF
% 3. NofGM: number of GMs in simulated catalog
% 4. NofMCS: Number of simulated catalog
% Output:
% 1. NofMCS synthetic Catalogs, each catalog has NofGM GM
% 2. The simulated GM features of the synthetic catalogs

%% 1. Use UQLab toolbox to construct PDF model
clear InputOpts

% 1.1 Copula Model
Copula_Type = GMPar_PDF.Copula.Type;
if strcmp(Copula_Type,'Gaussian')
    InputOpts.Copula.Type = Copula_Type; %'auto','Gaussian','Dvine','Cvine'
    InputOpts.Copula.Parameters = GMPar_PDF.Copula.Parameters;
elseif strcmp(Copula_Type,'DVine')||strcmp(Copula_Type,'CVine')
%     InputOpts.Copula.Type = Copula_Type; 
%     InputOpts.Copula.Structure = GMPar_PDF.Copula.Structure;
%     InputOpts.Copula.Families = GMPar_PDF.Copula.Families;
%     InputOpts.Copula.Rotations = GMPar_PDF.Copula.Rotations;
%     InputOpts.Copula.Parameters = GMPar_PDF.Copula.Parameters;
%     InputOpts.Copula.Truncation = GMPar_PDF.Copula.Truncation;
        InputOpts.Copula = GMPar_PDF.Copula;
end

% 1.2 Marginal Models
MarginalType = GMPar_PDF.Marginals.Types;
MarginalPar = GMPar_PDF.Marginals.Par;
DisBounds = GMPar_PDF.DisBounds;
Ndim = size(MarginalType,2);

for nn = 1:Ndim
    InputOpts.Marginals(nn).Type = MarginalType{nn};
    InputOpts.Marginals(nn).Bounds = DisBounds(:,nn)';
    InputOpts.Marginals(nn).Parameters = MarginalPar(nn,1:2);
end

%#1
% clear myInput
myInput = uq_createInput(InputOpts);

% #2: use old marginals
% load temp_InputOpts
% tempMarginals = myInput.Marginals;
% clear myInput
% myInput = uq_createInput(InputOpts);
% myInput.Marginals = tempMarginals;

% #3: use copula model
% load temp_InputOpts
% tempCopula = myInput.Copula;
% clear myInput
% myInput = uq_createInput(InputOpts);
% myInput.Copula = tempCopula;

% #4:
% clear InputOpts
% InputOpts.Marginals(1).Type = 'Gaussian';
% InputOpts.Marginals(2).Type = 'Gamma';
% InputOpts.Marginals(3).Type = 'Laplace';
% InputOpts.Marginals(4).Type = 'Beta';
% InputOpts.Marginals(5).Type = 'LogNormal';
% InputOpts.Marginals(6).Type = 'LogNormal';
% InputOpts.Marginals(7).Type = 'Exponential';
% DisBounds = [ -inf, 0, -inf, 0.02, 0.5, 5,  0;
%     inf, inf, inf,   1,  40,  45, 2];
% InputOpts.Inference.Data = GMPar_PDF.XED;
% InputOpts.Copula.Type = 'auto';
% myInput = uq_createInput(InputOpts);

%% 2. Simulate GM given the GMPar PDF model
simGM=[];
simGMPar=[];
tstart = tic;
for nofMCS = 1:NofMCS

    X_MCS = uq_getSample(myInput,NofGM*1,'MC');
    simGMPar{nofMCS} = X_MCS;

    for ir = 1:size(X_MCS,1)   
        if strcmp(SGMPar.Tmodulate.type,'spline')
            % % Parameters#1
            Par.wg_tmid = X_MCS(ir,2);
            Par.wg_rate = X_MCS(ir,3);
            Par.zeta_tmid = X_MCS(ir,4);
            Par.zeta_rate = 0;
            Filter_opt_par = [Par.wg_tmid,Par.zeta_tmid;
                Par.wg_rate,Par.zeta_rate];
            SGMPar.Fmodulate.opt_par = Filter_opt_par;
            % % Paramters#2
            SGMPar.Fmodulate.Par.tmid = [];
            AI_MCS = exp(X_MCS(ir,1));
            SGMPar.Tmodulate.Par.AI = AI_MCS;
            SGMPar.Tmodulate.Par.t_AI_p = [];
            Dt_MCS = X_MCS(ir,5:10);
            t_AI_p_MCS = [0,cumsum(Dt_MCS)];
            SGMPar.Tmodulate.Par.t_AI_p = t_AI_p_MCS;
            t_AI_p = SGMPar.Tmodulate.Par.t_AI_p;
            D5_95= t_AI_p(6)- t_AI_p(2);
            % % Parameter#3
            SGMPar.HighPassFilter.Wf = 2*pi*X_MCS(ir,11);
        elseif strcmp(SGMPar.Tmodulate.type,'Gamma')
            % % Parameters#1
            Par.wg_tmid = X_MCS(ir,2);
            Par.wg_rate = X_MCS(ir,3);
            Par.zeta_tmid = X_MCS(ir,4);
            Par.zeta_rate = 0;
            Filter_opt_par = [Par.wg_tmid,Par.zeta_tmid;
                Par.wg_rate,Par.zeta_rate];
            SGMPar.Fmodulate.opt_par = Filter_opt_par;
            % % Paramters#2
            SGMPar.Fmodulate.Par.tmid = [];
            AI_MCS = exp(X_MCS(ir,1));
            SGMPar.Tmodulate.Par.AI = AI_MCS;
             SGMPar.Tmodulate.Par.tmid = X_MCS(ir,5);
            SGMPar.Tmodulate.Par.D5_95 = X_MCS(ir,6);
            D5_95= SGMPar.Tmodulate.Par.D5_95;
            % % Parameter#3
            SGMPar.HighPassFilter.Wf = 2*pi*X_MCS(ir,7);
      
        end

        % ---- Set lower-bound to fc ---------
%         fc = SGMPar.HighPassFilter.Wf/2/pi;
%         fc_min = 1/D5_95;
%         fc = max([fc,fc_min]);
%         SGMPar.HighPassFilter.Wf = fc*2*pi;
        % ------------------------------------

        % % Start GM generation
        Option.NofSGM = 1;
        Option.ComputeSa = 'on';
        Option.T_t = Tt;
        Option.GenMethod = 'spectral';    % spectral or temporal
        Option.Sa_Zeta = 0.05;

        dt = 0.02;   
        [simGM{ir}(nofMCS)] = SGMgeneration(SGMPar,dt,[],Option);
%         [simGM{ir}(nofMCS)] = SGMgeneration(SGMPar,dt,[],Option);
        tend = toc(tstart);
        disp([num2str(nofMCS),'/',num2str(NofMCS),'--',num2str(ir),'/',num2str(size(X_MCS,1)),':','Elapsed time = ', ...
            num2str((tend)/60,'%.1f'),'mins']);
    end

end
return

