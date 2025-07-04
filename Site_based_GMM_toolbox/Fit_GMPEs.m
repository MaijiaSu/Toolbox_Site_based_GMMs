function SiteGMM = Fit_GMPEs(SiteGMM,GMPE_BasisFuns,ScenVar,EQID,initial_explainTerm)

NofGM = SiteGMM.NofGM;
Ndim = SiteGMM.Ndim;
%% Input XX in linear random-effect regression
% [XX,Z_corr,explainTerm] = GetGMPEs_input(SiteGMM.MetaData);
XX = GMPE_BasisFuns(ScenVar);

%% inter-event
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


%% Output theta_Z in linear regressions
% Ground Motion Parameters
% Fit the joint PDF of GMM parameters
SiteGMM = Fit_GMMPars_JointPDF(SiteGMM);

% Transformation of model parameters to the normal space
MagianlPDFs = SiteGMM.Pars_JointPDF.UqLabInput.Marginals;
CDF_XED = uq_all_cdf(SiteGMM.theta_GM_bounded, MagianlPDFs);
theta_Z = norminv(CDF_XED,0,1);
for nn = 1:Ndim
    id = theta_Z(:,nn)==inf|theta_Z(:,nn)==-inf;
    theta_Z(id,nn) = mean(theta_Z(~id,nn));
end

%%
ndim =  SiteGMM.Ndim/SiteGMM.NofComp;
% criteria#0
for nn = 1:SiteGMM.NofComp
    explainTerm_comp{nn} = initial_explainTerm;
end

%  criteria#1: the terms with Pearson correlation < 0.05 are removed
for nn = 1:SiteGMM.NofComp
    rho_comp{nn} = corr(theta_Z(:,(nn-1)*ndim+1:nn*ndim),XX,'type','Spearman');
    abs_rho{nn} = abs(rho_comp{nn});
    abs_rho{nn}(isnan(abs_rho{nn})) = 999;
    explainTerm_comp{nn}(abs_rho{nn}<0.05) = 0;
end


%  criteria#2: the terms with variance contribution < 0.05 are removed
explainTerm = [];
for nn = 1:SiteGMM.NofComp
    explainTerm = [explainTerm;explainTerm_comp{nn}];
end
IF_Beta_n = zeros(size(explainTerm));
id_fit = find(any(IF_Beta_n < 0.01, 2));

%% Fit a linear regression model
beta = [];
Sigma2_T = [];
std_XX = std(XX);
errY = [];

while ~isempty(id_fit)
    [beta,Sigma2_T,IF_Beta_n,errY]=lmeModel_fit(theta_Z,XX,explainTerm,Z_corr, ...
                              std_XX,id_fit,beta,Sigma2_T,IF_Beta_n,errY);

    IF_Beta_n(explainTerm==0) = 999; IF_Beta_n(:,1) = 999;
    id_fit  = find(any(IF_Beta_n < 0.01, 2));
    explainTerm(IF_Beta_n<0.01) = 0;
end

%%  criteria#3: Exclude the terms w.r.t. two horizontal components with different signs
for nn = 1:SiteGMM.NofComp
    explainTerm_comp{nn}= explainTerm((nn-1)*ndim+1:nn*ndim,:);
end
if SiteGMM.NofComp>=2
   explainTerm_H1 =   explainTerm_comp{1};
   explainTerm_H2 =   explainTerm_comp{2};
   id_differ = explainTerm_H1~=explainTerm_H2;
   explainTerm_H1(id_differ) = 0; explainTerm_H2(id_differ) = 0;
   explainTerm(1:ndim,:) = explainTerm_H1;
   explainTerm(ndim+1:2*ndim,:) = explainTerm_H2;
   id_fit = [];
   id_fit = find(any(id_differ==1, 2));
   id_fit = [id_fit;id_fit+ndim];
   [beta,Sigma2_T,IF_Beta_n,errY] = lmeModel_fit(theta_Z,XX,explainTerm,Z_corr, ...
       std_XX,id_fit,beta,Sigma2_T,IF_Beta_n,errY);
end

for nn = 1:SiteGMM.NofComp
    explainTerm_comp{nn} = explainTerm((nn-1)*ndim+1:nn*ndim,:);
end
%% Fit the resudual marginals
corrMat = corr(errY,errY,'type','Pearson');
% clear InputOpts
% InputOpts.Inference.Data = errY;
% InputOpts.Copula.Type = 'Gaussian';
% UqLabInput_GMPEerr = uq_createInput(InputOpts);


%% Save Data
SiteGMM.GMPEs.beta = beta;
SiteGMM.GMPEs.Sigma2_T = Sigma2_T;
SiteGMM.GMPEs.err_corrMat = corrMat;
SiteGMM.GMPEs.basisfuns = GMPE_BasisFuns;

% save SiteGMM_607NGArecords GM SiteGMM




