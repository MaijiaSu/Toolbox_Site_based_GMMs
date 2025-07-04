function SiteGMM = Fit_GMMPars_JointPDF(SiteGMM)

% Set the Supports of each paramters
DisBounds = GetSGMPars_support(SiteGMM);
tempdata = SiteGMM.theta_GM;
for nn = 1:size(DisBounds,2)
    id1 = find(tempdata(:,nn)<DisBounds(1,nn));
    id2 = find(tempdata(:,nn)>DisBounds(2,nn));
    tempdata(id1,nn)=DisBounds(1,nn);
    tempdata(id2,nn)=DisBounds(2,nn);
end
SiteGMM.theta_GM_bounded = tempdata;


% Fit the joint PDF using UQLab
InputOpts.Inference.Data = tempdata;
% for nd = 1:3
%     warning('need to adjust the marginal types')
%     InputOpts.Marginals(11*(nd-1)+1).Type = ['Gaussian'];
%     InputOpts.Marginals(11*(nd-1)+2).Type = 'LogNormal';
%     InputOpts.Marginals(11*(nd-1)+3).Type = 'LogNormal';
%     InputOpts.Marginals(11*(nd-1)+4).Type = 'LogNormal';
%     InputOpts.Marginals(11*(nd-1)+5).Type = 'LogNormal';
%     InputOpts.Marginals(11*(nd-1)+6).Type = 'LogNormal';
%     InputOpts.Marginals(11*(nd-1)+7).Type = {'LogNormal'};
%     InputOpts.Marginals(11*(nd-1)+8).Type = {'LogNormal'};
%     InputOpts.Marginals(11*(nd-1)+9).Type = 'Laplace';
%     InputOpts.Marginals(11*(nd-1)+10).Type = [];
%     InputOpts.Marginals(11*(nd-1)+11).Type = [];
% end
InputOpts.Copula.Type = 'independent';
UqLabInput = uq_createInput(InputOpts);




% Save Data
for nn = 1:size(DisBounds,2)
    nPar = numel(UqLabInput.Marginals(nn).Parameters);
    MarginalPar(nn,1:nPar) = UqLabInput.Marginals(nn).Parameters;
    nPar = numel(UqLabInput.Marginals(nn).Moments);
    MarginalMoment(nn,1:nPar) = UqLabInput.Marginals(nn).Moments;
end
%
JointPDF.Marginals.Types = {UqLabInput.Marginals.Type};
JointPDF.Copula = UqLabInput.Copula;
JointPDF.Marginals.Par = MarginalPar;
JointPDF.Marginals.Moments = MarginalMoment;
JointPDF.UqLabInput = UqLabInput;
JointPDF.DisBounds = DisBounds;
SiteGMM.Pars_JointPDF = JointPDF;

return