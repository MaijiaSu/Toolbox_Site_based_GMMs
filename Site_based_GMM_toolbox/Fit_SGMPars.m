
function SiteGMM = Fit_SGMPars(FitOption,GM)

tstart = tic;

NofComp = size(GM,1);
NofGM = size(GM,2);
% NofGM = 10;
for n_Comp = 1:NofComp

for ir = 1:NofGM
    tend = toc(tstart);
    disp(['n_comp-',num2str(n_Comp),':',num2str(ir),'/',num2str(NofGM), ...
        ':','Elapsed time = ', num2str(sum(tend)/60,'%.1f'),'mins']);

    % fitted ground motion parameters
    SGMPar_ir = SGMParFitting(GM{n_Comp,ir},FitOption);

    % save the resutls
    SGMPar{n_Comp,ir} = SGMPar_ir;

    % save the data
    if rem(ir,100)==0
        save('SGMPars','SGMPar','FitOption');
    end

end
end

% Covert the data into a table
theta_GM = GetSGMPars(SGMPar);
% save data into a structure variables
SiteGMM.NofGM = NofGM;
SiteGMM.NofComp = NofComp;
SiteGMM.Ndim = size(theta_GM,2);
SiteGMM.FitOption = FitOption;
SiteGMM.theta_GM = theta_GM;
% SiteGMM.MetaData = MetaData;

