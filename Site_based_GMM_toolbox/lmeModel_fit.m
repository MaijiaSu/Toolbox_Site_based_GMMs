
function [beta,Sigma2_T,IF_Beta_n,errY]=lmeModel_fit(theta_Z,XX,explainTerm,Z_corr, ...
                              std_XX,id_fit,beta,Sigma2_T,IF_Beta_n,errY)
% 1. Initialization

% 
for nn = id_fit'
        % ouput observations
        Y = theta_Z(:,nn);
        idX = logical(explainTerm(nn,:));
        X = XX(:,idX);
        if sum(idX) == 0
            X = XX(:,1:2);
            idX = logical([1,1,0,0,0,0,0]);
        end
        lmeStr = {'Intercept'};
        LMformula = 'Y ~ 1';
        for kk = 1:size(X,2)-1
            lmeStr(kk+1) = {['X',num2str(kk)]};
            LMformula = [LMformula,'+',lmeStr(kk+1)];
        end
        LMformula = strjoin(LMformula);
        % Fit liner random-effect model
        lmeModel = fitlmematrix(X,Y,Z_corr,[], ...
            'FixedEffectPredictors',lmeStr,...
            'CovariancePattern','Isotropic');

        % Fit a linear model (without random effect)
        tbl = array2table([Y,X(:,2:end)], 'VariableNames',['Y',lmeStr(2:end)]);
        reducelme = fitlme(tbl,LMformula);
        Ftest_results = compare(reducelme, lmeModel);
        tau_pvalue(nn,1) = Ftest_results.pValue(2);

        % Results
        beta_LM = fixedEffects(lmeModel);
        b_LM = randomEffects(lmeModel);
        [psi,mse,~] = covarianceParameters(lmeModel);

        % Save results
        % Coefficeints
        beta(nn,idX) = beta_LM;
        % Inter variance
        sigma2_inter(nn,1) = mse;
        % Intra variance
        tau2_intra(nn,1) = psi{1,1}(1,1);
        % Total varince
        Sigma2_T(nn,1) = sqrt(sigma2_inter(nn,1)+ tau2_intra(nn,1))^2;
        % R2
        Rsquare(nn,1) = lmeModel.Rsquared.Ordinary;
        Rsquare_adjusted(nn,1) = lmeModel.Rsquared.Adjusted;
        % Mean of prediction
        muY_mixfit{nn} = X*beta_LM+Z_corr*b_LM;
        muY_fix{nn} = X*beta_LM;
        % Intra and inter residuals
        intrer_error(nn,:) = Z_corr*b_LM;
        intra_error(nn,:) = Y- muY_mixfit{nn};
        % intra + inter biases
        errY(:,nn) = Y- muY_fix{nn};
        Sigma2_errY(nn,1) = var(errY(:,nn));
        pvalue(nn,idX) = lmeModel.Coefficients.pValue;
        % Varaince contribution
        TotalV = sum((std_XX.*beta(nn,:)).^2);
        IF_Beta_n(nn,:) = (std_XX.*beta(nn,:)).^2/TotalV;
end