function FittedSGMPars = GetSGMPars(SGMPar_aug)
GMPar_T = [];
GMPar_F = [];
Wf = [];
NofGM = size(SGMPar_aug,2);
NofComp = size(SGMPar_aug,1);

FittedSGMPars = [];

for n_comp = 1:NofComp
    for ir = 1:NofGM
        SGMPar_ir = SGMPar_aug{n_comp,ir};
        % Time-modulating function
        if strcmp(SGMPar_ir.Tmodulate.type,'spline')
            nPar_T = 7;
            t_AI_p(ir,:) = SGMPar_ir.Tmodulate.Par.t_AI_p(1:end);
            Dt_AIp(ir,:) = t_AI_p(ir,2:end)-t_AI_p(ir,1:end-1);
            AI(ir,1) = SGMPar_ir.Tmodulate.Par.AI(end);
            AI_log(ir,1) = log(AI(ir,1));
            GMPar_T(ir,:) = [AI_log(ir,1), Dt_AIp(ir,:)];
            Parstr_T = {'$\log(I_a)$','$D_{0-5}$','$D_{5-30}$','$D_{30-45}$',...
                '$D_{45-75}$','$D_{75-95}$','$D_{95-100}$'};

        elseif strcmp(SGMPar_ir.Tmodulate.type,'Gamma')
            nPar_T = 3;
            AI(ir,1) = SGMPar_ir.Tmodulate.Par.AI(end);
            AI_log(ir,1) = log(AI(ir,1));
            D5_95(ir,1)  = SGMPar_ir.Tmodulate.Par.D5_95;
            tmid(ir,1)  = SGMPar_ir.Tmodulate.Par.tmid;
            GMPar_T(ir,:) = [AI_log(ir,1),D5_95(ir,:),tmid(ir,1)];
            Parstr_T = {'$\log(I_a)$','$t_{mid}$','$D_{5-95}$'};
        end

        % Frequency-modulating function
        trendmodel = SGMPar_ir.Fmodulate.Paramtericmodel;
        if ischar(trendmodel)
            if strcmp(trendmodel,'ReducedLinearmodel')
                trendmodel = [2,1];
            end
        end
        [GMFilter,FilterNofPar] = SetFilterModel(SGMPar_aug{ir}.Fmodulate.type);
        nofTrendPars = zeros(1,FilterNofPar);
        nofTrendPars(find(trendmodel==1)) = 1;
        nofTrendPars(find(trendmodel==2)) = 2;
        nofTrendPars(find(trendmodel==3)) = 5;
        tempPar = [];
        for n_trend = 1:FilterNofPar-1
            tempPar =  [tempPar;SGMPar_ir.Fmodulate.opt_par(1:nofTrendPars(n_trend),n_trend)];
        end
        GMPar_F(ir,:) = tempPar';

        % High-pass filter
        Wf(ir,1) = SGMPar_ir.HighPassFilter.Wf/2/pi;

    end

    % Numerical issue
    % Wf(find(Wf<=0)) = 0.01;
    temp = [GMPar_T,GMPar_F,Wf];
    FittedSGMPars = [FittedSGMPars,temp];
    % ParStr = {'$\log(I_a)$','$\omega(t_{mid})$','$\omega^{\prime}(t_{mid})$',...
    %     '$\zeta(t_{mid})$','$D_{0-5}$','$D_{5-30}$','$D_{30-45}$',...
    %     '$D_{45-75}$','$D_{75-95}$','$D_{95-100}$','$f_{c}$'};
end
return