% CUT_off_FREQ = 25; 
% dteq = dt_H1;
% eq = acc_major;

function [eq_d,dt_d] = MyGMdownsampling(CUT_off_FREQ,dteq,eq)

if isempty(dteq)
    eq_d=[];
    dt_d=[];
else
    dt_target = 1/2/CUT_off_FREQ;       % downsampled rate
    d_rate = dt_target/dteq;          % downsampling rate
    
    if d_rate < 1  % return the input singnal
         eq_d = eq;
         dt_d = dteq;
    else
        if floor(d_rate) == d_rate
            eq_d = decimate(eq,d_rate);              % downsampled time series
            dt_d = dt_target;
        else % the dowsampling rate is not a integer
            d_rate = floor(dt_target/dteq);          % downsampling rate
            eq_d = decimate(eq,d_rate);
            dt_d = dteq*d_rate;
        end
    end
end
end

%%
% clc,clear,close all
% MyDefaultSet
% %% Show the sampling interval
% load GM_H1
% dtH1 = [GM(:).dt];
% load GM_H2
% dtH2 = [GM(:).dt];
% load GM_V
% tmp_dtV = [GM(:).dt];
% load EQPar
% NofGM = numel(dtH2);
% dtV = zeros(1,NofGM);
% V_flag = [EQPar(:).V_flag];
% dtV(V_flag) = tmp_dtV;
% 
% figure('Position',[100,100,800,400])
% plot(1:NofGM,dtH1,'o')
% hold on
% plot(1:NofGM,dtH2,'x')
% hold on
% plot(1:NofGM,dtV,'s')
% grid on
% xlabel('No.')
% ylabel('$dt$')
% xlim([1,NofGM])
% legend('H1','H2','V','box','off')
% saveas(gcf,'fig1.png')
% 
% %% number of records with different dt
% dtr = sort(unique(dtH1),'descend')';
% 
% for ii = 1:numel(dtr)
%    N(ii) = numel(find(dtH1==dtr(ii)));
% end
% 
% tab = [dtr';N]
% % N = histcounts(dtH1,[0;dtr]')
% 
% 
% %% Show the number of data points
% load GM_H1
% LH1 = [GM(:).L];
% load GM_H2
% LH2 = [GM(:).L];
% load GM_V
% tmp_LV = [GM(:).L];
% load EQPar
% NofGM = numel(LH2);
% LV = zeros(1,NofGM);
% V_flag = [EQPar(:).V_flag];
% LV(V_flag) = tmp_LV;
% 
% figure('Position',[100,100,800,400])
% plot(1:NofGM,LH1,'o')
% hold on
% plot(1:NofGM,LH2,'x')
% hold on
% plot(1:NofGM,LV,'s')
% grid on
% xlabel('No.')
% ylabel('$L$')
% xlim([1,NofGM])
% legend('H1','H2','V','box','off')
% saveas(gcf,'fig2.png')
% 
% 
