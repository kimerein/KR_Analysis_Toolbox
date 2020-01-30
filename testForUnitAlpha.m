function [correlates,lowVals,highVals]=testForUnitAlpha(low,high,t,f)

timeWindow=[3 3.6];
% freqWindow=[10 14];
freqWindow=[20 40];
% freqWindow=[0 5];

lowVals=nan(1,size(low,3));
highVals=nan(1,size(low,3));
correlates=zeros(1,size(low,3));
for i=1:size(low,3)
    lowVals(i)=reshape(nanmean(nanmean(low(t>=timeWindow(1) & t<=timeWindow(2),f>=freqWindow(1) & f<=freqWindow(2),i),1),2),1,1);
    highVals(i)=reshape(nanmean(nanmean(high(t>=timeWindow(1) & t<=timeWindow(2),f>=freqWindow(1) & f<=freqWindow(2),i),1),2),1,1);
    if lowVals(i)<highVals(i)
        correlates(i)=1;
    end
end


% all_units_noTheta_lowF1_Ntsr1=backup_all_units_noTheta_lowF1_Ntsr1;
% all_units_noTheta_highF1_Ntsr1=backup_all_units_noTheta_highF1_Ntsr1;
% all_units_theta_lowF1_Ntsr1=backup_all_units_theta_lowF1_Ntsr1;
% all_units_theta_highF1_Ntsr1=backup_all_units_theta_highF1_Ntsr1;
% all_units_noTheta_lowF1_Ntsr1=all_units_noTheta_lowF1_Ntsr1(:,:,(highVals-lowVals>0.05  & highVals_above-lowVals_above<0.05 & highVals_above-lowVals_above>-0.05));
% all_units_noTheta_highF1_Ntsr1=all_units_noTheta_highF1_Ntsr1(:,:,(highVals-lowVals>0.05  & highVals_above-lowVals_above<0.05 & highVals_above-lowVals_above>-0.05));
% test=reshape(nanmean(all_units_noTheta_lowF1_Ntsr1,3),size(all_units_noTheta_lowF1_Ntsr1,1),size(all_units_noTheta_lowF1_Ntsr1,2));
% test_high=reshape(nanmean(all_units_noTheta_highF1_Ntsr1,3),size(all_units_noTheta_highF1_Ntsr1,1),size(all_units_noTheta_highF1_Ntsr1,2));
% figure(); imagesc(noTheta_V1.allS.t(noTheta_V1.allS.t<maxTime),noTheta_V1.allS.f(noTheta_V1.allS.f<=maxFreq),test_high(noTheta_V1.allS.t<maxTime,noTheta_V1.allS.f<=maxFreq)'-test(noTheta_V1.allS.t<12,noTheta_V1.allS.f<=maxFreq)'); title('high minus low');
% figure(); imagesc(noTheta_V1.allS.t(noTheta_V1.allS.t<maxTime),noTheta_V1.allS.f(noTheta_V1.allS.f<=maxFreq),[test_high(noTheta_V1.allS.t<maxTime,noTheta_V1.allS.f<=maxFreq) test(noTheta_V1.allS.t<12,noTheta_V1.allS.f<=maxFreq)]'); title('high then low');
% figure(); imagesc(noTheta_V1.allS.t(noTheta_V1.allS.t<maxTime),noTheta_V1.allS.f(noTheta_V1.allS.f<=maxFreq),[test_high(noTheta_V1.allS.t<maxTime,noTheta_V1.allS.f<=maxFreq)]'); title('high');
% figure(); imagesc(noTheta_V1.allS.t(noTheta_V1.allS.t<maxTime),noTheta_V1.allS.f(noTheta_V1.allS.f<=maxFreq),[test(noTheta_V1.allS.t<maxTime,noTheta_V1.allS.f<=maxFreq)]'); title('low');