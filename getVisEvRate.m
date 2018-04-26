function [visevrates_noLED_noTheta,visevrates_LED_noTheta,visevrates_noLED_theta,visevrates_LED_theta]=getVisEvRate(psth,noLED,LED,noThetaTrials)

stimWindow=[4 6.5];

t=psth.t;
l=psth.unitLED{1};

visevrates_noLED_noTheta=nan(1,length(psth.psths));
visevrates_LED_noTheta=nan(1,length(psth.psths));
visevrates_noLED_theta=nan(1,length(psth.psths));
visevrates_LED_theta=nan(1,length(psth.psths));
for i=1:length(psth.psths)
    p=psth.psths{i};
    
    taketri=noThetaTrials' & ismember(l,noLED);
    subp=p(taketri,:);
    visevrates_noLED_noTheta(i)=nanmean(nanmean(subp(:,t>=stimWindow(1) & t<=stimWindow(2)),2),1);
    
    taketri=noThetaTrials' & ismember(l,LED);
    subp=p(taketri,:);
    visevrates_LED_noTheta(i)=nanmean(nanmean(subp(:,t>=stimWindow(1) & t<=stimWindow(2)),2),1);
    
    taketri=~noThetaTrials' & ismember(l,noLED);
    subp=p(taketri,:);
    visevrates_noLED_theta(i)=nanmean(nanmean(subp(:,t>=stimWindow(1) & t<=stimWindow(2)),2),1);
    
    taketri=~noThetaTrials' & ismember(l,LED);
    subp=p(taketri,:);
    visevrates_LED_theta(i)=nanmean(nanmean(subp(:,t>=stimWindow(1) & t<=stimWindow(2)),2),1);
end
    