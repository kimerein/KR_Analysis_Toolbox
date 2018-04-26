function [lowF1Trials,highF1Trials]=findTrialsWithHighOrLowF1(noTheta_trialAv_noLED,psth,uses,uses_tri,noThetaTrials,isNoTheta)

useLED=0;
stimWindow=[4 6.5];
F1range=[2.5 3.5];

isHigh=nan(size(psth.psths{1},1),1);
s=psth.unitStimcond{1};
tri=psth.unitTrials{1};
l=psth.unitLED{1};

currt=noTheta_trialAv_noLED(1).allS.t;
currf=noTheta_trialAv_noLED(1).allS.f;
for j=1:length(noTheta_trialAv_noLED)
    snow=uses{j};
    trinow=uses_tri{j};
    isInCurrTrialSet=zeros(1,length(tri));
    for i=1:length(snow)
        isInCurrTrialSet(ismember(tri,trinow{i}) & ismember(s,snow{i}) & ismember(l,useLED) & noThetaTrials'==isNoTheta)=1;
    end
        
    F1s=nan(length(noTheta_trialAv_noLED(j).allS.S),size(noTheta_trialAv_noLED(j).allS.S{1},3));
    allFs=nan(length(noTheta_trialAv_noLED(j).allS.S),size(noTheta_trialAv_noLED(j).allS.S{1},3));
    for k=1:length(noTheta_trialAv_noLED(j).allS.S)
        currs=noTheta_trialAv_noLED(j).allS.S{k};
        F1s(k,:)=reshape(nanmean(nanmean(currs(currt>=stimWindow(1) & currt<=stimWindow(2),currf>=F1range(1) & currf<=F1range(2),:),2),1),1,size(currs,3));
        allFs(k,:)=reshape(nanmean(nanmean(currs(currt>=stimWindow(1) & currt<=stimWindow(2),currf>=0 & currf<=50,:),2),1),1,size(currs,3));
    end
    F1s=nanmean(F1s,1)./nanmean(allFs,1);
    medF1=median(F1s(~isnan(F1s)));
    
    isHigh(isInCurrTrialSet==1)=F1s>medF1;
end
highF1Trials=tri(isHigh==1);
lowF1Trials=tri(isHigh==0);