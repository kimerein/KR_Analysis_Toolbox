function getPeakResponseFromPSTH_trialAv(psth,TF,stimWindow,spontWindow,outDir,noThetaTrials,uses,uses_tri)

ledOn=5.05;
ledOff=0;

t=psth.t;
l=psth.unitLED{1};
TFbins=stimWindow(1):1/TF:stimWindow(2);
TFbins_spont=spontWindow(1):1/TF:spontWindow(2);

noLED_ev_noTheta=nan(1,length(psth.psths));
LED_ev_noTheta=nan(1,length(psth.psths));
noLED_ev_theta=nan(1,length(psth.psths));
LED_ev_theta=nan(1,length(psth.psths));

peakPSTHs_noLED_noTheta=cell(1,length(psth.psths));
peakPSTHs_LED_noTheta=cell(1,length(psth.psths));
peakPSTHs_noLED_theta=cell(1,length(psth.psths));
peakPSTHs_LED_theta=cell(1,length(psth.psths));

peakPSTHs_noLED_noTheta_av=nan(length(psth.psths),(length(TFbins_spont)-1)+(length(TFbins)-1));
peakPSTHs_LED_noTheta_av=nan(length(psth.psths),(length(TFbins_spont)-1)+(length(TFbins)-1));
peakPSTHs_noLED_theta_av=nan(length(psth.psths),(length(TFbins_spont)-1)+(length(TFbins)-1));
peakPSTHs_LED_theta_av=nan(length(psth.psths),(length(TFbins_spont)-1)+(length(TFbins)-1));

s=psth.unitStimcond{1};
tri=psth.unitTrials{1};

for i=1:length(psth.psths)
    p=psth.psths{i};
    
    [stimPeaks,spontPeaks,psthOfPeaks]=getPeaksPSTH(p,t,uses,uses_tri,TFbins_spont,TFbins,l,ledOff,s,tri,noThetaTrials,1);
    noLED_ev_noTheta(i)=nanmean(stimPeaks-spontPeaks);
    peakPSTHs_noLED_noTheta{i}=psthOfPeaks;
    peakPSTHs_noLED_noTheta_av(i,:)=nanmean(psthOfPeaks,1);
    
    [stimPeaks,spontPeaks,psthOfPeaks]=getPeaksPSTH(p,t,uses,uses_tri,TFbins_spont,TFbins,l,ledOn,s,tri,noThetaTrials,1);
    LED_ev_noTheta(i)=nanmean(stimPeaks-spontPeaks);
    peakPSTHs_LED_noTheta{i}=psthOfPeaks;
    peakPSTHs_LED_noTheta_av(i,:)=nanmean(psthOfPeaks,1);
    
    [stimPeaks,spontPeaks,psthOfPeaks]=getPeaksPSTH(p,t,uses,uses_tri,TFbins_spont,TFbins,l,ledOff,s,tri,noThetaTrials,0);
    noLED_ev_theta(i)=nanmean(stimPeaks-spontPeaks);
    peakPSTHs_noLED_theta{i}=psthOfPeaks;
    peakPSTHs_noLED_theta_av(i,:)=nanmean(psthOfPeaks,1);
    
    [stimPeaks,spontPeaks,psthOfPeaks]=getPeaksPSTH(p,t,uses,uses_tri,TFbins_spont,TFbins,l,ledOn,s,tri,noThetaTrials,0);
    LED_ev_theta(i)=nanmean(stimPeaks-spontPeaks);
    peakPSTHs_LED_theta{i}=psthOfPeaks;
    peakPSTHs_LED_theta_av(i,:)=nanmean(psthOfPeaks,1); 
end

save([outDir '\' 'psth.mat'],'psth');
save([outDir '\' 'noLED_ev_noTheta.mat'],'noLED_ev_noTheta');   
save([outDir '\' 'LED_ev_noTheta.mat'], 'LED_ev_noTheta');
save([outDir '\' 'noLED_ev_theta.mat'],'noLED_ev_theta');
save([outDir '\' 'LED_ev_theta.mat'],'LED_ev_theta');

save([outDir '\' 'peakPSTHs_noLED_noTheta.mat'],'peakPSTHs_noLED_noTheta');
save([outDir '\' 'peakPSTHs_LED_noTheta.mat'],'peakPSTHs_LED_noTheta');
save([outDir '\' 'peakPSTHs_noLED_theta.mat'],'peakPSTHs_noLED_theta');
save([outDir '\' 'peakPSTHs_LED_theta.mat'],'peakPSTHs_LED_theta');

save([outDir '\' 'peakPSTHs_noLED_noTheta_av.mat'],'peakPSTHs_noLED_noTheta_av');
save([outDir '\' 'peakPSTHs_LED_noTheta_av.mat'],'peakPSTHs_LED_noTheta_av');
save([outDir '\' 'peakPSTHs_noLED_theta_av.mat'],'peakPSTHs_noLED_theta_av');
save([outDir '\' 'peakPSTHs_LED_theta_av.mat'],'peakPSTHs_LED_theta_av');
end

function [stimPeaks,spontPeaks,psthOfPeaks]=getPeaksPSTH(p,backup_t,uses,uses_tri,TFbins_spont,TFbins,l,ledVal,s,tri,noThetaTrials,noThetaVal)

stimPeaks=nan(1,length(uses));
spontPeaks=nan(1,length(uses));
psthOfPeaks=nan(length(uses),(length(TFbins_spont)-1)+(length(TFbins)-1));
t=downSampAv(backup_t,floor((0.3333/4)/(backup_t(2)-backup_t(1))));
for j=1:length(uses)
    currs=uses{j};
    currtri=uses_tri{j};
    takeTrials=zeros(size(p,1),1);
    for i=1:length(currs)
        subcurrs=currs{i};
        subcurrtri=currtri{i};
        takeTrials(ismember(s,subcurrs) & ismember(tri,subcurrtri))=1;
    end
    currp=nanmean(p(ismember(l',ledVal) & takeTrials==1 & noThetaTrials==noThetaVal,:),1);
    currp=downSampAv(currp,floor((0.3333/4)/(backup_t(2)-backup_t(1))));
    trialStimPeaks=nan(1,length(TFbins)-1);
    trialSpontPeaks=nan(1,length(TFbins_spont)-1);
    for k=1:length(TFbins)-1
        trialStimPeaks(k)=max(currp(t>=TFbins(k) & t<=TFbins(k+1)));
    end
    for k=1:length(TFbins_spont)-1
        trialSpontPeaks(k)=max(currp(t>=TFbins_spont(k) & t<=TFbins_spont(k+1)));
    end
    psthOfPeaks(j,:)=[trialSpontPeaks trialStimPeaks];
    stimPeaks(j)=nanmean(trialStimPeaks);
    spontPeaks(j)=nanmean(trialSpontPeaks);
end
    
end