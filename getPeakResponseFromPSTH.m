function getPeakResponseFromPSTH(psth,TF,stimWindow,spontWindow,outDir,noThetaTrials)

ledOn=5.05;
ledOff=0;

t=psth.t;
l=psth.unitLED{1};
TFbins=stimWindow(1):1/TF:stimWindow(2);
TFbins_spont=spontWindow(1):1/TF:spontWindow(2);
noLED_stimPeak=nan(1,length(psth.psths));
LED_stimPeak=nan(1,length(psth.psths));
noLED_spontPeak=nan(1,length(psth.psths));
LED_spontPeak=nan(1,length(psth.psths));

noLED_stimPeak_trials=nan(length(psth.psths),sum(ismember(l,ledOff)));
LED_stimPeak_trials=nan(length(psth.psths),sum(ismember(l,ledOn)));
noLED_spontPeak_trials=nan(length(psth.psths),sum(ismember(l,ledOff)));
LED_spontPeak_trials=nan(length(psth.psths),sum(ismember(l,ledOn)));

peakPSTHs=cell(1,length(psth.psths));

for i=1:length(psth.psths)
    p=psth.psths{i};
    stimPeaks=nan(1,size(p,1));
    spontPeaks=nan(1,size(p,1));
    psthOfPeaks=nan(size(p,1),(length(TFbins_spont)-1)+(length(TFbins)-1));
    for j=1:size(p,1)
        currp=p(j,:);
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
    peakPSTHs{i}=psthOfPeaks;
    noLED_stimPeak(i)=nanmean(stimPeaks(ismember(l,ledOff)));
    LED_stimPeak(i)=nanmean(stimPeaks(ismember(l,ledOn)));
    noLED_spontPeak(i)=nanmean(spontPeaks(ismember(l,ledOff)));
    LED_spontPeak(i)=nanmean(spontPeaks(ismember(l,ledOn)));
    
    noLED_stimPeak_trials(i,:)=stimPeaks(ismember(l,ledOff));
    LED_stimPeak_trials(i,:)=stimPeaks(ismember(l,ledOn));
    noLED_spontPeak_trials(i,:)=spontPeaks(ismember(l,ledOff));
    LED_spontPeak_trials(i,:)=spontPeaks(ismember(l,ledOn));
end

save([outDir '\' 'psth.mat'],'psth');
save([outDir '\' 'peakPSTHs.mat'],'peakPSTHs');   
save([outDir '\' 'noLED_stimPeak.mat'], 'noLED_stimPeak');
save([outDir '\' 'LED_stimPeak.mat'],'LED_stimPeak');
save([outDir '\' 'noLED_spontPeak.mat'],'noLED_spontPeak');
save([outDir '\' 'LED_spontPeak.mat'],'LED_spontPeak');

save([outDir '\' 'noLED_stimPeak_trials.mat'],'noLED_stimPeak_trials');
save([outDir '\' 'LED_stimPeak_trials.mat'],'LED_stimPeak_trials');
save([outDir '\' 'noLED_spontPeak_trials.mat'],'noLED_spontPeak_trials');
save([outDir '\' 'LED_spontPeak_trials.mat'],'LED_spontPeak_trials');
save([outDir '\' 'l.mat'],'l');
