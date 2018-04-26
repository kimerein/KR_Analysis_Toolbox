function [Sled,Snoled,F1amp,HFamp,alphaRatio]=getResponseProperties_trialByTrial(psth,takeTrials,Sled,Snoled)

noLEDcond=0;
LEDcond=5.05;
stimWindow=[4 6.5];
preStimWindow=[3.2 4];
F1freq=3;
F2freq=6;
HFalpha_band=[11.5 12.5];
LFalpha_band=[4 6];

useStimcond=[1:10000];

params.tapers=[3 10];
params.Fs=1/(psth.t(2)-psth.t(1));
params.fpass=[1 80];
params.trialave=0;
movingwin=[0.8 0.05];

trialscond=psth.unitTrials{1};
stimcond=psth.unitStimcond{1};
ledcond=psth.unitLED{1};

if isempty(Sled)
    Sled.S=cell(1,length(psth.psths));
    Sled.t=cell(1,length(psth.psths));
    Sled.f=cell(1,length(psth.psths));
    redo=1;
else
    redo=0;
end
if isempty(Snoled)
    Snoled.S=cell(1,length(psth.psths));
    Snoled.t=cell(1,length(psth.psths));
    Snoled.f=cell(1,length(psth.psths));
    redo=1;
else
    redo=0;
end
F1amp=cell(1,length(psth.psths));
F2amp=cell(1,length(psth.psths));
HFamp=cell(1,length(psth.psths));
LFamp=cell(1,length(psth.psths));
meanEvokedRate=cell(1,length(psth.psths));
meanPrestimRate=cell(1,length(psth.psths));
for i=1:length(psth.psths)
    p=psth.psths{i};
    disp(i);
    if redo==1
        [S,t,f]=mtspecgrampb(p(ismember(trialscond,takeTrials) & ismember(ledcond,LEDcond) & ismember(stimcond,useStimcond),:)',movingwin,params);
        Sled.S{i}=S;
        Sled.t{i}=t;
        Sled.f{i}=f;
    else
        S=Sled.S{i};
        t=Sled.t{i};
        f=Sled.f{i};
    end
    F1amp{i}=sqrt(nanmean(nanmean(S(t>=stimWindow(1) & t<=stimWindow(2),f>=F1freq-0.5 & f<=F1freq+0.5,:),1),2));
    F2amp{i}=sqrt(nanmean(nanmean(S(t>=stimWindow(1) & t<=stimWindow(2),f>=F2freq-0.5 & f<=F2freq+0.5,:),1),2));
    meanEvokedRate{i}=nanmean(p(ismember(trialscond,takeTrials) & ismember(ledcond,LEDcond) & ismember(stimcond,useStimcond),t>=stimWindow(1) & t<=stimWindow(2)),2);
    
    if redo==1
        [S,t,f]=mtspecgrampb(p(ismember(trialscond,takeTrials) & ismember(ledcond,noLEDcond) & ismember(stimcond,useStimcond),:)',movingwin,params);
        Snoled.S{i}=S;
        Snoled.t{i}=t;
        Snoled.f{i}=f;
    else
        S=Snoled.S{i};
        t=Snoled.t{i};
        f=Snoled.f{i};
    end
    HFamp{i}=sqrt(nanmean(nanmean(S(t>=preStimWindow(1) & t<=preStimWindow(2),f>=HFalpha_band(1) & f<=HFalpha_band(2),:),1),2));
    LFamp{i}=sqrt(nanmean(nanmean(S(t>=preStimWindow(1) & t<=preStimWindow(2),f>=LFalpha_band(1) & f<=LFalpha_band(2),:),1),2));
    meanPrestimRate{i}=nanmean(p(ismember(trialscond,takeTrials) & ismember(ledcond,noLEDcond) & ismember(stimcond,useStimcond),t>=preStimWindow(1) & t<=preStimWindow(2)),2);
end

alphaRatio=cell(1,length(LFamp));
for i=1:length(LFamp)
    alphaRatio{i}=LFamp{i}./HFamp{i};
end