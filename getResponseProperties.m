function [Sled,Snoled,F1amp,HFamp,alphaRatio]=getResponseProperties(psth,takeTrials,Sled,Snoled)

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
params.trialave=1;
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
F1amp=zeros(1,length(psth.psths));
F2amp=zeros(1,length(psth.psths));
HFamp=zeros(1,length(psth.psths));
LFamp=zeros(1,length(psth.psths));
meanEvokedRate=zeros(1,length(psth.psths));
meanPrestimRate=zeros(1,length(psth.psths));
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
    F1amp(i)=sqrt(nanmean(nanmean(S(t>=stimWindow(1) & t<=stimWindow(2),f>=F1freq-0.5 & f<=F1freq+0.5),1),2));
    F2amp(i)=sqrt(nanmean(nanmean(S(t>=stimWindow(1) & t<=stimWindow(2),f>=F2freq-0.5 & f<=F2freq+0.5),1),2));
    meanEvokedRate(i)=nanmean(nanmean(p(ismember(trialscond,takeTrials) & ismember(ledcond,LEDcond) & ismember(stimcond,useStimcond),t>=stimWindow(1) & t<=stimWindow(2)),1),2);
    
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
    HFamp(i)=sqrt(nanmean(nanmean(S(t>=preStimWindow(1) & t<=preStimWindow(2),f>=HFalpha_band(1) & f<=HFalpha_band(2)),1),2));
    LFamp(i)=sqrt(nanmean(nanmean(S(t>=preStimWindow(1) & t<=preStimWindow(2),f>=LFalpha_band(1) & f<=LFalpha_band(2)),1),2));
    meanPrestimRate(i)=nanmean(nanmean(p(ismember(trialscond,takeTrials) & ismember(ledcond,noLEDcond) & ismember(stimcond,useStimcond),t>=preStimWindow(1) & t<=preStimWindow(2)),1),2);
end

figure(); 
scatter(F1amp,HFamp);
[r,p]=corrcoef(F1amp,HFamp);
disp(r);
disp(p);

figure(); 
scatter(meanEvokedRate,meanPrestimRate);
[r,p]=corrcoef(meanEvokedRate,meanPrestimRate);
disp(r);
disp(p);

alphaRatio=LFamp./HFamp;
figure(); 
scatter(F1amp,alphaRatio);
[r,p]=corrcoef(F1amp,alphaRatio);
disp(r);
disp(p);

% % Regress off effect of mean firing rate
% Y=F1amp';
% X=[ones(length(meanEvokedRate),1) meanEvokedRate'];
% B=X\Y;
% yintercept=B(1);
% slope=B(2);
% ycalc=yintercept+slope*X(:,2);
% resid=Y-ycalc;
% F1amp_withoutMean=resid;
% 
% Y=HFamp';
% X=[ones(length(meanPrestimRate),1) meanPrestimRate'];
% B=X\Y;
% yintercept=B(1);
% slope=B(2);
% ycalc=yintercept+slope*X(:,2);
% resid=Y-ycalc;
% HFamp_withoutMean=resid;
% 
% Y=LFamp';
% X=[ones(length(meanPrestimRate),1) meanPrestimRate'];
% B=X\Y;
% yintercept=B(1);
% slope=B(2);
% ycalc=yintercept+slope*X(:,2);
% resid=Y-ycalc;
% LFamp_withoutMean=resid;
% 
% figure(); 
% scatter(F1amp_withoutMean,HFamp_withoutMean);
% [r,p]=corrcoef(F1amp_withoutMean,HFamp_withoutMean);
% disp(r);
% disp(p);
% 
% figure(); 
% scatter(F1amp,HFamp_withoutMean);
% [r,p]=corrcoef(F1amp,HFamp_withoutMean);
% disp(r);
% disp(p);
