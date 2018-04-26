function [pvals,likes_psth1]=getSigDiffs(psth1,psth2)

uses1=[1];
uses2=[1];
stimWindow=[4 6.5];
baseWindow=[0 4];
usel=[5.05];
F1freq=3;

params.tapers=[5 6];
params.Fs=1/(psth1.t(2)-psth1.t(1));
params.fpass=[F1freq-0.5 F1freq+0.5];
params.trialave=0;

pvals=nan(1,length(psth1.psths));
likes_psth1=nan(1,length(psth1.psths));
for i=1:length(psth1.psths)
    p1=psth1.psths{i};
    l1=psth1.unitLED{i};
    s1=psth1.unitStimcond{i};
    p2=psth2.psths{i};
    l2=psth2.unitLED{i};
    s2=psth2.unitStimcond{i};
    [F1amp_psth1,f1_F1]=mtspectrumpb(p1(ismember(l1,usel) & ismember(s1,uses1),psth1.t>=stimWindow(1) & psth1.t<=stimWindow(2))',params);
    [baseamp_psth1,f1_base]=mtspectrumpb(p1(ismember(l1,usel) & ismember(s1,uses1),psth1.t>=baseWindow(1) & psth1.t<=baseWindow(2))',params);
    [F1amp_psth2,f2_F1]=mtspectrumpb(p2(ismember(l2,usel) & ismember(s2,uses2),psth2.t>=stimWindow(1) & psth2.t<=stimWindow(2))',params);
    [baseamp_psth2,f2_base]=mtspectrumpb(p2(ismember(l2,usel) & ismember(s2,uses2),psth2.t>=baseWindow(1) & psth2.t<=baseWindow(2))',params);
    trialF1amps_psth1=nanmean(F1amp_psth1,1)-nanmean(baseamp_psth1,1);
    trialF1amps_psth2=nanmean(F1amp_psth2,1)-nanmean(baseamp_psth2,1);
    pvals(i)=ranksum(trialF1amps_psth1,trialF1amps_psth2);
    likes_psth1(i)=nanmean(trialF1amps_psth1)>nanmean(trialF1amps_psth2);
%     pvals(i)=ranksum(nanmean(F1amp_psth1,1),nanmean(baseamp_psth1,1));
%     likes_psth1(i)=nanmean(nanmean(F1amp_psth1,1))>nanmean(nanmean(baseamp_psth1,1));
end

    
    