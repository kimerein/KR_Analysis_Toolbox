function [topAv,bottomAv,pref,nonpref,ampOrder,F1phases,highestAmps,takeTop]=plotPrefVersusNonpref(psth,useUnit,refFreqBand,useWindow,useWindowForAmp,ampOrdering,F1phases_input,usel)

% takeN=3;
% stimconds=1:12;
% usel=[0.05];
% takeN=5;
% stimconds=1:12;
% takeN=12;
% stimconds=1:16;
% usel=[0];
takeN=1;
stimconds=1;
showFigs=0;

psths=psth.psths{useUnit};
trialLED=psth.unitLED{useUnit};
trialStimcond=psth.unitStimcond{useUnit};
trials=psth.unitTrials{useUnit};

% OLD
% params.Fs=1/(psth.t(2)-psth.t(1));
% params.tapers=[2 5];
% params.trialave=1;
% % params.fpass=[2.5 30];
% params.fpass=[1 30];

params.Fs=1/(psth.t(2)-psth.t(1));
params.tapers=[3 10];
params.trialave=1;
% params.fpass=[2.5 30];
params.fpass=[1 30];

% from freqsynch
% trialDuration=14.5;
% stimWindow=[4 5.5];
% ds=1;
% times=linspace(0,trialDuration,size(synchPSTH,2));
% 
% params.tapers=[3 10];
% params.Fs=1./(trialDuration/size(synchPSTH,2));
% params.fpass=[1 20];
% params.trialave=0;


F1amps=nan(length(stimconds),1);
F1phases=nan(length(stimconds),1);
F1coherence=nan(length(stimconds),1);
subt=psth.t(psth.t>=useWindowForAmp(1) & psth.t<=useWindowForAmp(2));
refx=subt;
refy=sin(2*pi*nanmean(refFreqBand).*refx);
for i=1:length(stimconds)
    currs=stimconds(i);
    useTrials=ismember(trialLED,usel) & ismember(trialStimcond,currs);
    temp=psths(useTrials,psth.t>=useWindowForAmp(1) & psth.t<=useWindowForAmp(2));
    % Get F1 amp
    [S,f]=mtspectrumpb(temp',params);
%     [S,f]=mtspectrumpb((temp-repmat(nanmean(temp,2),1,size(temp,2)))',params);
%     F1amps(i)=nanmean(S(f>=refFreqBand(1) & f<=refFreqBand(2)))-nanmean([nanmean(S(f<refFreqBand(1))) nanmean(S(f>refFreqBand(2)))]);
%     F1amps(i)=nanmean(S(f>=refFreqBand(1) & f<=refFreqBand(2)))./nanmean(S);
%     DCchange=nanmean(nanmean(temp,1),2);
%     F1amps(i)=(sqrt(nanmean(S(f>=refFreqBand(1) & f<=refFreqBand(2))))/25.05)-DCchange;
    F1amps(i)=nanmean(S(f>=refFreqBand(1) & f<=refFreqBand(2)));
    % Get F1 phase
    [C,phi,~,~,~,f]=coherencycpb(repmat(refy',1,size(temp,1)),temp',params);
    F1coherence(i)=nanmean(C(f>=refFreqBand(1) & f<=refFreqBand(2)));
    F1phases(i)=nanmean(phi(f>=refFreqBand(1) & f<=refFreqBand(2)));
end

[~,ampOrder]=sort(F1amps);
if ~isempty(ampOrdering)
    ampOrder=ampOrdering;
end
if ~isempty(F1phases_input)
    F1phases=F1phases_input;
end
% takeTop=find(F1phases>=-pi/2 & F1phases<0);
takeTop=ampOrder(end-takeN+1:end);
% highestAmps=nanmean(F1amps(ampOrder(end-takeN+1:end)));
highestAmps=nanmean(F1amps(ampOrder(1:takeN)));
% takeBottom=ampOrder(1:takeN-1);
takeBottom=ampOrder(1:takeN);
topAv=nan(takeN,2*length(psth.t));
bottomAv=nan(takeN,2*length(psth.t));
refPeriod=1/nanmean(refFreqBand);
spontTop=nan(takeN,length(psth.t(psth.t<=useWindow(1))));
spontBottom=nan(takeN,length(psth.t(psth.t<=useWindow(1))));
subt=psth.t(psth.t>=useWindow(1) & psth.t<=useWindow(2));
for i=1:length(takeTop)
    currs=takeTop(i);
    useTrials=ismember(trialLED,usel) & ismember(trialStimcond,currs);
    temp=psths(useTrials,psth.t>=useWindow(1) & psth.t<=useWindow(2));
    spontTop(i,:)=nanmean(psths(useTrials,psth.t<=useWindow(1)),1);
    shiftInTime=((F1phases(currs)+pi)/(2*pi))*refPeriod;
    startInd=length(psth.t)-ceil(shiftInTime/(psth.t(2)-psth.t(1)));
    endInd=startInd+length(subt)-1;
    if isnan(F1phases(currs))
        topAv(i,:)=nan;
    else
        topAv(i,startInd:endInd)=nanmean(temp,1);
    end
end
for i=1:length(takeBottom)
    currs=takeBottom(i);
    useTrials=ismember(trialLED,usel) & ismember(trialStimcond,currs);
    temp=psths(useTrials,psth.t>=useWindow(1) & psth.t<=useWindow(2));
    spontBottom(i,:)=nanmean(psths(useTrials,psth.t<=useWindow(1)),1);
    shiftInTime=((F1phases(currs)+pi)/(2*pi))*refPeriod;
    startInd=length(psth.t)-ceil(shiftInTime/(psth.t(2)-psth.t(1)));
    endInd=startInd+length(subt)-1;
    if isnan(F1phases(currs))
        bottomAv(i,:)=nan;
    else
        bottomAv(i,startInd:endInd)=nanmean(temp,1);
    end
end
  
ds=2;
if showFigs==1
    figure();
end
tAv=nanmean(topAv,1);
tAv=tAv(~isnan(tAv));
bAv=nanmean(bottomAv,1);
bAv=bAv(~isnan(bAv));
x1=[0:psth.t(2)-psth.t(1):(psth.t(2)-psth.t(1))*(size(spontTop,2)-1) (psth.t(2)-psth.t(1))*(size(spontTop,2)):psth.t(2)-psth.t(1):(psth.t(2)-psth.t(1))*(size(spontTop,2))+(psth.t(2)-psth.t(1))*(length(tAv)-1)];
y1=[nanmean(spontTop,1) tAv];
if showFigs==1
    plot(downSampAv(x1,ds),downSampAv(y1,ds),'Color','r');
    hold on;
end
x2=[0:psth.t(2)-psth.t(1):(psth.t(2)-psth.t(1))*(size(spontBottom,2)-1) (psth.t(2)-psth.t(1))*(size(spontBottom,2)):psth.t(2)-psth.t(1):(psth.t(2)-psth.t(1))*(size(spontBottom,2))+(psth.t(2)-psth.t(1))*(length(bAv)-1)];
y2=[nanmean(spontBottom,1) bAv];
if showFigs==1
    plot(downSampAv(x2,ds),downSampAv(y2,ds),'Color','b');
end
pref.x=x1;
pref.y=y1;
nonpref.x=x2;
nonpref.y=y2;

% plot(0:psth.t(2)-psth.t(1):(psth.t(2)-psth.t(1))*(length(bAv)-1),bAv,'Color','b');
    
    