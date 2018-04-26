function [LFPbySweep,Fs,ledConds,stimConds]=plotLFPFollowing(expt,fileInd,LFPbySweep,Fs)

timeWindow=[1 3];
stimWindow=timeWindow;
baseWindow=[0 1];
doCrossCorr=1;

dataDir='W:\New Acquisition Computer\';
if isempty(Fs)
    Fs=32000;
end
downSampFactor=10;
physChs=[23 15 22 14 19 11 20 12 17 9 16 8 18 10 21 13];
choosePhysCh=[6 7 8 9 10];
physCh=physChs(choosePhysCh);

useLEDcond1=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
useLEDcond2=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60]+0.05;
useTheseStimConds=[1];

daqFileNames=expt.files.names(fileInd);
disp('Using these daq files');
disp(daqFileNames);

if isempty(LFPbySweep)
    [LFPbySweep,Fs]=getJustLFPbySweep(dataDir,daqFileNames,Fs,downSampFactor,physCh);
end

% Match acquired LFP sweeps to LED conditions
ledConds=expt.sweeps.led(ismember(expt.sweeps.fileInd,fileInd));
if length(ledConds)~=size(LFPbySweep{1},1)
    disp('ledConds and acquired LFP sweeps do not match');
    return
end
% Get stimconds
stimConds=expt.sweeps.stimcond(ismember(expt.sweeps.fileInd,fileInd));
% Get LFP sweeps with chosen LED cond and stim. cond
for i=1:length(LFPbySweep)
    temp=LFPbySweep{i};
    L1{i}=temp(ismember(ledConds,useLEDcond1) & ismember(stimConds,useTheseStimConds),:);
    L2{i}=temp(ismember(ledConds,useLEDcond2) & ismember(stimConds,useTheseStimConds),:);
end
l=LFPbySweep{1};

figure(); 
xpoints=0:1/Fs:(size(L1{1},2)-1)*(1/Fs);
plot(xpoints,mean(L1{1},1),'Color','k');
hold on; 
plot(xpoints,mean(L2{1},1),'Color','r');
l1=L1{1};
return
[f1,avSpec1]=makePowerSpectrum(l1(:,xpoints>=timeWindow(1) & xpoints<=timeWindow(2)),Fs);
l2=L2{1};
[f2,avSpec2]=makePowerSpectrum(l2(:,xpoints>=timeWindow(1) & xpoints<=timeWindow(2)),Fs);
figure(); 
plot(f1,avSpec1,'Color','k');
hold on; 
plot(f2,avSpec2,'Color','r');

freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];

diffs=zeros(length(freqs),1);
cons=zeros(length(freqs),1);
leds=zeros(length(freqs),1);

for i=1:length(freqs)
    if doCrossCorr==1
        avCon=mean(l(ismember(ledConds,freqs(i)) & ismember(stimConds,useTheseStimConds),:),1);
        avLed=mean(l(ismember(ledConds,freqs(i)+0.05) & ismember(stimConds,useTheseStimConds),:),1);
        x=xpoints;
        psthEv.x=x(x>=stimWindow(1) & x<=stimWindow(2))-x(find(x>=stimWindow(1),1,'first'));
        % First control
        y=avCon;
        psthEv.y=y(x>=stimWindow(1) & x<=stimWindow(2));
        psthSp.x=x(x>=baseWindow(1) & x<baseWindow(2))-x(find(x>=baseWindow(1),1,'first'));
        psthSp.y=y(x>=baseWindow(1) & x<baseWindow(2));
        stim.x=psthEv.x;
        stim.y=sin(2*pi*floor(freqs(i)).*stim.x);
        [cc,cc_lags]=xcorr(stim.y-mean(stim.y),psthEv.y-mean(psthEv.y),[],'coeff'); % first input lags second input
        if length(stim.y(1:floor(end/2)))<length(psthSp.y)
            psthSp.y=psthSp.y(1:length(stim.y(1:floor(end/2))));
        end
        [ccSpont,cc_lagsSpont]=xcorr(stim.y(1:floor(end/2))-mean(stim.y(1:floor(end/2))),psthSp.y-mean(psthSp.y),[],'coeff');
        % neural response must lag stimulus
        % only consider as valid lags > 0
        new_ccX=cc_lags.*(stim.x(2)-stim.x(1));
        new_ccXSpont=cc_lagsSpont.*(stim.x(2)-stim.x(1));
        sampRate=new_ccX(2)-new_ccX(1);
        spontSampRate=new_ccXSpont(2)-new_ccXSpont(1);
        [f,avSpec]=makePowerSpectrum(cc(new_ccX>0 & new_ccX<=1),1/sampRate);
        [fSpont,avSpecSpont]=makePowerSpectrum(ccSpont(new_ccXSpont>0 & new_ccXSpont<=1),1/spontSampRate);
        cons(i)=mean(avSpec(f>freqs(i)-0.5 & f<freqs(i)+0.5))-mean(avSpecSpont(fSpont>freqs(i)-0.5 & fSpont<freqs(i)+0.5));
        % Then LED
        y=avLed;
        psthEv.y=y(x>=stimWindow(1) & x<=stimWindow(2));
        psthSp.x=x(x>=baseWindow(1) & x<baseWindow(2))-x(find(x>=baseWindow(1),1,'first'));
        psthSp.y=y(x>=baseWindow(1) & x<baseWindow(2));
        stim.x=psthEv.x;
        stim.y=sin(2*pi*floor(freqs(i)).*stim.x);
        [cc,cc_lags]=xcorr(stim.y-mean(stim.y),psthEv.y-mean(psthEv.y),[],'coeff'); % first input lags second input
        if length(stim.y(1:floor(end/2)))<length(psthSp.y)
            psthSp.y=psthSp.y(1:length(stim.y(1:floor(end/2))));
        end
        [ccSpont,cc_lagsSpont]=xcorr(stim.y(1:floor(end/2))-mean(stim.y(1:floor(end/2))),psthSp.y-mean(psthSp.y),[],'coeff');
        % neural response must lag stimulus
        % only consider as valid lags > 0
        new_ccX=cc_lags.*(stim.x(2)-stim.x(1));
        new_ccXSpont=cc_lagsSpont.*(stim.x(2)-stim.x(1));
        sampRate=new_ccX(2)-new_ccX(1);
        spontSampRate=new_ccXSpont(2)-new_ccXSpont(1);
        [f,avSpec]=makePowerSpectrum(cc(new_ccX>0 & new_ccX<=1),1/sampRate);
        [fSpont,avSpecSpont]=makePowerSpectrum(ccSpont(new_ccXSpont>0 & new_ccXSpont<=1),1/spontSampRate);
        leds(i)=mean(avSpec(f>freqs(i)-0.5 & f<freqs(i)+0.5))-mean(avSpecSpont(fSpont>freqs(i)-0.5 & fSpont<freqs(i)+0.5));
        diffs(i)=cons(i)-leds(i);
    else
        %     [f,avSpec]=makePowerSpectrum(l(ismember(ledConds,useLEDcond1) & ismember(stimConds,useTheseStimConds),xpoints>=timeWindow(1) & xpoints<=timeWindow(2)),Fs);
        %     [f2,avSpec2]=makePowerSpectrum(l(ismember(ledConds,useLEDcond2) & ismember(stimConds,useTheseStimConds),xpoints>=timeWindow(1) & xpoints<=timeWindow(2)),Fs);
        [f,avSpec]=makePowerSpectrum(l(ismember(ledConds,freqs(i)) & ismember(stimConds,useTheseStimConds),xpoints>=timeWindow(1) & xpoints<=timeWindow(2)),Fs);
        [f2,avSpec2]=makePowerSpectrum(l(ismember(ledConds,freqs(i)+0.05) & ismember(stimConds,useTheseStimConds),xpoints>=timeWindow(1) & xpoints<=timeWindow(2)),Fs);
%             [f,avSpec]=makePowerSpectrum(mean(l(ismember(ledConds,freqs(i)) & ismember(stimConds,useTheseStimConds),xpoints>=timeWindow(1) & xpoints<=timeWindow(2)),1),Fs);
%             [f2,avSpec2]=makePowerSpectrum(mean(l(ismember(ledConds,freqs(i)+0.05) & ismember(stimConds,useTheseStimConds),xpoints>=timeWindow(1) & xpoints<=timeWindow(2)),1),Fs);
        %     nonSpecificCon=mean(avSpec(~(f>freqs(i)-0.5 & f<freqs(i)+0.5)));
        %     nonSpecificLed=mean(avSpec2(~(f2>freqs(i)-0.5 & f2<freqs(i)+0.5)));
        nonSpecificCon=mean([avSpec(f>freqs(i)-1.5 & f<freqs(i)-0.5); avSpec(f>freqs(i)+0.5 & f<freqs(i)+1.5)]);
        nonSpecificLed=mean([avSpec2(f2>freqs(i)-1.5 & f2<freqs(i)-0.5); avSpec2(f2>freqs(i)+0.5 & f2<freqs(i)+1.5)]);
        conAmp=mean(avSpec(f>freqs(i)-0.5 & f<freqs(i)+0.5))-nonSpecificCon;
        ledAmp=mean(avSpec2(f2>freqs(i)-0.5 & f2<freqs(i)+0.5))-nonSpecificLed;
        cons(i)=conAmp;
        leds(i)=ledAmp;
        diffs(i)=conAmp-ledAmp;
    end
end

figure();
plot(freqs,cons,'Color','k');
hold on; 
plot(freqs,leds,'Color','r');
plot(freqs,diffs,'Color','b');