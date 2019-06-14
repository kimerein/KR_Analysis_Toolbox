function [altogether,LEDalign_top_first,LEDalign_bottom_first,LEDalign_top_second,LEDalign_bottom_second,firstF1amps,secondF1amps,ampOrder,F1phases]=plotPref_wrapper(psth,refFreqBand,useWindow,useWindowForAmp,ampOrder,F1phases)

ledOn=5.05;
ledOff=0;
trialDuration=14.5;

u=1:length(psth.psths);

ampOrder_backup=ampOrder;

nTrials=size(psth.psths{1},1);
firstN=floor(nTrials/2);
% lastN=nTrials-firstN;
lastN=firstN;
psthFirst.t=psth.t;
psthFirst.psths=cell(length(u),1);
psthFirst.unitTrials=cell(1,length(u));
psthFirst.unitStimcond=cell(1,length(u));
psthFirst.unitLED=cell(1,length(u));
psthSecond.t=psth.t;
psthSecond.psths=cell(length(u),1);
psthSecond.unitTrials=cell(1,length(u));
psthSecond.unitStimcond=cell(1,length(u));
psthSecond.unitLED=cell(1,length(u));
for i=1:length(u)
    temp=psth.psths{i};
    psthFirst.psths{i}=temp(1:firstN,:);
    psthSecond.psths{i}=temp(end-lastN+1:end,:);
    temp=psth.unitTrials{i};
    psthFirst.unitTrials{i}=temp(1:firstN);
    psthSecond.unitTrials{i}=temp(end-lastN+1:end);
    temp=psth.unitStimcond{i};
    psthFirst.unitStimcond{i}=temp(1:firstN);
    psthSecond.unitStimcond{i}=temp(end-lastN+1:end);
    temp=psth.unitLED{i};
    psthFirst.unitLED{i}=temp(1:firstN);
    psthSecond.unitLED{i}=temp(end-lastN+1:end);
end

nTimes=length(0:0.01:trialDuration);
allTop_PSTH.x=nan(length(u),nTimes);
allTop_PSTH.y=nan(length(u),nTimes);
allBottom_PSTH.x=nan(length(u),nTimes);
allBottom_PSTH.y=nan(length(u),nTimes);
firstTop_PSTH.x=nan(length(u),nTimes);
firstTop_PSTH.y=nan(length(u),nTimes);
firstBottom_PSTH.x=nan(length(u),nTimes);
firstBottom_PSTH.y=nan(length(u),nTimes);
secondTop_PSTH.x=nan(length(u),nTimes);
secondTop_PSTH.y=nan(length(u),nTimes);
secondBottom_PSTH.x=nan(length(u),nTimes);
secondBottom_PSTH.y=nan(length(u),nTimes);

allTop_PSTH_LED.x=nan(length(u),nTimes);
allTop_PSTH_LED.y=nan(length(u),nTimes);
allBottom_PSTH_LED.x=nan(length(u),nTimes);
allBottom_PSTH_LED.y=nan(length(u),nTimes);
firstTop_PSTH_LED.x=nan(length(u),nTimes);
firstTop_PSTH_LED.y=nan(length(u),nTimes);
firstBottom_PSTH_LED.x=nan(length(u),nTimes);
firstBottom_PSTH_LED.y=nan(length(u),nTimes);
secondTop_PSTH_LED.x=nan(length(u),nTimes);
secondTop_PSTH_LED.y=nan(length(u),nTimes);
secondBottom_PSTH_LED.x=nan(length(u),nTimes);
secondBottom_PSTH_LED.y=nan(length(u),nTimes);
firstF1amps=nan(length(u),1);
secondF1amps=nan(length(u),1);
for i=1:length(u)
    ampOrder=ampOrder_backup;
    if isempty(ampOrder)
        [topAv,bottomAv,pref,nonpref,ampOrder,F1phases,~,takeTop]=plotPrefVersusNonpref(psth,i,refFreqBand,useWindow,useWindowForAmp,[],[],ledOff);
    else
        [topAv,bottomAv,pref,nonpref,ampOrder,F1phases,~,takeTop]=plotPrefVersusNonpref(psth,i,refFreqBand,useWindow,useWindowForAmp,ampOrder,F1phases,ledOff);
    end
    temp=psthFirst.psths{i};
    uStimcond=psthFirst.unitStimcond{i};
    uLed=psthFirst.unitLED{i};
    takeN=1;
    LEDalign_top_first(i,:)=nanmean(temp(ismember(uStimcond,takeTop) & ismember(uLed,ledOff),:),1);
%     LEDalign_top_first(i,:)=nanmean(temp(ismember(uStimcond,takeTop) & ismember(uLed,ledOn),:),1);
%     LEDalign_top_first(i,:)=nanmean(temp(ismember(uStimcond,ampOrder(end-takeN+1:end)) & ismember(uLed,ledOn),:),1);
    LEDalign_bottom_first(i,:)=nanmean(temp(ismember(uStimcond,ampOrder(1:takeN)) & ismember(uLed,ledOn),:),1);
    temp=psthSecond.psths{i};
    uStimcond=psthSecond.unitStimcond{i};
    uLed=psthSecond.unitLED{i};
    LEDalign_top_second(i,:)=nanmean(temp(ismember(uStimcond,ampOrder(end-takeN+1:end)) & ismember(uLed,ledOn),:),1);
    LEDalign_bottom_second(i,:)=nanmean(temp(ismember(uStimcond,ampOrder(1:takeN)) & ismember(uLed,ledOn),:),1);
    allTop_PSTH.x(i,:)=interp1(pref.x,pref.x,0:0.01:trialDuration);
    allTop_PSTH.y(i,:)=interp1(pref.x,pref.y,0:0.01:trialDuration);
    allBottom_PSTH.x(i,:)=interp1(nonpref.x,nonpref.x,0:0.01:trialDuration);
    allBottom_PSTH.y(i,:)=interp1(nonpref.x,nonpref.y,0:0.01:trialDuration);
    [topAv,bottomAv,pref,nonpref,~,~,highestAmp]=plotPrefVersusNonpref(psthFirst,i,refFreqBand,useWindow,useWindowForAmp,ampOrder,F1phases,ledOff);
    firstF1amps(i)=highestAmp;
    firstTop_PSTH.x(i,:)=interp1(pref.x,pref.x,0:0.01:trialDuration);
    firstTop_PSTH.y(i,:)=interp1(pref.x,pref.y,0:0.01:trialDuration);
    firstBottom_PSTH.x(i,:)=interp1(nonpref.x,nonpref.x,0:0.01:trialDuration);
    firstBottom_PSTH.y(i,:)=interp1(nonpref.x,nonpref.y,0:0.01:trialDuration);
    [topAv,bottomAv,pref,nonpref,~,~,highestAmp]=plotPrefVersusNonpref(psthSecond,i,refFreqBand,useWindow,useWindowForAmp,ampOrder,F1phases,ledOff);
    secondF1amps(i)=highestAmp;
    secondTop_PSTH.x(i,:)=interp1(pref.x,pref.x,0:0.01:trialDuration);
    secondTop_PSTH.y(i,:)=interp1(pref.x,pref.y,0:0.01:trialDuration);
    secondBottom_PSTH.x(i,:)=interp1(nonpref.x,nonpref.x,0:0.01:trialDuration);
    secondBottom_PSTH.y(i,:)=interp1(nonpref.x,nonpref.y,0:0.01:trialDuration);
    [topAv,bottomAv,pref,nonpref]=plotPrefVersusNonpref(psth,i,refFreqBand,useWindow,useWindowForAmp,ampOrder,F1phases,ledOn);
    allTop_PSTH_LED.x(i,:)=interp1(pref.x,pref.x,0:0.01:trialDuration);
    allTop_PSTH_LED.y(i,:)=interp1(pref.x,pref.y,0:0.01:trialDuration);
    allBottom_PSTH_LED.x(i,:)=interp1(nonpref.x,nonpref.x,0:0.01:trialDuration);
    allBottom_PSTH_LED.y(i,:)=interp1(nonpref.x,nonpref.y,0:0.01:trialDuration);
    [topAv,bottomAv,pref,nonpref]=plotPrefVersusNonpref(psthFirst,i,refFreqBand,useWindow,useWindowForAmp,ampOrder,F1phases,ledOn);
    firstTop_PSTH_LED.x(i,:)=interp1(pref.x,pref.x,0:0.01:trialDuration);
    firstTop_PSTH_LED.y(i,:)=interp1(pref.x,pref.y,0:0.01:trialDuration);
    firstBottom_PSTH_LED.x(i,:)=interp1(nonpref.x,nonpref.x,0:0.01:trialDuration);
    firstBottom_PSTH_LED.y(i,:)=interp1(nonpref.x,nonpref.y,0:0.01:trialDuration);
    [topAv,bottomAv,pref,nonpref]=plotPrefVersusNonpref(psthSecond,i,refFreqBand,useWindow,useWindowForAmp,ampOrder,F1phases,ledOn);
    secondTop_PSTH_LED.x(i,:)=interp1(pref.x,pref.x,0:0.01:trialDuration);
    secondTop_PSTH_LED.y(i,:)=interp1(pref.x,pref.y,0:0.01:trialDuration);
    secondBottom_PSTH_LED.x(i,:)=interp1(nonpref.x,nonpref.x,0:0.01:trialDuration);
    secondBottom_PSTH_LED.y(i,:)=interp1(nonpref.x,nonpref.y,0:0.01:trialDuration);
end
    
altogether.allTop_PSTH=allTop_PSTH;
altogether.allBottom_PSTH=allBottom_PSTH;
altogether.firstTop_PSTH=firstTop_PSTH;
altogether.firstBottom_PSTH=firstBottom_PSTH;
altogether.secondTop_PSTH=secondTop_PSTH;
altogether.secondBottom_PSTH=secondBottom_PSTH;
altogether.allTop_PSTH_LED=allTop_PSTH_LED;
altogether.allBottom_PSTH_LED=allBottom_PSTH_LED;
altogether.firstTop_PSTH_LED=firstTop_PSTH_LED;
altogether.firstBottom_PSTH_LED=firstBottom_PSTH_LED;
altogether.secondTop_PSTH_LED=secondTop_PSTH_LED;
altogether.secondBottom_PSTH_LED=secondBottom_PSTH_LED;

figure(); 
subplot(2,1,1);
plot(nanmean(allTop_PSTH.x,1),nanmean(allTop_PSTH.y,1),'Color','r');
hold on;
plot(nanmean(allBottom_PSTH.x,1),nanmean(allBottom_PSTH.y,1),'Color','b');
subplot(2,1,2);
plot(nanmean(allTop_PSTH_LED.x,1),nanmean(allTop_PSTH_LED.y,1),'Color','r');
hold on;
plot(nanmean(allBottom_PSTH_LED.x,1),nanmean(allBottom_PSTH_LED.y,1),'Color','b');
title('All Trials');

figure(); 
subplot(2,1,1);
plot(nanmean(firstTop_PSTH.x,1),nanmean(firstTop_PSTH.y,1),'Color','r');
hold on;
plot(nanmean(firstBottom_PSTH.x,1),nanmean(firstBottom_PSTH.y,1),'Color','b');
subplot(2,1,2);
plot(nanmean(firstTop_PSTH_LED.x,1),nanmean(firstTop_PSTH_LED.y,1),'Color','r');
hold on;
plot(nanmean(firstBottom_PSTH_LED.x,1),nanmean(firstBottom_PSTH_LED.y,1),'Color','b');
title('First Half');

figure(); 
subplot(2,1,1);
plot(nanmean(secondTop_PSTH.x,1),nanmean(secondTop_PSTH.y,1),'Color','r');
hold on;
plot(nanmean(secondBottom_PSTH.x,1),nanmean(secondBottom_PSTH.y,1),'Color','b');
subplot(2,1,2);
plot(nanmean(secondTop_PSTH_LED.x,1),nanmean(secondTop_PSTH_LED.y,1),'Color','r');
hold on;
plot(nanmean(secondBottom_PSTH_LED.x,1),nanmean(secondBottom_PSTH_LED.y,1),'Color','b');
title('Second Half');