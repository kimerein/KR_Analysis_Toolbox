function plotStimsAndLEDCondsSeparateley(baselineWindow,LFPbySweep,nStimParam1Vals,nStimParam2Vals,ledForSweeps,stimsForSweeps,LFP_Fs,totalTrialLength)

% Discard first and last trials
LFPbySweep=LFPbySweep(2:end-1,:);
ledForSweeps=ledForSweeps(2:end-1);
stimsForSweeps=stimsForSweeps(2:end-1);

% Align trials
% baseInds=floor(baselineWindow(1)*LFP_Fs)+1:floor(baselineWindow(2)*LFP_Fs);
% for i=1:size(LFPBySweep,1)
%     baseVal=mean(LFPbySweep(i,baseInds));
%     LFPbySweep(i,:)=LFPbySweep(i,:)-baseVal;
% end

figure;

% Collapse over stim. variable 1
noLEDtrialsForConds=cell(length(nStimParam2Vals),1);
LEDtrialsForConds=cell(length(nStimParam2Vals),1);
for i=1:nStimParam2Vals
    theseConds=i:nStimParam2Vals:nStimParam1Vals*nStimParam2Vals;
    t=ismember(stimsForSweeps,theseConds);
    noLEDtrialsForConds{i}=find(t&(~(ledForSweeps>0)));
    LEDtrialsForConds{i}=find(t&(ledForSweeps>0));
end
noLEDtrialsForConds=matchManyNumberOfTrials(noLEDtrialsForConds);
LEDtrialsForConds=matchManyNumberOfTrials(LEDtrialsForConds);
disp(noLEDtrialsForConds);
disp(LEDtrialsForConds);
% LED OFF
subplot(2,1,1);
baseInds=floor(baselineWindow(1)*LFP_Fs)+1:floor(baselineWindow(2)*LFP_Fs);
for i=1:nStimParam2Vals
    baseVal=mean(mean(LFPbySweep(noLEDtrialsForConds{i},baseInds),1),2);
    LFPbySweep(noLEDtrialsForConds{i},:)=LFPbySweep(noLEDtrialsForConds{i},:)-baseVal;
    plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,mean(LFPbySweep(noLEDtrialsForConds{i},:),1),'Color','r');
    hold all;
end
axis([0 totalTrialLength  -0.3 0.05]);
hold off
% LED ON
subplot(2,1,2);
for i=1:nStimParam2Vals
    baseVal=mean(mean(LFPbySweep(LEDtrialsForConds{i},baseInds),1),2);
    LFPbySweep(LEDtrialsForConds{i},:)=LFPbySweep(LEDtrialsForConds{i},:)-baseVal;
    plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,mean(LFPbySweep(LEDtrialsForConds{i},:),1),'Color','b');
    hold all;
end
axis([0 totalTrialLength  -0.3 0.05]);
%title('Collapsed Over Stim. Variable 1');

figure;

% Collapse over stim. variable 2
noLEDtrialsForConds=cell(length(nStimParam2Vals),1);
LEDtrialsForConds=cell(length(nStimParam2Vals),1);
for i=1:nStimParam1Vals
    theseConds=i:1:i+nStimParam2Vals-1;
    t=ismember(stimsForSweeps,theseConds);
    noLEDtrialsForConds{i}=find(t&(~(ledForSweeps>0)));
    LEDtrialsForConds{i}=find(t&(ledForSweeps>0));
end
noLEDtrialsForConds=matchManyNumberOfTrials(noLEDtrialsForConds);
LEDtrialsForConds=matchManyNumberOfTrials(LEDtrialsForConds);
disp(noLEDtrialsForConds);
disp(LEDtrialsForConds);
% LED OFF
subplot(2,1,1);
for i=1:nStimParam1Vals
     baseVal=mean(mean(LFPbySweep(noLEDtrialsForConds{i},baseInds),1),2);
    LFPbySweep(noLEDtrialsForConds{i},:)=LFPbySweep(noLEDtrialsForConds{i},:)-baseVal;
    plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,mean(LFPbySweep(noLEDtrialsForConds{i},:),1),'Color','r');
    hold all;
end
axis([0 totalTrialLength  -0.3 0.05]);
hold off

% LED ON
subplot(2,1,2);
for i=1:nStimParam1Vals
    baseVal=mean(mean(LFPbySweep(LEDtrialsForConds{i},baseInds),1),2);
    LFPbySweep(LEDtrialsForConds{i},:)=LFPbySweep(LEDtrialsForConds{i},:)-baseVal;
    plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,mean(LFPbySweep(LEDtrialsForConds{i},:),1),'Color','b');
    hold all;
end
axis([0 totalTrialLength  -0.3 0.05]);
title('Collapsed Over Stim. Variable 2');