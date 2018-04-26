function [spontChangeWithinBlock,trialTally]=scriptForMeasuringSpontChange(psth,trialTally,isTheta)

% trialTally=1;
firstHalf=1:floor(size(psth.psths{1},1)./2); 
secondHalf=floor(size(psth.psths{1},1)./2)+1:size(psth.psths{1},1);
firstHalf=1:size(psth.psths{1},1);
for i=1:length(psth.psths)
p=psth.psths{i};
l=psth.unitLED{i};
subtheta=isTheta(trialTally:trialTally+size(psth.psths{1},1)-1);
p=p(firstHalf,:);
l=l(firstHalf);
sT=subtheta(firstHalf);
currSpont(i)=nanmean(nanmean(p(ismember(l,0) & ~sT,psth.t<=4),2),1);
end
firstHalf=1:floor(size(psth.psths{1},1)./2); secondHalf=floor(size(psth.psths{1},1)./2)+1:size(psth.psths{1},1);
for i=1:length(psth.psths)
p=psth.psths{i};
l=psth.unitLED{i};
subtheta=isTheta(trialTally:trialTally+size(psth.psths{1},1)-1);
p=p(secondHalf,:);
l=l(secondHalf);
sT=subtheta(secondHalf);
currSpont_secondHalf(i)=nanmean(nanmean(p(ismember(l,0) & ~sT,psth.t<=4),2),1);
end
trialTally=trialTally+size(psth.psths{1},1);
spontChangeWithinBlock=currSpont-currSpont_secondHalf;
spontChangeWithinBlock=currSpont;