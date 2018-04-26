function outpsth=getBurstFraction(psth,burstpsth)

outpsth.t=psth.t;
outpsth.unitTrials=psth.unitTrials;
outpsth.unitStimcond=psth.unitStimcond;
outpsth.unitLED=psth.unitLED;

for i=1:length(psth.psths)
    p=burstpsth.psths{i}./psth.psths{i};
%     p(ismember(burstpsth.psths{i},0))=0;
%     p(~ismember(burstpsth.psths{i},0) & ~ismember(psth.psths{i},
%     p(isnan(p))=0;
    outpsth.psths{i}=p;
end