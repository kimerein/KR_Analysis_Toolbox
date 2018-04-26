function outpsth=filtPSTH(psth,takeTrials)


outpsth.t=psth.t;
outpsth.psths=cell(length(psth.psths),1);
outpsth.unitTrials=cell(1,length(psth.psths));
outpsth.unitStimcond=cell(1,length(psth.psths));
outpsth.unitLED=cell(1,length(psth.psths));

tri=psth.unitTrials{1};
for i=1:length(psth.psths)
    p=psth.psths{i};
    outpsth.psths{i}=p(ismember(tri,takeTrials),:);
    unitTrials=psth.unitTrials{i};
    outpsth.unitTrials{i}=unitTrials(ismember(tri,takeTrials));
    unitStimcond=psth.unitStimcond{i};
    outpsth.unitStimcond{i}=unitStimcond(ismember(tri,takeTrials));
    unitLED=psth.unitLED{i};
    outpsth.unitLED{i}=unitLED(ismember(tri,takeTrials));
end
    