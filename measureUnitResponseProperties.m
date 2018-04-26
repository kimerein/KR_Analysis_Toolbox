function measureUnitResponseProperties(spikes,unitAssigns)




end

function [x,psths_t,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED]=getTrialByTrialUnitPSTH_sub(spikes,allAssigns,useLED,bin,trialDuration,useStimcond)

% useStimcond={[1:10000]};
% % useLED=[0 5];
% useLED=[0 5.05 5 -10];
% bin=1.2;
% trialDuration=10.5;

% Get trials for each unit
unitByUnitTrials=cell(1,length(allAssigns));
unitByUnitStimcond=cell(1,length(allAssigns));
unitByUnitLED=cell(1,length(allAssigns));
tt=unique(spikes.sweeps.trials);
if tt(1)>1 && length(spikes.sweeps.trials)==length(spikes.sweeps.stimcond)
    spikes.sweeps.trials=spikes.sweeps.trials-tt(1)+1;
    spikes.trials=spikes.trials-tt(1)+1;
end
if any(isnan(spikes.sweeps.led))
    spikes.sweeps.led(isnan(spikes.sweeps.led))=-10;
    spikes.led(isnan(spikes.led))=-10;
end
if any(isnan(spikes.sweeps.stimcond))
    spikes.sweeps.stimcond(isnan(spikes.sweeps.stimcond))=-10;
    spikes.stimcond(isnan(spikes.stimcond))=-10;
end

for i=1:length(allAssigns)
    unitByUnitTrials{i}=unique(spikes.sweeps.trials);
    unitByUnitStimcond{i}=spikes.sweeps.stimcond(unique(spikes.sweeps.trials));
    unitByUnitLED{i}=spikes.sweeps.led(unique(spikes.sweeps.trials));
end

psths_t=cell(length(allAssigns),length(useStimcond));
x=[];
for i=1:length(allAssigns)
    for j=1:length(useStimcond)
        useSpikes=filtspikes(spikes,0,'assigns',allAssigns(i),'stimcond',useStimcond{j});
        currTrialCon=unitByUnitTrials{i};
        currStimCon=unitByUnitStimcond{i};
        currLEDCon=unitByUnitLED{i};
        if ~isempty(useLED)
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{j}) & ismember(currLEDCon,useLED));
        else
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{j}));
        end
        if isempty(unitByUnitConsensus)
            disp('problem with unitByUnitConsensus');
        end
        [~,~,~,x,~,psths_t{i,j}]=psth_wStd_trialByTrial(useSpikes,bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
    end
end

end