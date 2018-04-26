function [allSpikes,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED]=prepare_compareUnitAverageTuningCurves(allSpikes,registerSpikes)

mapping1=1:32;
mapping2=[1 1 1 1 2 2 2 2 3 3 3 3 4 4 4 4 5 5 5 5 6 6 6 6 7 7 7 7 8 8 8 8];
% mapping1=1:8;
% mapping2=1:8;
shiftSpiketimes=0;
shiftSpiketimesBy=0.04; % in seconds

allSpikes.led=floor(allSpikes.led);
allSpikes.sweeps.led=floor(allSpikes.sweeps.led);
for i=1:length(mapping1)
    allSpikes.stimcond(allSpikes.stimcond==mapping1(i))=ones(size(allSpikes.stimcond(allSpikes.stimcond==mapping1(i)))).*mapping2(i);
    allSpikes.sweeps.stimcond(allSpikes.sweeps.stimcond==mapping1(i))=ones(size(allSpikes.sweeps.stimcond(allSpikes.sweeps.stimcond==mapping1(i)))).*mapping2(i);
end

if shiftSpiketimes==1
    allSpikes.spiketimes=allSpikes.spiketimes+shiftSpiketimesBy;
end

a=unique(allSpikes.assigns);
bytrials=registerSpikes.sweeps.trials;
for i=1:length(a)
    unitByUnitTrials{i}=bytrials;
end
bystimcond=registerSpikes.sweeps.stimcond;
for i=1:length(mapping1)
    bystimcond(bystimcond==mapping1(i))=ones(size(bystimcond(bystimcond==mapping1(i)))).*mapping2(i);
end
for i=1:length(a)
    unitByUnitStimcond{i}=bystimcond;
end
byled=registerSpikes.sweeps.led;
byled=floor(byled);
for i=1:length(a)
    unitByUnitLED{i}=byled;
end