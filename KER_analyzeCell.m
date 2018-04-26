function KER_analyzeCell(expt,spikes)
% Problem with this function: this function calls the non-existent function
% KER_cellAnalysisParams.m
% Analyze a unit of expt and make figures about the firing of this unit
% 
% PARAMETERS:
% expt: the expt structure from which this unit is pulled
% spikes: the spikes structure

global spikesDir;

spikesDir='C:\Documents and Settings\Admin\My Documents\MATLAB\Data\Analyzed\SortedSpikes';


% Initialize parameters for analysis
[assigns,whichDaqFiles,cellName,ONlength,ONstart,totalTrialLength,dateName,nOrients,nPhases]=KER_cellAnalysisParams();
daqFileNames=expt.files.names(whichDaqFiles);
saveToDir=strcat(strcat('C:\Documents and Settings\Admin\My Documents\Recurrent Activity - Data\Data from Expt. on',strcat(dateName,'\')),cellName);
if ~exist(saveToDir,'dir')
    mkdir(saveToDir);
end

LFPbySweep=[];
tempspikes=filtspikes(spikes,0,'assigns',assigns);

% Make PSTH plots
KER_PSTH(tempspikes,totalTrialLength,nOrients,'Orientation',nPhases,'Phase',ONstart,ONstart+ONlength,ONstart+ONlength,ONstart+ONlength+1);

% Show histogram of OFF period spiking from trial to trial - idea is to
% search for global oscillations in brain state
countTrialByTrialOFFspikes=zeros(max(tempspikes.trials),1);
for i=1:max(tempspikes.trials)
    someSpikes=filtspikes(tempspikes,0,'trials',i);
    countTrialByTrialOFFspikes(i)=length(someSpikes.spiketimes(someSpikes.spiketimes>ONstart+ONlength & someSpikes.spiketimes<=totalTrialLength));
end
[h,x]=hist(countTrialByTrialOFFspikes,max(countTrialByTrialOFFspikes)+1);
figure;
bar(x,h);
redThresh=input('Please enter OFF period spiking threshold for red trials. All trials with >= this number of OFF spikes will be red.'); 

% Plot raster sorted by stimulus condition
KER_SortedRaster(tempspikes,[],nOrients,nPhases,totalTrialLength,1);

% Find red and blue trials for separate LFP analysis (idea is to separate
% stimulus-triggered LFP traces by global brain state)
redTrials=[];
blueTrials=[];
for i=1:max(spikes.trials)
    someSpikes=filtspikes(tempspikes,0,'trials',i);
    countOFFspikes=length(someSpikes.spiketimes(someSpikes.spiketimes>ONstart+ONlength & someSpikes.spiketimes<=totalTrialLength));
    if countOFFspikes>=redThresh
        redTrials=[redTrials i];
    else
        blueTrials=[blueTrials i];
    end
end

% Only compare averages for red and blue brain states for the same number
% of trials in each average
sweepInds=randperm(size(LFPbySweep,1));
LFPbySweep=LFPbySweep(sweepInds(1:chooseNSweeps),:);
if length(redTrials)>length(blueTrials)
    takeInds=randperm(length(redTrials));
    redTrials=redTrials(sort(takeInds(1:length(blueTrials))));
elseif length(blueTrials)>length(redTrials)
    takeInds=randperm(length(blueTrials));
    blueTrials=blueTrials(sort(takeInds(1:length(redTrials))));
end

% Compare stimulus-triggered LFP traces for red and blue brain states
figure;
[LFPbySweep,Fs,bandPassedLFPbySweep]=KER_analyzeLFP(tempspikes,daqFileNames,LFPbySweep,redTrials,blueTrials,ONlength,ONstart,totalTrialLength);

% Draw stimulus-triggered LFPs for all spikes
KER_spikeTriggeredLFP(tempspikes,LFPbySweep,Fs,0,totalTrialLength);
KER_spikeTriggeredLFP(tempspikes,bandPassedLFPbySweep,Fs,0,totalTrialLength);
% Draw stimulus-triggered LFPs for ON period spikes
KER_spikeTriggeredLFP(tempspikes,LFPbySweep,Fs,ONstart,ONlength);
KER_spikeTriggeredLFP(tempspikes,bandPassedLFPbySweep,Fs,ONstart,ONlength);
% Draw stimulus-triggered LFPs for OFF period spikes
KER_spikeTriggeredLFP(tempspikes,LFPbySweep,Fs,ONstart+ONlength,totalTrialLength);
KER_spikeTriggeredLFP(tempspikes,bandPassedLFPbySweep,Fs,ONstart+ONlength,totalTrialLength);

% Save results of analysis
saveToDir=strcat(saveToDir,'\');
saveas(figure(1),strcat(saveToDir,'PSTH_Total.bmp'));
saveas(figure(2),strcat(saveToDir,'PSTH_Orientation&Phase.bmp'));
saveas(figure(3),strcat(saveToDir,'PSTH_Orientation.bmp'));
saveas(figure(4),strcat(saveToDir,'orientationTuning.bmp'));
saveas(figure(5),strcat(saveToDir,'hist_OFFspikes.bmp'));
saveas(figure(6),strcat(saveToDir,'colored_raster_redGreaterThanEqualTo',num2str(redThresh),'OFFspikes.bmp'));
saveas(figure(7),strcat(saveToDir,'stimTriggeredLFP.bmp'));
saveas(figure(9),strcat(saveToDir,'stimTriggeredLFP_significanceOfGreaterOFFspiking.bmp'));
saveas(figure(10),strcat(saveToDir,'stimTriggeredLFP_filtered30to80Hz.bmp'));
saveas(figure(13),strcat(saveToDir,'stimTriggeredLFP_filtered30to80Hz_ampDiffWithGreaterThanEqualTo',num2str(redThresh),'OFFSpikes.bmp'));
saveas(figure(14),strcat(saveToDir,'spikeTriggeredLFP.bmp'));
saveas(figure(15),strcat(saveToDir,'spikeTriggered_bandPassed30to80LFP.bmp'));
saveas(figure(16),strcat(saveToDir,'spikeTriggeredLFP_ONperiod.bmp'));
saveas(figure(17),strcat(saveToDir,'spikeTriggeredLFP_bandPassed30to80_ONperiod.bmp'));
saveas(figure(18),strcat(saveToDir,'spikeTriggeredLFP_OFFperiod.bmp'));
saveas(figure(19),strcat(saveToDir,'spikeTriggeredLFP_bandPassed30to80_OFFperiod.bmp'));
saveas(figure(20),strcat(saveToDir,'raster_sortedByStimCond.bmp'));
end