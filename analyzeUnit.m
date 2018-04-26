function [figs,unitAnalysisInfo]=analyzeUnit(expt,spikes,params,LFP,bandPassedLFP,LFP_Fs,unitAnalysisInfo)
% Analyzes a unit
% 
% RETURNS:
% figs: handles to analysis figures
% unitAnalysisInfo: a structure with information about analysis choices
% 
% PARAMETERS:
% expt: the expt structure currently being analyzed
% spikes: the spikes structure for this unit
% params: a structure with information about the setup of a trial
% LFP: a matrix with samples as columns and rows as trials
% bandPassedLFP: same as LFP, but represents band-passed data
% LFP_Fs: sampling rate of LFP, bandPassedLFP
% unitAnalysisInfo: a structure containing information about analysis
% choices for this unit

figs=[];

% Make and return handles to PSTH plots (also polar plots for oriented
% stimuli)
quantResponse='integral';       % choose method for quantifying ON and OFF responses
ONresponseWindow=[params.ONstart params.ONstart+params.ONlength];
OFFresponseWindow=[params.ONstart+params.ONlength params.totalTrialLength];
[psth_handles,psthAxes_handles,colorMapCode]=plotPSTH(spikes,params,ONresponseWindow,OFFresponseWindow,quantResponse);
% Order of returned figures in psth_handles:
% Fig. 1: Total PSTH; PSTH collapsed over Var 1; PSTH collapsed over Var 2
% Fig. 2: Separate PSTH plots for each stimulus conditions, also containing
% ON and OFF responses for different stimulus conditions
% Fig. 3: Orientation polar plots for ON and OFF responses
figs=[figs; psth_handles];

% Show histogram of OFF period spiking from trial to trial - idea is to
% distinguish global oscillations in brain state
% Then ask user to select a threshold of OFF period spiking to distinguish
% 2 brain states (e.g., "UP", that is, more spiking, and "DOWN", less spiking, states)
[spikingThresh,redTrials,blueTrials]=defineBrainStates(spikes,params);
unitAnalysisInfo.definedRedTrials=redTrials;
unitAnalysisInfo.definedBlueTrials=blueTrials;

% Make and return handle to raster plots:
%  - Raster sorted by stimulus condition
%  - Raster in chronological trial order 
%    with red and blue trials, as defined in defineBrainStates
fig=figure;
axesmatrix(2,1,1);
r=sortRaster(spikes,params,0,colorMapCode,fig);
axesmatrix(2,1,2);
r2=twoColorRaster(spikes,redTrials,blueTrials);
figs=[figs; r];

% Only compare LFP averages for red and blue brain states 
% taking same number of trials for each average 
% so as to be able to compare variability
if isempty(unitAnalysisInfo.comparedRedTrials) || isempty(unitAnalysisInfo.comparedBlueTrials)
    if length(redTrials)>length(blueTrials)
        takeInds=randperm(length(redTrials));
        redTrials=redTrials(sort(takeInds(1:length(blueTrials))));
    else
        takeInds=randperm(length(blueTrials));
        blueTrials=blueTrials(sort(takeInds(1:length(redTrials))));
    end
else
    redTrials=unitAnalysisInfo.comparedRedTrials;
    blueTrials=unitAnalysisInfo.comparedBlueTrials;
end
    
% Spike-triggered LFP
% Make a sheet of spike-triggered LFP graphs
% All spikes vs. ON period spikes vs. OFF period spikes
% Raw LFP vs. band-passed LFP
[spikeTrigLFP_figs]=plotSpikeTrigLFPs(spikes,LFP,bandPassedLFP,LFP_Fs,params,redTrials,blueTrials);
figs=[figs; spikeTrigLFP_figs];
% Returns one figure containing ON and OFF spike-triggered LFP (for all,
% red, blue trials) & one figure containing Stimulus Var 1-specific
% spike-triggered LFPs

% Try different LFP phases?

% Return unit-specific analysis details and figures
unitAnalysisInfo.redTrial_spikingThresh=spikingThresh;
unitAnalysisInfo.responseQuantMethod=quantResponse;
unitAnalysisInfo.stimConds_colorCode=colorMapCode;
unitAnalysisInfo.comparedRedTrials=redTrials;
unitAnalysisInfo.comparedBlueTrials=blueTrials;
