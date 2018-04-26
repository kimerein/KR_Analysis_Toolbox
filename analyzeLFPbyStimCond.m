function [LFPfigs,analysisInfo]=analyzeLFPbyStimCond(LFPdata,stimsForSweeps,ledForSweeps,alignToInitialOrMean,LFP_Fs,params,Av_Amp_Compare,compareSumOrAv,compareByStim_1,compareByStim_2,colorMapCode,rectifyForStats,statsBinWidth,showPhoto,useLED,photoSignal,ledSignal,showSpecgrams,useNTrials)
% Part of the automatic experiment analysis package
% Wrapper for this package is allAnalyses.m
% This function analyzes LFPdata and makes figures
% RETURNS:
% LFPfigs: handles to analysis figures
% analysisInfo: information structure about analysis choices
% 
% PARAMETERS:
% LFPdata: columns are samples, rows are trials
% stimsForSweeps: vector with stimulus condition number for each trial
% ledForSweeps: vector with LED value (e.g., step amplitude) for each trial
% alignToInitialOrMean: if 'initial', aligns all LFP trials to initial
%     values, if 'mean', aligns all trials to means
% LFP_Fs: sampling rate of LFPdata
% params: structure with information about stimulus periods 
% Av_Amp_Compare: three-element vector. if Av_Amp_Compare(3)=1, run
% compareLFPconds: runs statistics on LFP signals for 2 different
% stimulus conditions (need to implement code for pairwise stim. cond. comparisons;
% currently Av_Amp_Compare(1:2) do nothing
% compareSumOrAv: decides whether to sum over trials (pass in 'sum') 
% or average over trials (pass in 'average')
% compareByStim_1: two-element vector, specifies which stimulus condition
% to compare (statistically) with compareByStim_2
% first element is the stim. condition # for the first varying stimulus 
% parameter; second element is stim. cond. # for second varying stimulus
% parameter (e.g., if spatial freq. and phase are changing,
% compareByStim_1=[1 4] will specify the first spatial freq. cond., as indexed in
% expt, and the fourth phase condition)
% compareByStim_2: same as compareByStim_1 but is compared with
% compareByStim_1
% colorMapCode is a structure array specifying the plotting color associated with
% each stimulus/LED condition, see makeColorCode.m
% rectifyForStats: passed to compareLFPconds.m (which is a little obsolete)
% statsBinWidth: the binWidth for each ttest in indices (not currently
% used)
% showPhoto: if 1, show photoSignal trace at bottom of each graph
% useLED: if 1, separate stimulus conditions by whether LED on or off for
% analysis and plotting
% photoSignal: show this signal at bottom of each graph is useLED==1
% ledSignal: if useLED==1, show at bottom of each graph
% showSpecgrams: if 1, make spectrograms
% useNTrials: if showSpecgrams==1, useNTrials is number of trials to
% average to make each spectrogram

LFPfigs=[];
analysisInfo=[];

if strcmp(alignToInitialOrMean,'initial') % Align to initial value of each trial
    for i=1:size(LFPdata,1)
        initVal=LFPdata(i,1);
        LFPdata(i,:)=LFPdata(i,:)-initVal;
    end
elseif strcmp(alignToInitialOrMean,'mean') % Align to mean of each trial
    mVals=mean(LFPdata,2);
    for i=1:size(LFPdata,2)
        LFPdata(:,i)=LFPdata(:,i)-mVals;
    end
else % 'none'
end

if strcmp(compareSumOrAv,'average')
    [l,axes_handles,colorMapCode]=plotLFPconds(LFPdata,LFP_Fs,stimsForSweeps,ledForSweeps,params,'average',colorMapCode,1,1,statsBinWidth,showPhoto,useLED,photoSignal,ledSignal);
    LFPfigs=[LFPfigs; l];
else
    [l,axes_handles,colorMapCode]=plotLFPconds(LFPdata,LFP_Fs,stimsForSweeps,ledForSweeps,params,'add',colorMapCode,1,1,statsBinWidth,showPhoto,useLED,photoSignal,ledSignal);
    LFPfigs=[LFPfigs; l];
end

% Make spectrograms if requested
if showSpecgrams==1
    l=plotSpectrogramsByCondition(LFPdata,LFP_Fs,stimsForSweeps,ledForSweeps,params,1,useNTrials,[1 0 0]);
    LFPfigs=[LFPfigs; l];
end
    
% Compare signal responses for different stimulus conditions
% Currently, only compares two stimulus conditions
% Need to replace with code for pairwise comparisons (e.g., 
% calculateStimCondStats)
if Av_Amp_Compare(3)==1
    i1=compareByStim_1(1);
    j1=compareByStim_1(2);
    i2=compareByStim_2(1);
    j2=compareByStim_2(2);
    group1_values=(i1-1)*length(params.Var2_values)+j1;
    group2_values=(i2-1)*length(params.Var2_values)+j2;
    [l,redTrials,blueTrials]=compareLFPconds(LFPdata,stimsForSweeps,LFP_Fs,compareSumOrAv,params,group1_values,group2_values,compareSumOrAv,rectifyForStats);
    LFPfigs=[LFPfigs; l];
end