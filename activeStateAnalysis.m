function [usingTrialsLED1,usingTrialsLED2,LFPbySweep]=activeStateAnalysis(expt,spikes,pr,UPstates_LFP,LFPbySweep,muadata)

% Look in each sub-function for specific conditions of this analysis!!
fileInd=[1:10000];
saveDir='W:\Analysis Computer\Active States Analysis for Figures\';
ledForGettingUPthresh=[1.01];
trialDuration=5;
muaBinSize=40; % in ms
ledCond1=[1.01 7.07];
ledCond2=[3.03 5.05];
% ledCond1=[0];
% ledCond2=[5];
allStimcond=1:32;
% LEDwindow=[1.2 1.325];
% LEDwindow=[3.32 3.82];
% LEDwindow=[0.32 0.82];
LEDwindow=[1.34 1.84];
% LEDwindow=[1 2];
% LEDwindow=[3.22 3.345];
% LEDwindow=[3.221 3.221+0.125];

%% Step 1: Read in LFP and calculate power ratio
% find_LFP_UPstates_wrapper.m
% Look in this sub-function for Fs, physCh, usePhysCh_ind, dataDir
if isempty(pr)
    if ~isempty(LFPbySweep)
        [UPs,pr,muax,muay,ledToReturn,LFPbySweep]=find_LFP_UPstates_wrapper(expt,fileInd,spikes,LFPbySweep);
    else
        [UPs,pr,muax,muay,ledToReturn,LFPbySweep]=find_LFP_UPstates_wrapper(expt,fileInd,spikes,[]);
        save([saveDir 'LFPbySweep.mat'],'LFPbySweep');
        save([saveDir 'powerRatio.mat'],'pr');
    end
end

%% Step 2: Find best active/inactive state cut-off threshold
% find_UPDOWN_cutoff.m
% Look in this sub-function for windowWidth,windowOffset,windowTimewindow
% [fitted,~,~,~,~,UPthresh]=find_UPDOWN_cutoff(spikes,ledForGettingUPthresh,[],pr,[],[],allStimcond);
% % UPthresh=9.2;
UPthresh=2;

%% Step 3: Get LFP power ratio-defined active states using this best threshold
if isempty(UPstates_LFP)
    x=linspace(0,trialDuration,length(pr{1}));
    UPstates_LFP=returnLFPUPs(x,pr,UPthresh);
    save([saveDir 'UPstates_LFP.mat'],'UPstates_LFP');
end

%% Step 4: Get MUA data for these active states
if isempty(muadata)
    [muax,muay,trialLEDs,trialStims]=get_muax_muay(expt,filtspikes(spikes,0,'fileInd',fileInd),fileInd,muaBinSize,trialDuration,0);
    save([saveDir 'muaData_forUPs.mat'],'muax','muay','trialLEDs','trialStims');
else
    muax=muadata.muax;
    muay=muadata.muay;
    trialLEDs=muadata.trialLEDs;
    trialStims=muadata.trialStims;
end

%% For awake only, use MUA rather than power ratio to define UP state
% Comment out for anesthetized
% UPthresh=108;
% UPstates_LFP=returnLFPUPs(muax,muay,UPthresh);
% save([saveDir 'UPstates_MUA.mat'],'UPstates_LFP');

%% Step 5: Show average active state
% alignUPs_general.m
% Look in this sub-function for beforeLED, takeBefore, takeAfter
% [~,~,xled1,yled1,~,~,usingTrialsLED1]=alignUPs_general(UPstates_LFP,trialLEDs,ledCond1,trialStims,allStimcond,LEDwindow,muax,muay);
% [~,~,xled2,yled2,~,~,usingTrialsLED2]=alignUPs_general(UPstates_LFP,trialLEDs,ledCond2,trialStims,allStimcond,LEDwindow,muax,muay);
% alignedToLED.xled1=xled1;
% alignedToLED.yled1=yled1;
% alignedToLED.xled2=xled2;
% alignedToLED.yled2=yled2;
% save([saveDir 'alignedToLED.mat'],'alignedToLED');
[xled1,yled1,~,~,~,~,usingTrialsLED1]=alignUPs_general(UPstates_LFP,trialLEDs,ledCond1,trialStims,allStimcond,LEDwindow,muax,muay);
[xled2,yled2,~,~,~,~,usingTrialsLED2]=alignUPs_general(UPstates_LFP,trialLEDs,ledCond2,trialStims,allStimcond,LEDwindow,muax,muay);
alignedToOnset.xled1=xled1;
alignedToOnset.yled1=yled1;
alignedToOnset.xled2=xled2;
alignedToOnset.yled2=yled2;
save([saveDir 'alignedToOnset.mat'],'alignedToOnset');
