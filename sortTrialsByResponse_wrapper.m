function trialClassifications=sortTrialsByResponse_wrapper(LFPbySweep,bandPassedLFPbySweep,ledForSweeps,LFP_Fs,baselineWindow,earlyWindow,earlyThresh,lateWindow,lateThresh,compareGroup1,compareGroup2,excludeFirstAndLastTrials)

trialClassifications=sortSingleTrials(bandPassedLFPbySweep,ledForSweeps,LFP_Fs,baselineWindow,earlyWindow,earlyThresh,lateWindow,lateThresh);
disp(['Fraction of total trials that are blue (in compareGroup1): ' num2str(sum(ismember(trialClassifications,compareGroup1)))]);
disp(['Fraction of total trials that are red (in compareGroup2): ' num2str(sum(ismember(trialClassifications,compareGroup2)))]);

if excludeFirstAndLastTrials
    LFPbySweep=LFPbySweep(2:end-1,:);
    bandPassedLFPbySweep=bandPassedLFPbySweep(2:end-1,:);
    ledForSweeps=ledForSweeps(2:end-1);
    trialClass2=trialClassifications(2:end-1);
end
plotBestTrials_general(LFPbySweep,bandPassedLFPbySweep,LFP_Fs,trialClass2,baselineWindow,compareGroup1,compareGroup2);
