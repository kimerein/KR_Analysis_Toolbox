There are a few major analysis packages included in this folder:
1. Automatic experiment analysis -- look at raw LFP data, gamma, beta and alpha frequencies as a function of time for different stimulus and LED conditions; also consider unit activity as response over time and function of different stimulus and LED conditions 
Wrapper function is allAnalyses.m
2. Compare response windows between conditions
Main function is compareResponseWindows_wrapper.m
3. Histograms of response distributions by stim. condition, ideal observer analysis of trial-by-trial information in LFP, correlations within a response
Main functions are calculateStimCondResponseDistributions.m and correlateResponsePeriods.m

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FUNCTIONS IN PACKAGE 1: Automatic analysis of an experiment
alignToPhotodiode.m
allAnalyses.m
analyzeLFPbyStimCond.m
analyzeResponseToStim.m
anaylzeUnit.m
bandPassLFP.m
compareLFPconds.m (only compares 2 conditions, re-write to be more general)
Code from Ed's rotation for spectrograms
getAnalysisParams.m
getLEDConditions.m
getLFPbySweep.m

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FUNCTIONS IN PACKAGE 2: Compare responses (with specific windows) between conditions
compareResponseWindows.m
compareResponseWindows_wrapper.m

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
FUNCTIONS IN PACKAGE 3: Trial-by-trial response variability and correlations within a response
calculateStimCondResponseDistributions.m
calculateStimCondStats.m
calculateWindowStdev.m
correlateResponsePeriods.m
getBinomial_pValue.m

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
OTHER

FUNCTIONS FOR PRE-PROCESSING DATA:
alignLFPToPhotodiode.m
alignToPhotodiode.m
bandPassLFP.m
getLFPbySweep.m
calculateLEDartifact.m
calculateLEDartifact_wrapper.m
downSamp.m
getLEDConditions.m
getLFPbySweepSeveralFiles.m

GRAPHING FUNCTIONS:
axesmatrix.m
histnorm.m
K_psth.m
K_superPSTH.m

SPECTROGRAM CODE AND OTHER
CODE FROM ROTATION IN ED CALLAWAY'S LAB:
compare_freq_bands.m
create_gabormorlet.m
defineBrainStates.m
gabor_morlet_config.m
gabor_morlet_plot.m

FUNCTIONS FOR ANALYZING UNITS 
(EXCLUDING THOSE INCLUDED IN THE AUTOMATIC ANALYSIS PACKAGE):

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
OBSOLETE

OBSOLETE FUNCTIONS FOR ANALYZING LFP
(HAVE BEEN RE-WRITTEN OR INCORPORATED INTO AUTOMATIC EXPERIMENT ANALYSIS PACKAGE):
analyzeLFPduringStim.m
analyzeLFPSweeps.m
analyzeLFPSweepsForRinging.m
compareLFPblocks.m
getLFPAnalysisParams.m

OBSOLETE FUNCTIONS FOR ANALYZING UNITS:
KER_analyzeCell.m
KER_analyzeCell_backup.m

OBSOLETE PRE-PROCESSING FUNCTIONS:
calculateLEDartifactSimple.m


