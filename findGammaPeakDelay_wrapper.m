function [delays,LFPbySweep]=findGammaPeakDelay_wrapper(expt,fileInd,useLEDconds,useStimConds,timeWindow,freqBand,LFPbySweep,Fs)

downSampFactor=10;
physCh=[23 15 22 14 19 11 20 12 17 9 16 8 18 10 21 13];
trialDuration=5;
ledOnset=0.3;

physCh=physCh(5);
dataInt=[];
% dataDir='W:\New Acquisition Computer\';
dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData\';
daqFileNames=expt.files.names(fileInd);
disp('Using these daq files');
disp(daqFileNames);

if isempty(LFPbySweep)
    Fs=32000;
    [LFPbySweep,Fs]=getJustLFPbySweep(dataDir,daqFileNames,Fs,downSampFactor,physCh);
    if any(size(LFPbySweep{1})==0)
        disp('No data in this daq file.');
    end
end
% Match acquired LFP sweeps to LED conditions
ledConds=expt.sweeps.led(ismember(expt.sweeps.fileInd,fileInd));
if length(ledConds)~=size(LFPbySweep{1},1)
    disp('ledConds and acquired LFP sweeps do not match');
    return
end
% Get stimconds
stimConds=expt.sweeps.stimcond(ismember(expt.sweeps.fileInd,fileInd));
if length(stimConds)~=size(LFPbySweep{1},1)
    disp('stimConds and acquired LFP sweeps do not match');
    return
end
% Get trough delays
xpoints=linspace(0,trialDuration,size(LFPbySweep{1},2));
delays=[];
delays=findGammaPeakDelay(xpoints,LFPbySweep{1},ledConds,useLEDconds,stimConds,useStimConds,timeWindow,freqBand,Fs,ledOnset);

