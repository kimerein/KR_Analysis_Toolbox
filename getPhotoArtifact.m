function avLFP=getPhotoArtifact(expt,fileInd,LFPbySweep)

Fs=32000;
useLedCond=[2 4 6 8];
useWindowA=[1.29 1.5];
useWindowB=[3 3.5];
downSampFactor=5;
physCh=[23 15 22 14 19 11 20 12 17 9 16 8 18 10 21 13];
dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData\';
specgramFromThisChInd=3;
useTheseStimConds=9;

daqFileNames=expt.files.names(fileInd);
disp('Using these daq files');
disp(daqFileNames);

if isempty(LFPbySweep)
    % Get LFP
    [LFPbySweep,Fs,allCompleteSweeps]=getJustLFPbySweep(dataDir,daqFileNames,Fs,downSampFactor,physCh(specgramFromThisChInd));
    if any(size(LFPbySweep{1})==0)
        disp('No data in this daq file.');
    end
else
    Fs=Fs/downSampFactor;
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

% Get right LFP sweeps
temp=LFPbySweep{1};
LFPbySweep{1}=temp(ismember(ledConds,useLedCond) & ismember(stimConds,useTheseStimConds),:);

% Get average LFP
avLFP=mean(LFPbySweep{1},1);

% Zero everything but useWindows
x=0:1/Fs:(length(avLFP)-1)*(1/Fs);
avLFP(x<=useWindowA(1))=0;
avLFP(x>=useWindowA(2) & x<=useWindowB(1))=0;
avLFP(x>=useWindowB(2))=0;
avLFP(x>=useWindowA(1) & x<=useWindowA(2) & avLFP<=0)=0;
avLFP(x>=useWindowB(1) & x<=useWindowB(2) & avLFP>=0)=0;

% Plot photo-artifact
figure(); 
plot(x,avLFP);