function [UPstates,powerRatio,mua_xpoints,mua_ypoints,ledToReturn,LFPbySweep]=find_LFP_UPstates_wrapper(expt,fileInd,spikes,LFPbySweep)

powerRatio=[];
mua_xpoints=[];
mua_ypoints=[];
ledToReturn=[];
UPstates=[];


% lowerBand=30:50;
% upperBand=50:80;
% lowerBand=1:8;
% upperBand=20:50;
lowerBand=[5 30];
upperBand=[30 100];

Fs=32000;
downSampFactor=10;
noLEDcond=nan;
physCh=[23 15 22 14 19 11 20 12 17 9 16 8 18 10 21 13];
usePhysCh_ind=5;
% usePhysCh_ind=6;
dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData_older\';
% dataDir='W:\New Acquisition Computer\';
% useTheseStimConds=1;

daqFileNames=expt.files.names(fileInd);
disp('Using these daq files');
disp(daqFileNames); 

% Get LFP
% LFPbySweep=getJustLFPbySweep(dataDir,daqFileNames,Fs,downSampFactor,physCh);
if isempty(LFPbySweep)
    LFPbySweep=getJustLFPbySweep(dataDir,daqFileNames,Fs,downSampFactor,physCh(usePhysCh_ind));
end
% Line above me is right 
Fs=Fs/downSampFactor;
if any(size(LFPbySweep{1})==0)
    disp('No data in this daq file.');
end
% LFPbySweep=[]; % Only for case when not using LFP to define UP states

% Match acquired LFP sweeps to LED conditions
ledConds=expt.sweeps.led(ismember(expt.sweeps.fileInd,fileInd));
if length(ledConds)~=size(LFPbySweep{1},1)
    disp('ledConds and acquired LFP sweeps do not match');
    return
end

% Get stimconds
stimConds=expt.sweeps.stimcond(ismember(expt.sweeps.fileInd,fileInd));

% temp=LFPbySweep{usePhysCh_ind};
temp=LFPbySweep{1}; % Uncomment me if using LFP to define UP states
subLFPbySweep=LFPbySweep;
if ~isnan(noLEDcond)
    subLFPbySweep=temp(ismember(ledConds,noLEDcond),:);
    subSpikes=filtspikes(spikes,0,'led',noLEDcond);
    ledToReturn=ledConds(ismember(ledConds,noLEDcond));
else
    subLFPbySweep=temp;
    subSpikes=spikes;
    ledToReturn=nan;
end
disp('These numbers should be the same');
disp(size(subLFPbySweep,1));
% disp(length(unique(subSpikes.trials)));
disp(length(unique(subSpikes.sweeps.trials)));

[UPstates,powerRatio,mua_xpoints,mua_ypoints]=find_UPstates_with_LFP(subLFPbySweep,Fs,lowerBand,upperBand,subSpikes);
% [UPstates,mua_xpoints,mua_ypoints]=find_UPstates_withoutLFP(subSpikes);