function [avLFPbySweep,LFPbySweep,ledConds,stimConds]=measureLFPfollowing(expt,fileInd,LFPbySweep,ledConds,stimConds)

Fs=32000;
downSampFactor=10;
noLEDcond=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60 1.03 2.03 4.03 6.03 8.03 10.03 12.03 14.03 16.03 18.03 20.03 30.03 40.03 50.03 60.03];
useLEDcond=16;
physChInd=4;
physCh=[23 15 22 14 19 11 20 12 17 9 16 8 18 10 21 13];
dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData\';
useTheseStimConds=1:9;

if isempty(LFPbySweep)    
    daqFileNames=expt.files.names(fileInd);
    disp('Using these daq files');
    disp(daqFileNames);
    
    % Make output directory for CSD data
    % saveToDir='';
    % if ~exist(saveToDir,'dir')
    %     mkdir(saveToDir);
    % end
    
    % Get LFP
    [LFPbySweep,Fs,allCompleteSweeps]=getJustLFPbySweep(dataDir,daqFileNames,Fs,downSampFactor,physCh(physChInd));
    if any(size(LFPbySweep{1})==0)
        disp('No data in this daq file.');
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
    
    % Get LFP sweeps without LED on
    for i=1:length(LFPbySweep)
        temp=LFPbySweep{i};
        %     LFPbySweep{i}=temp(ledConds==noLEDcond & ismember(stimConds,useTheseStimConds),:);
        %     LFPbySweep{i}=temp(ledConds==noLEDcond,:);
        LFPbySweep{i}=temp(ismember(ledConds,noLEDcond) & ismember(stimConds,useTheseStimConds),:);
    end
    
    % Get average LFP
    for i=1:length(LFPbySweep)
        avLFPbySweep{i}=mean(LFPbySweep{i},1);
        % Shift
        avLFPbySweep{i}=avLFPbySweep{i}-mean(avLFPbySweep{i});
    end
    
    % Display average LFP
    figure();
    cumOffset=0;
    for i=1:length(avLFPbySweep)
        o=max(avLFPbySweep{i});
        cumOffset=cumOffset+o;
        plot(0:1/Fs:(length(avLFPbySweep{i})-1)*(1/Fs),avLFPbySweep{i}-cumOffset);
        hold on;
    end
end

% Fs=Fs/downSampFactor;
% % a=LFPbySweep{1};
% % a=a{1};
% a=LFPbySweep;
% subLFPbySweep=a(ismember(ledConds,useLEDcond),:);
% figure(); 
% plot(0:1/Fs:(size(subLFPbySweep,2)-1)*(1/Fs),mean(subLFPbySweep,1));

