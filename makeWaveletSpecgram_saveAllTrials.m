% function [LFPbySweep,Fs]=makeWaveletSpecgram_saveAllTrials(expt,fileInd,LFPbySweep,pa)
function [LFPbySweep,Fs,p,LFPspecgram,trialSpecgrams]=makeWaveletSpecgram_saveAllTrials(expt,fileInd,LFPbySweep,pa)

Fs=32000;
downSampFactor=5;
% noLEDcond=[3 5];
physCh=[23 15 22 14 19 11 20 12 17 9 16 8 18 10 21 13];
dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData\';
% useTheseStimConds=9;
specgramFromThisChInd=3;
removePhotoartifact=1;
photoLEDs=[2 4 6 8];

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
    
if removePhotoartifact==1
    c=0.45;
    temp=LFPbySweep{1};
    temp(ismember(ledConds,photoLEDs),:)=temp(ismember(ledConds,photoLEDs),:)-repmat(c.*pa,sum(ismember(ledConds,photoLEDs)),1);
    LFPbySweep{1}=temp;
end       

    % Get LFP sweeps without LED on
%     for i=1:length(LFPbySweep)
%         temp=LFPbySweep{i};
%         %     LFPbySweep{i}=temp(ledConds==noLEDcond & ismember(stimConds,useTheseStimConds),:);
%         %     LFPbySweep{i}=temp(ledConds==noLEDcond,:);
%         LFPbySweep{i}=temp(ismember(ledConds,noLEDcond) & ismember(stimConds,useTheseStimConds),:);
%     end


% Get average LFP
% for i=1:length(LFPbySweep)
%     avLFPbySweep{i}=mean(LFPbySweep{i},1);
%     % Shift
%     avLFPbySweep{i}=avLFPbySweep{i}-mean(avLFPbySweep{i});
% end

% Display average LFP
% figure();
% cumOffset=0;
% for i=1:length(avLFPbySweep)
%     o=max(avLFPbySweep{i});
%     cumOffset=cumOffset+o;
%     plot(0:1/Fs:(length(avLFPbySweep{i})-1)*(1/Fs),avLFPbySweep{i}-cumOffset);
%     hold on;
% end

[p,LFPspecgram,allTrialSpecgrams,otherDetails]=makeWaveletSpecgram(LFPbySweep{1},1:size(LFPbySweep{1},1),Fs);
trialSpecgrams.specgrams=allTrialSpecgrams;
trialSpecgrams.ledConds=ledConds;
trialSpecgrams.stimConds=stimConds;
trialSpecgrams.otherDetails=otherDetails;

for i=1:floor(length(trialSpecgrams.specgrams)/12)
    if (i-1)*12+12>length(trialSpecgrams.specgrams)
        a=trialSpecgrams.specgrams((i-1)*12+1:end);
        trialDetails.ledConds=trialSpecgrams.ledConds((i-1)*12+1:end);
        trialDetails.stimConds=trialSpecgrams.stimConds((i-1)*12+1:end);
    else
        a=trialSpecgrams.specgrams((i-1)*12+1:(i-1)*12+12);
        trialDetails.ledConds=trialSpecgrams.ledConds((i-1)*12+1:(i-1)*12+12);
        trialDetails.stimConds=trialSpecgrams.stimConds((i-1)*12+1:(i-1)*12+12);
    end
    trialDetails.otherDetails=trialSpecgrams.otherDetails;
    save(['Z:\Updating Data\Analysis Computer\All Specgram Data\Mawake54B\minus photoartifact\fileAIs_130to133_specgrams_' num2str(i) '.mat'],'a');
    save(['Z:\Updating Data\Analysis Computer\All Specgram Data\Mawake54B\minus photoartifact\fileAIs_130to133_details_' num2str(i) '.mat'],'trialDetails');
end