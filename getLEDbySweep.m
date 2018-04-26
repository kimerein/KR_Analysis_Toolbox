function [Fs,ledAv,ledData,gotDaqs,allCompleteSweeps,fromFileInds]=getLEDbySweep(dataDir,daqFileNames,Fs,downSampFactor,ledChannel)
% Extract data from daq files and organize by trials
% Will not include data from a new daq file if the length of the sweep does
% not match sweep length from the previously read daq file
%
% RETURNS:
% Fs: sampling rate of all returned data
% ledAv: the average over all trials of the LED signal from daq files
% ledData: a matrix of raw LED data, columns are samples at different
% times, rows are trials
% gotDaqs: a vector with length equal to the number of daq file names
% passed in; has 1 for each daq file included in the returned data, 0 for
% each daq file excluded from the returned data
% allCompleteSweeps: vector the length of the rows of LFPbySweep, 1 for
% each trial that was complete, 0 for each trial that had to be padded
% (0's) because it was shorter than a full trial
% 
% PARAMETERS:
% daqFileNames: a cell array of daq file names to be read and processed
% Fs: sampling rate of the data
% downSampFactor: factor by which to down-sample all data
% ledChannel: the Matlab index for LED signal


ledData=[];
ledAv=[];
allCompleteSweeps=[];
gotDaqs=zeros(length(daqFileNames),1);
fromFileInds=[];
for i=1:length(daqFileNames)
    daqFileName=daqFileNames{i};
    [Fs2,thisLEDData,completeSweeps]=readTrialTriggeredLED(i,dataDir,daqFileName,Fs,downSampFactor,ledChannel);
    if any(size(thisLEDData)==0)
        disp(['No data from daq file ' daqFileName '. Not using it.']);
        gotDaqs(i)=0;
%     elseif size(ledData,2)>0 && (size(thisLEDData,2)~=size(ledData,2))
%         disp(['Sweep length from ' daqFileName ' does not match other sweeps. Not using data from this file.']);
%         gotDaqs(i)=0;
    % Use below if calling getLEDbySweep from MakeExptSRO.m
    elseif size(ledData,2)>0 && (size(thisLEDData,2)~=size(ledData,2))
        if size(ledData,2)>size(thisLEDData,2)
            % Pad thisLEDData with zeros
            gotDaqs(i)=1;
            ledData=[ledData; [thisLEDData zeros(size(thisLEDData,1),size(ledData,2)-size(thisLEDData,2))]];
            allCompleteSweeps=[allCompleteSweeps; completeSweeps];
            fromFileInds=[fromFileInds; i];
            % Don't bother with LED average
        elseif size(thisLEDData,2)>size(ledData,2)
            % Pad ledData with zeros
            gotDaqs(i)=1;
            ledData=[ledData zeros(size(ledData,1),size(thisLEDData,2)-size(ledData,2))];
            ledData=[ledData; thisLEDData];
            allCompleteSweeps=[allCompleteSweeps; completeSweeps];
            fromFileInds=[fromFileInds; i];
            % Don't bother with LED average
        end
    else
        gotDaqs(i)=1;
        ledData=[ledData; thisLEDData];
        allCompleteSweeps=[allCompleteSweeps; completeSweeps];
        fromFileInds=[fromFileInds; i];
%         if i==1
%             ledAv=mean(thisLEDData((completeSweeps==1)&logical([0; ones(size(thisLEDData,1)-1,1)]),:),1);
%         else
%             ledAv=mean([ledAv; mean(thisLEDData((completeSweeps==1)&logical([0; ones(size(thisLEDData,1)-1,1)]),:),1)],1);
%         end
    end
end
Fs=Fs2;

gotDaqs=logical(gotDaqs);
end