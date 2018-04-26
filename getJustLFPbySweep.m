function [LFPbySweep,Fs,allCompleteSweeps]=getJustLFPbySweep(dataDir,daqFileNames,Fs,downSampFactor,physChannels)
% Extract data from daq files and organize by trials
% Will not include data from a new daq file if the length of the sweep does
% not match sweep length from the previously read daq file
%
% RETURNS:
% LFPbySweep: matrix of raw physiology data from daq files, columns are samples at
% different times, rows are trials (separated by NaNs in original .daq
% data), each row has been band-passed between def_LowerCutoff and
% def_UpperCutoff
% Fs: sampling rate of all returned data
% bandPassedLFPbySweep: each row of LFPbySweep is band-passed between lower_cutoff and 
% upper_cutoff, same set-up as LFPbySweep 
% photoData: matrix of raw photodiode data from daq files, columns are
% samples at different times, rows are trials
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
% def_LowerCutoff: lower cutoff for band-passing LFPbySweep
% def_UpperCutoff: upper cutoff for band-passing LFPbySweep
% lower_cutoff: lower cutoff for band-passing bandPassedLFPbySweep
% upper_cutoff: upper cutoff for band-passing bandPassedLFPbySweep
% downSampFactor: factor by which to down-sample all data
% physChannel: the Matlab index for the physiology channel on the daq
% board (that is, when daq file was written, a channel number was definedf or the physiology channel,
% this number is what is referenced by phyChannel)
% photoChannel: same as physChannel for photodiode signal
% ledChannel: same as physChannel for LED signal

for i=1:length(physChannels)
    LFPbySweep{i}=[];
end

allCompleteSweeps=[];
gotDaqs=zeros(length(daqFileNames),1);
for i=1:length(daqFileNames)
    daqFileName=daqFileNames{i};
    [thisLFPbySweep,Fs2,completeSweeps]=readJustStimTriggeredLFP(dataDir,daqFileName,Fs,downSampFactor,physChannels);
    if i==1
        Fs3=Fs2;
    end
    if any(size(thisLFPbySweep{1})==0)
        disp(['No data from daq file ' daqFileName '. Not using it.']);
        gotDaqs(i)=0;
    elseif size(LFPbySweep{1},2)>0 && (size(thisLFPbySweep{1},2)~=size(LFPbySweep{1},2))
        disp(['Sweep length from ' daqFileName ' does not match other sweeps. Not using data from this file.']);
        gotDaqs(i)=0;
    else
        gotDaqs(i)=1;
        for j=1:length(physChannels)
            LFPbySweep{j}=[LFPbySweep{j}; thisLFPbySweep{j}];
        end
        allCompleteSweeps=[allCompleteSweeps; completeSweeps];
    end
end
Fs=Fs3;

gotDaqs=logical(gotDaqs);
end