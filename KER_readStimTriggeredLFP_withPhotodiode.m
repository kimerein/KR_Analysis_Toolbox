function [LFPbySweep,Fs,photoData]=KER_readStimTriggeredLFP_withPhotodiode(daqFileName,Fs,LPcutoff,HPcutoff,physChannel,photoChannel)
% Reads the daq file containing the LFP from daqFileName
% Plots the LFP averaged across all sweeps
% physChannel is the physiology channel of the rig (NI-DAQ channel number
% in Matlab)
% photoChannel is the photodiode channel of the rig (NI-DAQ channel number
% in Matlab)
% Fs is the sampling frequency of the data in daqFileName
% Returns an array of low-pass filtered (below LPcutoff) and high-pass filtered
% (above HPcutoff) data and the new sampling frequency (newFs) of LFPbySweep
% Each row of LFPbySweep is the data from one sweep
% Different rows are different sweeps
% Different columns are different samples

data=daqread(daqFileName);
physData=data(:,physChannel);
photoData=data(:,photoChannel);

sweepSeps=find(isnan(physData));
if numel(sweepSeps)~=0
    sweepSepData=zeros(length(sweepSeps)+1,sweepSeps(1)-1);
    for i=1:length(sweepSeps)-1
        if i==100
            disp('processing\n');
        end
        sweepSepData(i,:)=physData(sweepSeps(i)+1:sweepSeps(i+1)-1)';
    end
else
    sweepSepData=physData;
end

DownSampFactor = 10;
bDownSamp = 1;
if bDownSamp
    [filtSweepSepData Fs]=KER_DownSamp(sweepSepData,Fs,DownSampFactor);
end

%photoData
photoData=photoData(1:DownSampFactor:end);
%photoData

% Low-pass-filter data
disp('LP-filtering');
filtSweepSepData=fftFilter(filtSweepSepData',Fs,LPcutoff,1);
filtSweepSepData=filtSweepSepData';
disp('Done LP-filtering');

% High-pass-filter data
disp('HP-filtering');
filtSweepSepData=fftFilter(filtSweepSepData',Fs,HPcutoff,2);
filtSweepSepData=filtSweepSepData';
disp('Done HP-filtering');

trialLength=size(filtSweepSepData,2)*(1/Fs);
plot(0:trialLength/(size(filtSweepSepData,2)-1):trialLength,mean(filtSweepSepData,1));
title('Average Stimulus-Triggered LFP');

LFPbySweep=filtSweepSepData;
