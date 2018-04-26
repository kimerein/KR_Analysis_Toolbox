function ledTraces=checkLEDonsets(expt,fileInd,dataDir)

Fs=32000;
downSampFactor=3;
ledCh=2;
% dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData_older\';
% dataDir='W:\New Acquisition Computer\';

daqFileNames=expt.files.names(fileInd);
disp('Using these daq files');
disp(daqFileNames);

% Get LED input
[ledbySweep,Fs]=getJustLFPbySweep(dataDir,daqFileNames,Fs,downSampFactor,ledCh);
if any(size(ledbySweep{1})==0)
    disp('No data in this daq file.');
end

figure(); 
plot(0:1/Fs:(length(ledbySweep{1})-1)*(1/Fs),ledbySweep{1});
ledTraces.xpoints=0:1/Fs:(length(ledbySweep{1})-1)*(1/Fs);
ledTraces.ypoints=ledbySweep{1};