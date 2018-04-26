function [runningTrials,velocity,runningData,ledData,Fs]=getRunningTrials(expt, fileInd)

% dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData_older\';
dataDir='W:\New Acquisition Computer\';

Fs=32000;
runningCh=5;
ledCh=4;
downSampFactor=10;
scaleFactor=0.9392;
discontinuityThresh=0.2;
fractionRunningThresh=0.25;

daqFileNames=expt.files.names(fileInd);

% Get running data
[runningData,ledData,Fs]=getRunningDataBySweep(dataDir,daqFileNames,Fs,downSampFactor,runningCh,ledCh);
if any(size(runningData)==0)
    disp('No data in these daq files');
end

if length(expt.sweeps.trials(ismember(expt.sweeps.fileInd,fileInd)))~=size(runningData,1)
    disp('runningData and number of trials in expt do not match');
    return
end

% Determine animal's velocity
r=runningData-scaleFactor*ledData;
for i=1:size(r,1)
    r(i,1:10)=r(i,11)*ones(1,10);
end
% Remove discontinuities in position data
rdiff=diff(r');
rdiff=rdiff';
for i=1:size(rdiff,1)
    inds=find(abs(rdiff(i,:))>=discontinuityThresh);
    useinds=inds(inds>10);
    for j=useinds
        if j+4>size(r,2)
            r(i,end-4:end)=r(i,end-4)*ones(1,5);
        else
            r(i,:)=[r(i,1:j) r(i,j) r(i,j) r(i,j) r(i,j+4:size(r,2))-r(i,j+4)+r(i,j)];
        end
    end
end
der=r(:,1:end-1000)-r(:,1000+1:end); % 1000 fewer points than r per row
for i=1:size(der,1)
    velocity(i,:)=[der(i,:) der(i,end)*ones(1,1000)];
end
velocity=abs(velocity);
for i=1:size(velocity,1)
    velocity(i,:)=smooth(velocity(i,:)); % default window is 5
end

% Define a "running trial" as one for which more than half the trial has
% velocity > 0
runningTrials=zeros(size(velocity,1),1);
for i=1:size(velocity,1)
    if sum(velocity(i,:)>5*10^-3)>fractionRunningThresh*size(velocity,2)
        runningTrials(i)=1;
    end
end
end


function [dataBySweep,ledData,Fs]=getRunningDataBySweep(dataDir, daqFileNames, Fs, downSampFactor, runningCh, ledCh)
dataBySweep=[];
ledData=[];
for i=1:length(daqFileNames)
    daqFileName=daqFileNames{i};
    [thisDatabySweep,thisledDataBySweep,Fs2]=getRunningData(daqFileName,dataDir,runningCh,downSampFactor,Fs,ledCh);
    Fs=Fs2;
    if any(size(thisDatabySweep)==0)
        disp(['No data from daq file ' daqFileName '. Not using it.']);
    elseif size(dataBySweep,2)>0 && (size(thisDatabySweep,2)~=size(dataBySweep,2))
        disp(['Sweep length from ' daqFileName ' does not match other sweeps. Not using data from this file.']);
    else
        dataBySweep=[dataBySweep; thisDatabySweep];
        ledData=[ledData; thisledDataBySweep];
    end   
end
end

function [sepData,ledData,Fs]=getRunningData(daqFileName, dataDir, runningCh, downSampFactor, Fs, ledCh)
   sdata=daqread(strcat(dataDir,daqFileName),'Channels',[ledCh runningCh]);
   data=sdata(:,2);
   lData=sdata(:,1);
   sweepSeps=find(isnan(data));
   disp('Loading DAQ encoder data');
   if numel(sweepSeps)~=0
       sepData=zeros(length(sweepSeps)+1,sweepSeps(1)-1);
       lsepData=zeros(length(sweepSeps)+1,sweepSeps(1)-1);
       sepData(1,:)=data(1:sweepSeps(1)-1)';
       lsepData(1,:)=lData(1:sweepSeps(1)-1)';
       for i=1:length(sweepSeps)
           if mod(i,100)==0
               disp('Loading DAQ encoder data');
           end
           if i==length(sweepSeps)
               if length(data)-sweepSeps(end)==length(sepData(1,:))
                   sepData(i+1,:)=data(sweepSeps(i)+1:end)';
                   lsepData(i+1,:)=lData(sweepSeps(i)+1:end)';
               else
                   sepData(i+1,:)=[data(sweepSeps(i)+1:end); zeros(length(sepData(1,:))-(length(data)-sweepSeps(end)),1)]';
                   lsepData(i+1,:)=[lData(sweepSeps(i)+1:end); zeros(length(lsepData(1,:))-(length(lData)-sweepSeps(end)),1)]';
               end
           else
               sepData(i+1,:)=data(sweepSeps(i)+1:sweepSeps(i+1)-1)';
               lsepData(i+1,:)=lData(sweepSeps(i)+1:sweepSeps(i+1)-1)';
           end
       end
   else
       sepData=data';
       lsepData=data';
   end
   if downSampFactor~=1
       currFs=Fs;
       [filtSweepSepData Fs]=downSamp(sepData,currFs,downSampFactor);
       [lfiltSweepSepData Fs2]=downSamp(lsepData,currFs,downSampFactor);
   else
       filtSweepSepData=sepData;
       lfiltSweepSepData=lsepData;
   end
   sepData=filtSweepSepData;
   ledData=lfiltSweepSepData;
end