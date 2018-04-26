function [dataInt,avLFPbySweep,LFPbySweep,ledConds,stimConds]=makeCSDfig(expt,fileInd,noLEDcond,useTheseStimConds,dataDir)
% function [dataInt,avLFPbySweep,LFPbySweep,ledConds,stimConds]=makeCSDfig(expt,fileInd,noLEDcond,useTheseStimConds,dataDir)
% INPUT
% expt - expt structure from exv
% fileInd - the indices of daq files to use for CSD
% noLEDcond - the LED conditions to select trials for CSD (vector of
% values, e.g., [0 5])
% useTheseStimConds - the stimulus condition (stimcond) condition to select
% trials (vector of values, e.g., [1 2 3 4 5 6 7 8])
% dataDir - string specifying the directory containing .daq files
% 
% OUTPUT
% dataInt - a matrix containing the CSD, view with imagesc, for example
% avLFPbySweep - the average LFP
% LFPbySweep - a cell array of n dimensions where n is the number of
% channels, each element containing a matrix of LFP data where rows are
% trials and columns are samples
% ledConds - a vector of led conditions corresponding to each trial in
% LFPbySweep
% stimConds - a vector of led conditions corresponding to each trial in
% LFPbySweep

% Order of channels from RigDefs -- each user should customize to rig
physCh=[23 15 22 14 19 11 20 12 17 9 16 8 18 10 21 13];

Fs=32000;
downSampFactor=10;
dataDir=[dataDir '\'];

dataInt=[];

daqFileNames=expt.files.names(fileInd);
disp('Using these daq files');
disp(daqFileNames);

% Get LFP
[LFPbySweep,Fs,allCompleteSweeps]=getJustLFPbySweep(dataDir,daqFileNames,Fs,downSampFactor,physCh);
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
subplot(2,1,1);
cumOffset=0;
for i=1:length(avLFPbySweep)
    o=max(avLFPbySweep{i});
    cumOffset=cumOffset+o;
    plot(0:1/Fs:(length(avLFPbySweep{i})-1)*(1/Fs),avLFPbySweep{i}-cumOffset);
    hold on;
end

% Calculate CSD
csd=zeros(length(avLFPbySweep),size(avLFPbySweep{1},2));
csd(1,:)=(2*avLFPbySweep{1}+avLFPbySweep{2})/3;
csd(16,:)=(2*avLFPbySweep{16}+avLFPbySweep{15})/3;
for i=2:15
    csd(i,:)=(avLFPbySweep{i-1}+2*avLFPbySweep{i}+avLFPbySweep{i+1})/4;
end
for i=1:16
    csd(i,:)=smooth(csd(i,:),17);
end
data1=csd(1:2:16,:);
data2=csd(2:2:16,:);
data1=diff(data1,2,1);
data2=diff(data2,2,1);
data=zeros(length(data1),2*size(data1,1));
DataCounter=1;
for i=1:size(data1,1)
    data(:,DataCounter)=data1(i,:);
    DataCounter=DataCounter+1;
    data(:,DataCounter)=data2(i,:);
    DataCounter=DataCounter+1;
end
data=-data';

% 2D linear interpolation;
[x y]=size(data);
y=1:y;
x=1:x;
[xi yi]=meshgrid(1:0.1:max(x),1:0.05:max(y));
dataInt=interp2(x,y,data',xi,yi);

subplot(2,1,2);
dataInt=dataInt';
imagesc(dataInt);
drawnow;
