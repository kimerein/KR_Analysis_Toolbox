function [dataInt,avLFPbySweep,LFPbySweep,ledConds,stimConds]=makeCSD2(expt,fileInd,avlfp,lfpbs,noLEDcond)
% 17-20 fileInd

% Fs=32000;
Fs=25000;
downSampFactor=1;
noLEDcond=[1 1.05 2 2.05 4 4.05 6 6.05 8 8.05 10 10.05 12 12.05 14 14.05 16 16.05 18 18.05 20 20.05 30 30.05 40 40.05 50 50.05 60 60.05];
% noLEDcond=[noLEDcond 1.01 2.01 4.01 6.01 8.01 10.01 12.01 14.01 16.01 18.01 20.01 30.01 40.01 50.01 60.01];
% noLEDcond=[noLEDcond 1.03 2.03 4.03 6.03 8.03 10.03 12.03 14.03 16.03 18.03 20.03 30.03 40.03 50.03 60.03];
% noLEDcond=[noLEDcond 0.0100    0.0200    0.0400    0.0600    0.0800    0.1000    0.1200    0.1400 0.1600    0.1800    0.2000    0.3000    0.4000    0.5000    0.6000    5.0100 5.0200    5.0400    5.0600    5.0800    5.1000    5.1200    5.1400    5.1600 5.1800    5.2000    5.3000    5.4000    5.5000    5.6000];
% noLEDcond=[0 0.1515 0.5050 0.7575 1.010];
% noLEDcond=[0];
% noLEDcond=[ 1.05 2.05 4.05 6.05 8.05 10.05 12.05 14.05  16.05 18.05 20.05 30.05 40.05  50.05  60.05];
% noLEDcond=[0 0.05 5.00 5.05];
% noLEDcond=[0 0.1515 0.2525];
% noLEDcond=[0 0.05 0.15 0.25 0.35 0.55 2.05 5.00 5.05];
% noLEDcond=[0 0.7575 0.03 0.05 0.1515 0.2525 0.5050 0.15 0.25 0.35 0.55 2.05 5.00 5.05 1 1.05 2 2.05 4 4.05 6 6.05 8 8.05 10 10.05 12 12.05 14 14.05 16 16.05 18 18.05 20 20.05 30 30.05 40 40.05 50 50.05 60 60.05];
% noLEDcond=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
% physCh=[23 15 22 14 19 11 20 12 17 9 16 8 18 10 21 13];
% physCh=[25 17]; % dLGN, V1; bottom set of chs, top set
% physCh=[32]; % top of dLGN, nearest hippocampus
physCh=[31 24 30 5];
% physCh=[6]; 
% physCh=[23]; % more superficial channel in V1
% physCh=[23 32]; % both V1 and hippo

% physCh=physCh(6);
dataInt=[];
avLFPbySweep=[];

% dataDir='Z:\Updating Data\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData\';
% dataDir='E:\From new server PART 1\New Acquisition Computer\';
% dataDir='F:\From new server PART 2\New Acquisition Computer 2\';
dd=RigDefs();
dataDir=dd.Dir.Data;
% dataDir='F:\New RawData\';
% dataDir='C:\Users\Admin\Documents\MATLAB\RawData\';
useTheseStimConds=1:128;

daqFileNames=expt.files.names(fileInd);
disp('Using these daq files');
disp(daqFileNames);

% Make output directory for CSD data
% saveToDir='';
% if ~exist(saveToDir,'dir')
%     mkdir(saveToDir);
% end

if ~isempty(lfpbs)
    LFPbySweep=lfpbs;
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
elseif isempty(avlfp) && isempty(lfpbs)
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
else
    avLFPbySweep=avlfp;
    ledConds=[];
    stimConds=[];
end
return

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

% 2D linear interpolation
% [x y]=size(data);
% y=1:y;
% x=1:x;
% [xi yi]=meshgrid(1:1:max(x),1:0.05:max(y));
% dataInt=interp2(x,y,data',xi,yi);
[x y]=size(data);
y=1:y;
x=1:x;
[xi yi]=meshgrid(1:0.1:max(x),1:0.05:max(y));
dataInt=interp2(x,y,data',xi,yi);

subplot(2,1,2);
dataInt=dataInt';
imagesc(dataInt);
drawnow;

% Calculate CSD
% % csd=zeros(length(avLFPbySweep)-2,size(avLFPbySweep{1},2));
% % for i=2:length(avLFPbySweep)-1
% %     csd(i-1,:)=avLFPbySweep{i+1}-avLFPbySweep{i}+avLFPbySweep{i-1};
% % end
% % csd=zeros(length(avLFPbySweep)-4,size(avLFPbySweep{1},2));
% for i=3:length(avLFPbySweep)-2
%     csd(i-2,:)=avLFPbySweep{i+2}-avLFPbySweep{i}+avLFPbySweep{i-2};
% end
% 
% % Smoothing
% % for j=1:50
% %     for i = 1:size(csd,1)
% %         csd(i,:) = smooth(csd(i,:),17);
% %     end
% % end
%     
% % Plot CSD
% subplot(2,1,2);
% % heatmap(csd,'Standardize','column');
% heatmap(csd);

% subplot(3,1,3);
