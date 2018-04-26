function [dataInt,avLFPbySweep,LFPbySweep,ledConds,stimConds]=makeCSD_sarah(expt,fileInd)
% 17-20 fileInd

Fs=32000;
downSampFactor=10;
noLEDcond=[1 2 3];
physCh=[23 15 22 14 19 11 20 12 17 9 16 8 18 10 21 13]; % superficial to deep
dataDir='Z:\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData\';
% useTheseStimConds=1;

daqFileNames=expt.files.names(fileInd);
disp('Using these daq files');
disp(daqFileNames);

% Make output directory for CSD data
% saveToDir='';
% if ~exist(saveToDir,'dir')
%     mkdir(saveToDir);
% end

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

% Get LFP sweeps without LED on
for i=1:length(LFPbySweep)
    temp=LFPbySweep{i};
%     LFPbySweep{i}=temp(ledConds==noLEDcond & ismember(stimConds,useTheseStimConds),:);
%     LFPbySweep{i}=temp(ledConds==noLEDcond,:);
    LFPbySweep{i}=temp(ismember(ledConds,noLEDcond),:);
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
