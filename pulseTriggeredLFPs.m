% Show pulse-triggered LFP and filtered LFP

takeFractionOfTrials=0.05;
%takeEveryNTrials=40;

a=load('E:\Results\GammaLFPs\KR_2010-08-02\Stim60\PhotoAligned_Data\LFPbySweep.mat');
LFPbySweep=a.LFPbySweep;
a=load('E:\Results\GammaLFPs\KR_2010-08-02\Stim60\PhotoAligned_Data\ledForSweeps.mat');
ledForSweeps=a.ledForSweeps;
a=load('E:\Results\GammaLFPs\KR_2010-08-02\Stim60\PhotoAligned_Data\ledAv.mat');
ledAv=a.ledAv;
a=load('E:\Results\GammaLFPs\KR_2010-08-02\Stim60\PhotoAligned_Data\LFP_Fs.mat');
LFP_Fs=a.LFP_Fs;
a=load('E:\Results\GammaLFPs\KR_2010-08-02\Stim60\PhotoAligned_Data\bandPassedLFPbySweep.mat');
bandPassedLFPbySweep=a.bandPassedLFPbySweep;

% Align to initial values
%LFPbySweep=LFPbySweep-LFPbySweep(:,1)*ones(1,size(LFPbySweep,2));

% Align all trials to mean for each trial
LFPbySweep=LFPbySweep-mean(LFPbySweep,2)*ones(1,size(LFPbySweep,2));
bandPassedLFPbySweep=bandPassedLFPbySweep-mean(bandPassedLFPbySweep,2)*ones(1,size(bandPassedLFPbySweep,2));

% Normalize ledForSweeps to make it a logical array
% Assumes a step LED (that is, only 2 LED values)
minLED=min(ledForSweeps);
maxLED=max(ledForSweeps);
if minLED==maxLED
    ledForSweeps(ledForSweeps==0)=0;
    ledForSweeps(ledForSweeps~=0)=1;
else
    ledForSweeps(ledForSweeps==minLED)=0;
    ledForSweeps(ledForSweeps==maxLED)=1;
end
set1=logical(ledForSweeps);
set2=~logical(ledForSweeps);

% Only use LED ON trials
LFPbySweep=LFPbySweep(set1,:);
bandPassedLFPbySweep=bandPassedLFPbySweep(set1,:);

t=floor(size(LFPbySweep,1)*takeFractionOfTrials);
all=randperm(size(LFPbySweep,1));
ts=all(1:t);

cs=[];
figure;
subplot(3,1,1);
%for i=1:takeEveryNTrials:size(LFPbySweep,1)
for i=ts
    randColor=[rand(1,1) rand(1,1) rand(1,1)];
    cs=[cs; randColor];
    plot((1/LFP_Fs):(1/LFP_Fs):(1/LFP_Fs)*size(LFPbySweep,2),LFPbySweep(i,:),'Color',randColor);
    hold on;
end

subplot(3,1,2);
j=1;
for i=ts
    %randColor=[rand(1,1) rand(1,1) rand(1,1)];
    plot((1/LFP_Fs):(1/LFP_Fs):(1/LFP_Fs)*size(bandPassedLFPbySweep,2),bandPassedLFPbySweep(i,:),'Color',cs(j,:));
    hold on;
    j=j+1;
end

% Smooth, sum and rectify gamma-filtered LFP
bandPassedLFPbySweep(bandPassedLFPbySweep<0)=-bandPassedLFPbySweep(bandPassedLFPbySweep<0);
for i=1:size(bandPassedLFPbySweep,1)
    bandPassedLFPbySweep(i,:)=smooth(bandPassedLFPbySweep(i,:)',250)'; % Approx. 9 ms window for smoothing, given Fs=32000
end

maxes=zeros(length(ts),1);
j=1;
for i=ts
    maxes(j)=max(bandPassedLFPbySweep(i,floor(1/(1/LFP_Fs))+1:floor(1.1/(1/LFP_Fs))+1));
    j=j+1;
end
[y,newInds]=sort(maxes);
ts=ts(newInds);
    
color1='red';
color2='blue';

subplot(3,1,3);
j=1;
for i=ts
    if j>floor(length(ts)/2)
        c=color2;
    else
        c=color1;
    end
    %randColor=[rand(1,1) rand(1,1) rand(1,1)];
    plot((1/LFP_Fs):(1/LFP_Fs):(1/LFP_Fs)*size(bandPassedLFPbySweep,2),bandPassedLFPbySweep(i,:),'Color',c);
    hold on;
    j=j+1;
end