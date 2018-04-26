function [window,led]=find_UPDOWN_cutoff_compareWLFPpower(spikes,led,window,expt)

% window={[0 0.5]; ...
%         [0.25 0.75]; ...
%         [0.5 1]; ...
%         [0.75 1.25]; ...
%         [1 1.5]; ...
%         [1.25 1.75]; ...
%         [1.5 2]; ...
%         [1.75 2.25]; ...
%         [2 2.5]; ...
%         [2.25 2.75]; ...
%         [2.5 3]};
% led={[1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]; ...
%         [1 3 8 9]};
% window={[1 1.5]+0.05; ...
%         [1.1 1.6]+0.05; ...
%         [1.2 1.7]+0.05; ...
%         [1.3 1.8]+0.05};
% led={[2 7]; ...
%         [2 7]; ...
%         [2 7]; ...
%         [2 7]};

windowWidth=0.5;
windowOffset=0.25;
windowTimewindow=[0 3];
trialTimewindow=[0 6];
ledVals=[1 3];
bins=100;
UP_spiking_thresh=404.5;
fileInd=unique(spikes.fileInd);

lowerBand=[0 30];
upperBand=[30 100];

a=windowTimewindow(1);
b=windowTimewindow(2);
i=1;
while a<b
    if a+windowWidth>b
        break
    end
    window{i}=[a a+windowWidth];
    led{i}=ledVals;
    a=a+windowOffset;
    i=i+1;
end

Fs=32000;
downSampFactor=10;
physCh=[23 15 22 14 19 11 20 12 17 9 16 8 18 10 21 13];
usePhysCh_ind=1;
% dataDir='Z:\New Acquisition Computer - Starting 11-2011\Raw Data Backups\KR Data\Data\RawData\';
dataDir='W:\New Acquisition Computer\';
daqFileNames=expt.files.names(fileInd);
LFPbySweep=getJustLFPbySweep(dataDir,daqFileNames,Fs,downSampFactor,physCh(usePhysCh_ind));
Fs=Fs/downSampFactor;
if any(size(LFPbySweep{1})==0)
    disp('No data in this daq file.');
end
ledConds=expt.sweeps.led(ismember(expt.sweeps.fileInd,fileInd));
if length(ledConds)~=size(LFPbySweep{1},1)
    disp('ledConds and acquired LFP sweeps do not match');
    return
end
temp=LFPbySweep{1};
subLFPbySweep=temp(ismember(ledConds,ledVals),:);
subSpikes=filtspikes(spikes,0,'led',ledVals);
disp('These numbers should be the same');
disp(size(subLFPbySweep,1));
disp(length(unique(subSpikes.trials)));

alln=[];
times=trialTimewindow(1):1/Fs:trialTimewindow(2)-(1/Fs);
if length(times)~=size(subLFPbySweep,2)
    disp('times length is wrong');
    return
end
j_up=1;
j_down=1;
up_mua=[];
down_mua=[];
for i=1:length(led)
    [m,s,n]=calcMeanAndStdDuringWindow(filtspikes(spikes,0,'led',led{i}),window{i});
    if length(n)~=size(subLFPbySweep,1)
        disp('Trial counts do not match');
        break
    end
    up_inds=n>=UP_spiking_thresh;
    down_inds=n<UP_spiking_thresh;
    up_mua=[up_mua; n(up_inds)];
    down_mua=[down_mua; n(down_inds)];
    w=window{i};
    up_LFP(j_up:j_up+sum(up_inds)-1,:)=subLFPbySweep(up_inds,times>=w(1) & times<w(2));
    j_up=j_up+sum(up_inds);
    down_LFP(j_down:j_down+sum(down_inds)-1,:)=subLFPbySweep(down_inds,times>=w(1) & times<w(2));
    j_down=j_down+sum(down_inds);
    alln=[alln; n];
    disp(i);
end

[heights,centers]=hist(alln,bins);
centers=[centers(1)-(centers(2)-centers(1)) centers centers(end)+(centers(2)-centers(1))];
heights=[0 heights 0];
figure(); 
hist(alln,bins);
hold on;
plot(centers,heights);
% fitted=fit(centers',heights','gauss2');
% fitx=centers(1)-(centers(2)-centers(1)):1:centers(end)+(centers(2)-centers(1));
% y1=fitted.a1.*exp(-((fitx-fitted.b1)./fitted.c1).^2);
% plot(fitx,y1,'Color','r');
% y2=fitted.a2.*exp(-((fitx-fitted.b2)./fitted.c2).^2);
% plot(fitx,y2,'Color','g');

up_pr=makePowerSpectrum_returnPowerRatio(up_LFP,Fs,lowerBand,upperBand);
down_pr=makePowerSpectrum_returnPowerRatio(down_LFP,Fs,lowerBand,upperBand);

figure(); 
scatter(up_mua,up_pr,[],'g');
figure(); 
scatter(down_mua,down_pr,[],'r');
figure(); 
scatter(down_mua,down_pr,[],'r');
hold on;
scatter(up_mua,up_pr,[],'g');

all_pr=[up_pr; down_pr];
[heights,centers]=hist(all_pr,bins);
centers=[centers(1)-(centers(2)-centers(1)) centers centers(end)+(centers(2)-centers(1))];
heights=[0 heights 0];
figure(); 
hist(all_pr,bins);
hold on;
plot(centers,heights);
fitted=fit(centers',heights','gauss2');
fitx=centers(1)-(centers(2)-centers(1)):1:centers(end)+(centers(2)-centers(1));
y1=fitted.a1.*exp(-((fitx-fitted.b1)./fitted.c1).^2);
plot(fitx,y1,'Color','r');
y2=fitted.a2.*exp(-((fitx-fitted.b2)./fitted.c2).^2);
plot(fitx,y2,'Color','g');