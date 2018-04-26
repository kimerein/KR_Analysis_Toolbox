function plotGammaPowerProxy(bandPassedLFPbySweep,ledForSweeps)
% Plots bandPassedLFPbySweep

backup=bandPassedLFPbySweep;

totalTrialLength=3.2;

% figure;
% plot(0:totalTrialLength/(size(bandPassedLFPbySweep,2)-1):totalTrialLength,mean(bandPassedLFPbySweep(ledForSweeps>0,:),1),'Color','b');
% hold on
% plot(0:totalTrialLength/(size(bandPassedLFPbySweep,2)-1):totalTrialLength,mean(bandPassedLFPbySweep(~(ledForSweeps>0),:),1),'Color','k');
% axis tight

% Rectified gamma-filtered LFP
figure;
bandPassedLFPbySweep(bandPassedLFPbySweep<0)=-bandPassedLFPbySweep(bandPassedLFPbySweep<0);
plot(0:totalTrialLength/(size(bandPassedLFPbySweep,2)-1):totalTrialLength,mean(bandPassedLFPbySweep(ledForSweeps>0,:),1),'Color','b');
hold on
plot(0:totalTrialLength/(size(bandPassedLFPbySweep,2)-1):totalTrialLength,mean(bandPassedLFPbySweep(~(ledForSweeps>0),:),1),'Color','k');
axis tight

% Smoothed gamma-filtered LFP integral
for i=1:size(bandPassedLFPbySweep,1)
    bandPassedLFPbySweep(i,:)=smooth(bandPassedLFPbySweep(i,:)',250)'; % Approx. 9 ms window for smoothing, given Fs=32000
end

figure;
subplot(3,1,1);
plot(0:totalTrialLength/(size(backup,2)-1):totalTrialLength,mean(backup(ledForSweeps>0,:),1),'Color','b');
hold on
plot(0:totalTrialLength/(size(backup,2)-1):totalTrialLength,mean(backup(~(ledForSweeps>0),:),1),'Color','k');
axis tight
subplot(3,1,2);
plot(0:totalTrialLength/(size(bandPassedLFPbySweep,2)-1):totalTrialLength,mean(bandPassedLFPbySweep(ledForSweeps>0,:),1),'Color','b');
hold on
plot(0:totalTrialLength/(size(bandPassedLFPbySweep,2)-1):totalTrialLength,mean(bandPassedLFPbySweep(~(ledForSweeps>0),:),1),'Color','k');
axis tight