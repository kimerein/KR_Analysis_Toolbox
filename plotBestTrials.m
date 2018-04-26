function plotBestTrials(LFPbySweep,bandPassedLFPbySweep,LFP_Fs,trialClassifications,baselineWindow)

% Separate figures for trials with classification 8, 3 and 4
% Use sortSingleTrials.m to classify single trials
% Will plot overlaid trial pairs with 
% 1 LED on trial and 1 LED off trial
% The LED on trial will have classification 8,
% that is, high gamma preceding the LED pulse and low gamma (silencing)
% during/just after the LED pulse
% The LED off trial will have high gamma preceding the LED pulse
% and either high or low gamma during/just after the LED

totalTrialLength=(1/LFP_Fs)*size(LFPbySweep,2);
baselineInds=floor(baselineWindow(1)*LFP_Fs)+1:floor(baselineWindow(2)*LFP_Fs);
% Absolute value of gamma-filtered signal
bandPassedLFPbySweep(bandPassedLFPbySweep<0)=-bandPassedLFPbySweep(bandPassedLFPbySweep<0);
% Smoothed gamma-filtered LFP integral
smoothedBPLFPbySweep=zeros(size(bandPassedLFPbySweep));
for i=1:size(bandPassedLFPbySweep,1)
    smoothedBPLFPbySweep(i,:)=smooth(bandPassedLFPbySweep(i,:)',20)'; % Approx. 4.5 ms window for smoothing, given Fs=32000
end
for i=1:size(bandPassedLFPbySweep,1)
    for j=1:10
        smoothedBPLFPbySweep(i,:)=smooth(smoothedBPLFPbySweep(i,:)',20)'; % Approx. 4.5 ms window for smoothing, given Fs=32000
    end
end

% Plot averages
c=[8 4 3];
for i=c
    if any(trialClassifications==i)
        figure;
        title(['Average with Classification ' num2str(i)]);
        subplot(2,1,1);
        plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,mean(LFPbySweep(trialClassifications==i,:),1),'k');
        subplot(2,1,2);
        plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,mean(smoothedBPLFPbySweep(trialClassifications==i,:),1),'k');
    end
end

% For big experiments, can't look at all trials, so randomly sub-sample
% r = 1 + (size(LFPbySweep,1)+1-1).*rand(100,1);
% subSample=sort(unique(floor(r)));
% LFPbySweep=LFPbySweep(subSample,:);
% bandPassedLFPbySweep=bandPassedLFPbySweep(subSample,:);
% trialClassifications=trialClassifications(subSample);

% Plot single trials
c=[8 4 3];
for i=c
    if any(trialClassifications==i)
        figure;
        % Plan figure and trial spacing
        n=sum(trialClassifications==i);
        if n>50
            n=50;
        end
        theEnd=n*0.6;
        trialOffsets=0:0.6:theEnd;
        trialOffsets=trialOffsets';
        trialOffsets=trialOffsets*ones(1,size(LFPbySweep,2));
        theseInds=find(trialClassifications==i);
        plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,LFPbySweep(theseInds(1:n),:)-mean(LFPbySweep(theseInds(1:n),:),2)*ones(1,size(LFPbySweep,2))+trialOffsets(1:n,:),'Color','k');
        title(['Trials with Classification ' num2str(i)]);
    end
end
        
% Plot overlapping trials for comparison
n8=sum(trialClassifications==8);
n34=sum((trialClassifications==3)|(trialClassifications==4));
n=min(n8,n34);
if n>100
    n=100;
    a=find(trialClassifications==8);
    plotThese8=a(1:n);
    a=find((trialClassifications==3)|(trialClassifications==4));
    plotThese34=a(1:n);
else
    if n8>n34
        a=find(trialClassifications==8);
        plotThese8=a(1:n);
        plotThese34=find((trialClassifications==3)|(trialClassifications==4));
    else
        a=find((trialClassifications==3)|(trialClassifications==4));
        plotThese34=a(1:n);
        plotThese8=find(trialClassifications==8);
    end
end

if any(trialClassifications==8)
    figure;
    % Plan figure and trial spacing
    theEnd=n*0.6;
    trialOffsets=0:0.6:theEnd;
    trialOffsets=trialOffsets';
    trialOffsets=trialOffsets*ones(1,size(LFPbySweep,2));
    plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,LFPbySweep(plotThese8,:)-mean(LFPbySweep(plotThese8,:),2)*ones(1,size(LFPbySweep,2))+trialOffsets(1:n,:),'Color','b');
    hold on;
    plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,LFPbySweep(plotThese34,:)-mean(LFPbySweep(plotThese34,:),2)*ones(1,size(LFPbySweep,2))+trialOffsets(1:n,:),'Color','r');
end
% Gamma integral
if any(trialClassifications==8)
    figure;
    % Plan figure and trial spacing
    theEnd=n*0.1;
    trialOffsets=0:0.1:theEnd;
    trialOffsets=trialOffsets';
    trialOffsets=trialOffsets*ones(1,size(LFPbySweep,2));
    plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,smoothedBPLFPbySweep(plotThese8,:)-mean(smoothedBPLFPbySweep(plotThese8,:),2)*ones(1,size(LFPbySweep,2))+trialOffsets(1:n,:),'Color','b');
    hold on;
    plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,smoothedBPLFPbySweep(plotThese34,:)-mean(smoothedBPLFPbySweep(plotThese34,:),2)*ones(1,size(LFPbySweep,2))+trialOffsets(1:n,:),'Color','r');
end

% Averages comparison
figure;
subplot(2,1,1);
plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,mean(LFPbySweep(plotThese8,:),1)-mean(mean(LFPbySweep(plotThese8,baselineInds),1),2)*ones(1,size(LFPbySweep,2)),'Color','b');
hold on;
plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,mean(LFPbySweep(plotThese34,:),1)-mean(mean(LFPbySweep(plotThese34,baselineInds),1),2)*ones(1,size(LFPbySweep,2)),'Color','r');
subplot(2,1,2);
plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,mean(smoothedBPLFPbySweep(plotThese8,:),1)-mean(mean(smoothedBPLFPbySweep(plotThese8,baselineInds),1),2)*ones(1,size(LFPbySweep,2)),'Color','b');
hold on;
plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,mean(smoothedBPLFPbySweep(plotThese34,:),1)-mean(mean(smoothedBPLFPbySweep(plotThese34,baselineInds),1),2)*ones(1,size(LFPbySweep,2)),'Color','r');
% subplot(2,1,2);
% plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,mean(smoothedBPLFPbySweep(plotThese8,:),1),'Color','b');
% hold on;
% plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,mean(smoothedBPLFPbySweep(plotThese34,:),1),'Color','r');