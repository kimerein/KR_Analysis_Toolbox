function plotHighGammaTrials(LFPbySweep,bandPassedLFPbySweep,ledForSweeps,LFP_Fs,totalTrialLength,baselineWindow,gammaAssessWindow)

% Discard first and last trials
LFPbySweep=LFPbySweep(2:end-1,:);
bandPassedLFPbySweep=bandPassedLFPbySweep(2:end-1,:);
ledForSweeps=ledForSweeps(2:end-1);

assessInds=floor(gammaAssessWindow(1)*LFP_Fs)+1:floor(gammaAssessWindow(2)*LFP_Fs);

% Align trials to baseline
baseInds=floor(baselineWindow(1)*LFP_Fs)+1:floor(baselineWindow(2)*LFP_Fs);
for i=1:size(LFPbySweep,1)
    baseVal=mean(mean(LFPbySweep(i,baseInds),1),2);
    LFPbySweep(i,:)=LFPbySweep(i,:)-baseVal;
end


% Get gamma power proxy
bandPassedLFPbySweep(bandPassedLFPbySweep<0)=-bandPassedLFPbySweep(bandPassedLFPbySweep<0);

% Gamma power proxy stats.
disp('mean');
disp(mean(mean(bandPassedLFPbySweep(:,assessInds),2),1));
disp('stdev');
disp(std(mean(bandPassedLFPbySweep(:,assessInds),2)));

% Threshold for "high gamma trial" classification
%highGammaThresh=0.040;
highGammaThresh=0.022;

highGammaTrials=zeros(1,size(LFPbySweep,1));
for i=1:size(LFPbySweep,1)
    if mean(bandPassedLFPbySweep(i,assessInds))>highGammaThresh
        highGammaTrials(i)=1;
    else
        highGammaTrials(i)=0;
    end
end

theEnd=(size(LFPbySweep,1)+1)*0.6;
trialOffsets=0:0.6:theEnd;
trialOffsets=trialOffsets';
trialOffsets=trialOffsets*ones(1,size(LFPbySweep,2));

if any(highGammaTrials)
    figure;
    %subplot(3,1,1);
    n=length(find(highGammaTrials&(~(ledForSweeps>0))));
    plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,LFPbySweep(highGammaTrials&(~(ledForSweeps>0)),:)+trialOffsets(1:n,:),'Color','r');
    figure;
    %subplot(3,1,2);
    n=length(find(highGammaTrials&(ledForSweeps>0)));
    plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,LFPbySweep(highGammaTrials&(ledForSweeps>0),:)+trialOffsets(1:n,:),'Color','b');
    figure;
    %subplot(3,1,3);
    n=length(find(highGammaTrials&(~(ledForSweeps>0))));
    plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,LFPbySweep(highGammaTrials&(~(ledForSweeps>0)),:)+trialOffsets(1:n,:),'Color','r');
    hold on;
    n=length(find(highGammaTrials&((ledForSweeps>0))));
    plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,LFPbySweep(highGammaTrials&((ledForSweeps>0)),:)+trialOffsets(1:n,:),'Color','b');
end
   
% if any(~highGammaTrials)
%     figure;
%     subplot(3,1,1);
%     n=length(find(~highGammaTrials&(~(ledForSweeps>0))));
%     plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,LFPbySweep(~highGammaTrials&(~(ledForSweeps>0)),:)+trialOffsets(1:n,:),'Color','r');
%     subplot(3,1,2);
%     n=length(find(~highGammaTrials&(ledForSweeps>0)));
%     plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,LFPbySweep(~highGammaTrials&(ledForSweeps>0),:)+trialOffsets(1:n,:),'Color','b');
%     subplot(3,1,3);
%     n=length(find(~highGammaTrials&(~(ledForSweeps>0))));
%     plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,LFPbySweep(~highGammaTrials&(~(ledForSweeps>0)),:)+trialOffsets(1:n,:),'Color','r');
%     hold on;
%     n=length(find(~highGammaTrials&((ledForSweeps>0))));
%     plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,LFPbySweep(~highGammaTrials&((ledForSweeps>0)),:)+trialOffsets(1:n,:),'Color','b');
% end

% Smoothed gamma-filtered LFP integral
for i=1:size(bandPassedLFPbySweep,1)
    for j=1:3
        bandPassedLFPbySweep(i,:)=smooth(bandPassedLFPbySweep(i,:)',40)'; % Approx. 9 ms window for smoothing, given Fs=32000
    end
end
% Align trials to baseline
baseInds=floor(0.2*LFP_Fs)+1:floor(0.4*LFP_Fs);
for i=1:size(LFPbySweep,1)
    baseVal=mean(mean(bandPassedLFPbySweep(i,baseInds),1),2);
    bandPassedLFPbySweep(i,:)=bandPassedLFPbySweep(i,:)-baseVal;
end

% n1=length(find(highGammaTrials&(~(ledForSweeps>0))));
% n2=length(find(highGammaTrials&((ledForSweeps>0))));
% if n1>n2
%     n2

if any(highGammaTrials)
    figure;
    n=length(find(highGammaTrials&(~(ledForSweeps>0))));
    plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,bandPassedLFPbySweep(highGammaTrials&(~(ledForSweeps>0)),:),'Color','r');
    hold on;
    useThese=find(highGammaTrials&((ledForSweeps>0)));
    for i=1:n
        plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,bandPassedLFPbySweep(useThese(i),:),'Color','b');
        hold on;
    end
end



