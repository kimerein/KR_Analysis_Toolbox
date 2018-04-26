function plotLEDtriggeredLFP(data,ledBySweep,ledForSweeps,LEDthresh,Fs,trialClassifications,compareGroup1,compareGroup2,secondsToPlot,bandPass,lowerCutoff,upperCutoff,alignToBaseline,baselineWindow)

% Align to baseline?
baselineInds=floor(baselineWindow(1)*Fs)+1:floor(baselineWindow(2)*Fs);
if alignToBaseline
    data=data-mean(data(:,baselineInds),2)*ones(1,size(data,2));
end

%timeWindowInds=floor(timeWindowToPlot(1)*Fs)+1:floor(timeWindowToPlot(2)*Fs);
indsToPlot=secondsToPlot*Fs;
% Only use the specified trials
data=data((ismember(trialClassifications,compareGroup1)|ismember(trialClassifications,compareGroup2))'&ledForSweeps>0,:);
ledBySweep=ledBySweep((ismember(trialClassifications,compareGroup1)|ismember(trialClassifications,compareGroup2))'&ledForSweeps>0,:);
trialClassifications=trialClassifications((ismember(trialClassifications,compareGroup1)|ismember(trialClassifications,compareGroup2))'&ledForSweeps>0);

% Put trial groups into cell array (this code could be extended for comparing more
% than two groups of trial types)
trialTypeGroups=cell(2,1);
trialTypeGroups{1}=compareGroup1;
trialTypeGroups{2}=compareGroup2;

% Check to be sure there are trials that have received these
% classifications
for i=1:length(trialTypeGroups)
    currGroup=trialTypeGroups{i};
    noneOfTheseTrials=1;
    for j=1:length(currGroup)
        if any(trialClassifications==currGroup(j))
            noneOfTheseTrials=0;
            break;
        end
    end
    if noneOfTheseTrials==1
        disp(['There are no trials in group' num2str(i)]);
        trialTypeGroups{i}=[];
    end
end

% Band-pass
if bandPass
    data=bandPassLFP(data,Fs,lowerCutoff,upperCutoff,0);
end

% Align in time to LED onset
% Find LED onsets in indices
onsetInds=zeros(size(data,1),1);
for i=1:length(onsetInds)
    onsetInds(i)=find(ledBySweep(i,:)>LEDthresh,1,'first');
end
plotData=zeros(size(data,1),indsToPlot);
plotLedBySweep=zeros(size(data,1),indsToPlot);
for i=1:size(data,1)
    offset=onsetInds(i)-floor(indsToPlot/2);
    if offset<=0
        disp(['Time-window-to-plot too big to align trial ' num2str(i)]);
        offset=1;
    end
    plotData(i,:)=data(i,offset:offset+indsToPlot-1);
    plotLedBySweep(i,:)=ledBySweep(i,offset:offset+indsToPlot-1);
end

LEDonsetTime=(1/Fs)*(offset+floor(indsToPlot/2));
meanLEDindex=floor(mean(onsetInds));
startTime=(1/Fs)*(meanLEDindex-floor(indsToPlot/2));
endTime=(1/Fs)*(meanLEDindex+floor(indsToPlot/2));
if startTime<0 || endTime>(1/Fs)*size(data,2)
    disp('Use a smaller window to plot. Window should fit in each trial.');
    return
end

% Equalize the number of trials in each group for better comparison
[set1,set2]=matchNumberOfTrials(find(ismember(trialClassifications,trialTypeGroups{1})),find(ismember(trialClassifications,trialTypeGroups{2})));
trialGroups={set1,set2};

% Check alignment to LED pulse
figure;
plot(plotLedBySweep');

% Plot averages for group 1 and group 2 trials separately
figure;
for i=1:length(trialGroups)
    trialsInThisGroup=trialGroups{i};
    if ~isempty(trialsInThisGroup)
        if i==1
            c='red';
        else
            c='black';
        end
        subplot(length(trialTypeGroups)*2+1,1,(i*2)-1);
        plot(startTime:(endTime-startTime)/(indsToPlot-1):endTime,mean(plotData(trialsInThisGroup,:),1),c);
        axis tight;
        hold on;
        minVal=min(mean(plotData(trialsInThisGroup,:),1));
        maxVal=max(mean(plotData(trialsInThisGroup,:),1));
        plot([LEDonsetTime LEDonsetTime],[minVal maxVal],'blue');
        axis tight;
        
        subplot(length(trialTypeGroups)*2+1,1,i*2);
        axis tight;
        absValData=plotData(trialsInThisGroup,:);
        absValData(absValData<0)=-absValData(absValData<0);
        plot(startTime:(endTime-startTime)/(indsToPlot-1):endTime,sum(absValData,1),c);
        hold on;
        minVal=min(sum(absValData,1));
        maxVal=max(sum(absValData,1));
        plot([LEDonsetTime LEDonsetTime],[minVal maxVal],'blue');
        axis tight;
        
        subplot(length(trialTypeGroups)*2+1,1,length(trialTypeGroups)*2+1);
        plot(startTime:(endTime-startTime)/(indsToPlot-1):endTime,sum(absValData,1),c);
        hold on;
        minVal=min(sum(absValData,1));
        maxVal=max(sum(absValData,1));
        plot([LEDonsetTime LEDonsetTime],[minVal maxVal],'blue');
        axis tight;
    end
end

% Averages comparison
figure;
for i=1:length(trialGroups)
    if i==1
        c='red';
    else
        c='black';
    end
    trialsInThisGroup=trialGroups{i};
    if ~isempty(trialsInThisGroup)
        plot(startTime:(endTime-startTime)/(indsToPlot-1):endTime,mean(plotData(trialsInThisGroup,:),1),c);
        hold all;
    end
    minVal=min(mean(plotData(trialsInThisGroup,:),1));
    maxVal=max(mean(plotData(trialsInThisGroup,:),1));
    plot([LEDonsetTime LEDonsetTime],[minVal maxVal],'blue');
    axis tight;
end
