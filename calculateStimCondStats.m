function calculateStimCondStats(data,stimsForSweeps,LFP_Fs,baselineWindow)
% Calculates significance of differences between responses for different
% stimulus conditions; pairwise comparisons of response distributions for
% all stimulus conditions; ttest stats., significance is plotted below
% each time window on comparison graph of superimposed means
%
% data: matrix with samples as columns and trials as rows; this is data
% from which response distributions will be drawn
% stimsForSweeps: vector with each stimulus condition for all trials; rows
% correspond to trial rows of data
% LFP_Fs: sampling rate of data
% baselineWindow: time window relative to trial beginning to use to
% calculate baseline

binWidth=100; % in indices
binOffset=50; % in indices

totalTrialLength=(1/LFP_Fs)*size(data,2);
times=0:totalTrialLength/(size(data,2)-1):totalTrialLength;
baselineInds=floor(baselineWindow(1)*LFP_Fs)+1:floor(baselineWindow(2)*LFP_Fs);

% Align all trials to baseline
data=data-mean(data(:,baselineInds),2)*ones(1,size(data,2));

allStims=sort(unique(stimsForSweeps));

for i=1:length(allStims)
    for j=i:length(allStims)
        figure;
        set1=stimsForSweeps==allStims(i);
        set2=stimsForSweeps==allStims(j);
        maxall=max(max(mean(data(set1,:),1)),max(mean(data(set2,:),1)));
        minall=min(min(mean(data(set1,:),1)),min(mean(data(set2,:),1)));
        plot(times,mean(data(set1,:),1),'Color','b');
        hold on;
        plot(times,mean(data(set2,:),1),'Color','r');
        [sigSignal_lessThanPoint05,sigSignal_lessThanPoint01,sigSignal_lessThanPoint001,indsForBins]=ttestSignals(data(set1,:),data(set2,:),binWidth,binOffset);
        for k=1:length(sigSignal_lessThanPoint05)
            if sigSignal_lessThanPoint001(i)~=0
                if k==length(sigSignal_lessThanPoint001)
                    line([times(indsForBins(k)) times(end)],[minall-(abs(maxall-minall))/30 minall-(abs(maxall)-abs(minall))/30],'Color',[0.1 0.1 0.1],'LineWidth',6);
                else
                    line([times(indsForBins(k)) times(indsForBins(k+1))],[minall-(abs(maxall)-abs(minall))/30 minall-(abs(maxall)-abs(minall))/30],'Color',[0.1 0.1 0.1],'LineWidth',6);
                end
                continue
            end
            if sigSignal_lessThanPoint01(k)~=0
                if k==length(sigSignal_lessThanPoint01)
                    line([times(indsForBins(k)) times(end)],[minall-(abs(maxall-minall))/30 minall-(abs(maxall)-abs(minall))/30],'Color',[0.5 0.5 0.5],'LineWidth',6);
                else
                    line([times(indsForBins(k)) times(indsForBins(k+1))],[minall-(abs(maxall)-abs(minall))/30 minall-(abs(maxall)-abs(minall))/30],'Color',[0.5 0.5 0.5],'LineWidth',6);
                end
                continue
            end
            if sigSignal_lessThanPoint05(k)~=0
                if k==length(sigSignal_lessThanPoint05)
                    line([times(indsForBins(k)) times(end)],[minall-(abs(maxall-minall))/30 minall-(abs(maxall)-abs(minall))/30],'Color',[0.75 0.75 0.75],'LineWidth',6);
                else
                    line([times(indsForBins(k)) times(indsForBins(k+1))],[minall-(abs(maxall)-abs(minall))/30 minall-(abs(maxall)-abs(minall))/30],'Color',[0.75 0.75 0.75],'LineWidth',6);
                end
            end
        end
    end
end
