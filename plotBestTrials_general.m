function plotBestTrials_general(LFPbySweep,bandPassedLFPbySweep,LFP_Fs,trialClassifications,baselineWindow,compareGroup1,compareGroup2)
% Compares different trial types (defined in terms of their response
% features)
% CompareGroup1 is plotted in blue
% CompareGroup2 is plotted in red
%
% Will make separate figures for trials with different classifications
% Use sortSingleTrials.m to classify single trials
% Will compare trials with classifications contained in compareGroup1 with
% trials with classifications contained in compareGroup2
% Will plot separate trial pairs, overlaid single trial pairs, averages, 
% and overlaid averages for the compared groups
% 
% PARAMETERS:
% LFPbySweep: matrix with samples as columns and trials as rows
% bandPassedLFPbySweep: matrix with samples as columns and trials as rows,
% band-passed data of LFPbySweep
% LFP_Fs: sampling rate of LFPbySweep and bandPassedLFPbySweep
% trialClassifications: a column vector with one entry for every row of
% LFPbySweep (every trial), each element i gives the trial classification
% (see sortSingleTrials.m) for the corresponding ith row (trial) of LFPbySweep/bandPassedLFPbySweep
% baselineWindow: a two-element vector, first element gives the start time
% in seconds
% for a baseline period (will be used to align trials) relative to the
% start of the trial, second element gives the end of this baseline period
% in seconds relative to start of trial
% compareGroup1: a vector containing the trial classification #s (see
% sortSingleTrials.m) to include in Trial Types Group 1
% compareGroup2: a vector containing the trial classification #s (see
% sortSingleTrials.m) to include in Trial Types Group 2
% Groups 1 and 2 will be compared

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

% Get total trial length and baseline window
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

% Plot averages for group 1 and group 2 trials separately
for i=1:length(trialTypeGroups)
    currGroup=trialTypeGroups{i};
    if ~isempty(currGroup)
        figure;
        title(['Average with Classification(s) ' num2str(currGroup)]);
        trialsInThisGroup=ismember(trialClassifications,currGroup);
        subplot(2,1,1);
        plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,mean(LFPbySweep(trialsInThisGroup,:),1),'k');
        subplot(2,1,2);
        plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,mean(smoothedBPLFPbySweep(trialsInThisGroup,:),1),'k');
    end
end

% For big experiments, can't look at all trials, so randomly sub-sample
% r = 1 + (size(LFPbySweep,1)+1-1).*rand(100,1);
% subSample=sort(unique(floor(r)));
% LFPbySweep=LFPbySweep(subSample,:);
% bandPassedLFPbySweep=bandPassedLFPbySweep(subSample,:);
% trialClassifications=trialClassifications(subSample);

% Plot single trials
for i=1:length(trialTypeGroups)
    currGroup=trialTypeGroups{i};
    if ~isempty(currGroup)
        figure;
        % Plan figure and trial spacing
        n=sum(ismember(trialClassifications,currGroup));
        if n>50
            n=50;
        end
        theEnd=n*0.6;
        trialOffsets=0:0.6:theEnd;
        trialOffsets=trialOffsets';
        trialOffsets=trialOffsets*ones(1,size(LFPbySweep,2));
        theseInds=find(ismember(trialClassifications,currGroup));
        plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,LFPbySweep(theseInds(1:n),:)-mean(LFPbySweep(theseInds(1:n),:),2)*ones(1,size(LFPbySweep,2))+trialOffsets(1:n,:),'Color','k');
        title(['Trials with Classification(s) ' num2str(currGroup)]);
    end
end
        
% Plot overlapping trials for comparison
n8=sum(ismember(trialClassifications,trialTypeGroups{1}));
n34=sum(ismember(trialClassifications,trialTypeGroups{2}));
n=min(n8,n34);
if n8==0
    n=n34;
elseif n34==0
    n=n8;
end
if n>100
    n=100;
    if ~isempty(trialTypeGroups{1})
        a=find(ismember(trialClassifications,trialTypeGroups{1}));
        plotThese8=a(1:n);
    else
        plotThese8=[];
    end
    if ~isempty(trialTypeGroups{2})
        a=find(ismember(trialClassifications,trialTypeGroups{2}));
        plotThese34=a(1:n);
    else
        plotThese34=[];
    end
else
    if n8>n34
        if ~isempty(trialTypeGroups{1})
            a=find(ismember(trialClassifications,trialTypeGroups{1}));
            plotThese8=a(1:n);
        else
            plotThese8=[];
        end
        if ~isempty(trialTypeGroups{2})
            plotThese34=find(ismember(trialClassifications,trialTypeGroups{2}));
        else
            plotThese34=[];
        end
    else
        if ~isempty(trialTypeGroups{2})
            a=find(ismember(trialClassifications,trialTypeGroups{2}));
            plotThese34=a(1:n);
        else
            plotThese34=[];
        end
        if ~isempty(trialTypeGroups{1})
            plotThese8=find(ismember(trialClassifications,trialTypeGroups{1}));
        else
            plotThese8=[];
        end
    end
end

if ~isempty(trialTypeGroups{1}) || ~isempty(trialTypeGroups{2})  
    figure;
    % Plan figure and trial spacing
    theEnd=n*0.6;
    trialOffsets=0:0.6:theEnd;
    trialOffsets=trialOffsets';
    trialOffsets=trialOffsets*ones(1,size(LFPbySweep,2));
    if ~isempty(trialTypeGroups{1})
        plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,LFPbySweep(plotThese8,:)-mean(LFPbySweep(plotThese8,:),2)*ones(1,size(LFPbySweep,2))+trialOffsets(1:n,:),'Color','b');
    end
    hold on;
    if ~isempty(trialTypeGroups{2})
        plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,LFPbySweep(plotThese34,:)-mean(LFPbySweep(plotThese34,:),2)*ones(1,size(LFPbySweep,2))+trialOffsets(1:n,:),'Color','r');
    end
end
% Gamma integral
if ~isempty(trialTypeGroups{1}) || ~isempty(trialTypeGroups{2})  
    figure;
    % Plan figure and trial spacing
    theEnd=n*0.1;
    trialOffsets=0:0.1:theEnd;
    trialOffsets=trialOffsets';
    trialOffsets=trialOffsets*ones(1,size(LFPbySweep,2));
    if ~isempty(trialTypeGroups{1})
        plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,smoothedBPLFPbySweep(plotThese8,:)-mean(smoothedBPLFPbySweep(plotThese8,:),2)*ones(1,size(LFPbySweep,2))+trialOffsets(1:n,:),'Color','b');
    end
    hold on;
    if ~isempty(trialTypeGroups{2})
        plot(0:totalTrialLength/(size(LFPbySweep,2)-1):totalTrialLength,smoothedBPLFPbySweep(plotThese34,:)-mean(smoothedBPLFPbySweep(plotThese34,:),2)*ones(1,size(LFPbySweep,2))+trialOffsets(1:n,:),'Color','r');
    end
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