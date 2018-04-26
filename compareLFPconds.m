function [LFPfigs,redTrials,blueTrials]=compareLFPconds(LFPdata,stimsForSweeps,LFP_Fs,ti,params,group1_values,group2_values,compareSumOrAv,rectifyForStats)
% Compares LFP responses for two specific stimulus conditions
% 
% RETURNS:
% LFPfigs: figure handles of analysis
% redTrials: indices corresponding to LFPData rows used in the redTrials
% set
% blueTrials: indices corresponding to LFPData rows used in the blueTrials
% set
% blueTrials and redTrials sets are compared, statistically, using ttests
% on response windows tiled across each trial
% 
% PARAMETERS:
% LFPdata: matrix with samples as columns and trials as rows
% stimsForSweeps: stimsForSweeps: vector with each stimulus condition for all trials; rows
% correspond to trial rows of data
% LFP_Fs: sampling rate of data
% ti: figure title
% params: a structure with information about the setup of a trial
% must contain the field: totalTrialLength, specifying the length of a
% trial in seconds
% group1_values: a vector of the stimulus conditions to include in the redTrials set;
% all trials with any of these stimulus conditions will be included in redTrials
% group2_values: a vector of the stimulus conditions to include in the blueTrials set
% all trials with any of these stimulus conditions will be included in
% blueTrials
% compareSumOrAverage: if 'sum', sum over trials; if 'average', average
% over trials
% rectifyForStats: if 1, take absolute value of each trial in LFPdata
% before plotting or calculating statistics

LFPfigs=[];

redTrials=ismember(stimsForSweeps,group1_values);
blueTrials=ismember(stimForSweeps,group2_values);

% Plan figure layout
nfigRows=3;
nfigCols=1;

f1=figure;
LFPfigs=[LFPfigs; f1];

% Plot average stimulus-triggered LFP for different stimulus conditions -- 
% i.e., red trials vs. blue trials
% Check that have same number of red and blue trials to compare
% Else randomly sample from red or blue trials to get same number for
% comparison
if length(redTrials)~=length(blueTrials)
    if length(redTrials)>length(blueTrials)
        takeInds=randperm(length(redTrials));
        redTrials=redTrials(sort(takeInds(1:length(blueTrials))));
    else
        takeInds=randperm(length(blueTrials));
        blueTrials=blueTrials(sort(takeInds(1:length(redTrials))));
    end
end

% Plot red vs. blue trials without standard deviations
redTrials=redTrials';
blueTrials=blueTrials';
axesmatrix(nfigRows,nfigCols,1);
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,combineData(averageOrAdd,LFPdata(redTrials,:),1),'Color','r');
hold on;
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,combineData(averageOrAdd,LFPdata(blueTrials,:),1),'Color','b');
title([ti  ' ' compareSumOrAv]);

% Plot red vs. blue trials with standard deviations
redTrials=redTrials';
blueTrials=blueTrials';
axesmatrix(nfigRows,nfigCols,2);
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,combineData(averageOrAdd,LFPdata(redTrials,:),1),'Color','r');
hold on;
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,combineData(averageOrAdd,LFPdata(redTrials,:),1)+std(LFPdata(redTrials,:),1),'Color','m');
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,combineData(averageOrAdd,LFPdata(redTrials,:),1)-std(LFPdata(redTrials,:),1),'Color','m');
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,combineData(averageOrAdd,LFPdata(blueTrials,:),1),'Color','b');
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,combineData(averageOrAdd,LFPdata(blueTrials,:),1)+std(LFPdata(blueTrials,:),1),'Color','c');
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,combineData(averageOrAdd,LFPdata(blueTrials,:),1)-std(LFPdata(blueTrials,:),1),'Color','c');
title([ti ' - With Std. Devs. ' compareSumOrAv]);

if rectifyForStats==0
    % Calculate red vs. blue regions with statistically significant differences
    % (t-test) and plot on figure
    axesmatrix(nfigRows,nfigCols,3);
    plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,mean(LFPdata(redTrials,:),1),'Color','r');
    hold on;
    plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,mean(LFPdata(blueTrials,:),1),'Color','b');
    minValForPlot=min(min([mean(LFPdata(redTrials,:),1); mean(LFPdata(blueTrials,:),1)],[],2),[],1);
    binSize=0.01;
    for i=binSize:binSize:params.totalTrialLength
        redComp=mean(LFPdata(redTrials,floor((i-binSize)/(params.totalTrialLength/(size(LFPdata,2)-1)))+1:floor(i/(params.totalTrialLength/(size(LFPdata,2)-1)))+1),2);
        blueComp=mean(LFPdata(blueTrials,floor((i-binSize)/(params.totalTrialLength/(size(LFPdata,2)-1)))+1:floor(i/(params.totalTrialLength/(size(LFPdata,2)-1)))+1),2);
        [h,p,ci]=ttest2(redComp,blueComp);
        if h==1
            if mean(redComp)>mean(blueComp)
                currCol='r';
            else
                currCol='b';
            end
            line([i-binSize i],[minValForPlot-(-minValForPlot)/5 minValForPlot-(-minValForPlot)/5],'Color',currCol);
        end
    end
    title([ti ' - Red v. Blue,Black=Significance,Non-Rectified Signal']);
else
    % For filtered LFP sweeps, get significance by comparing the
    % rectified signal in time bin 
    axesmatrix(nfigRows,nfigCols,3);
    plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,mean(LFPdata(redTrials,:),1),'Color','r');
    hold on;
    plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,mean(LFPdata(blueTrials,:),1),'Color','b');
    minValForPlot=min(min([mean(LFPdata(redTrials,:),1); mean(LFPdata(blueTrials,:),1)],[],2),[],1);
    binSize=0.04;
    for i=binSize:binSize:params.totalTrialLength
        thisRedSeg=LFPdata(redTrials,floor((i-binSize)/(params.totalTrialLength/(size(LFPdata,2)-1)))+1:floor(i/(params.totalTrialLength/(size(LFPdata,2)-1)))+1);
        thisRedSeg(thisRedSeg<0)=-thisRedSeg(thisRedSeg<0);
        thisRedSeg=sum(thisRedSeg,2);
        thisBlueSeg=LFPdata(blueTrials,floor((i-binSize)/(params.totalTrialLength/(size(LFPdata,2)-1)))+1:floor(i/(params.totalTrialLength/(size(LFPdata,2)-1)))+1);
        thisBlueSeg(thisBlueSeg<0)=-thisBlueSeg(thisBlueSeg<0);
        thisBlueSeg=sum(thisBlueSeg,2);
        [h,p,ci]=ttest2(thisRedSeg,thisBlueSeg);
        if h==1
            if mean(thisRedSeg)>mean(thisBlueSeg)
                currCol='r';
            else
                currCol='b';
            end
            line([i-binSize i],[minValForPlot-(-minValForPlot)/5 minValForPlot-(-minValForPlot)/5],'Color',currCol);
        end
    end
    title([ti ' - Red v. Blue,Black=Significance,Rectified Signal']);
end
end
    
function r=combineData(averageOrAdd,data,dim)
    if strcmp(averageOrAdd,'average')
        r=mean(data,dim);
    else
        r=sum(data,dim);
    end
end
