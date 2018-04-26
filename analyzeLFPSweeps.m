function [LFPfigs,usedRed,usedBlue]=analyzeLFPSweeps(LFPdata,bandPassedLFP,LFP_Fs,redTrials,blueTrials,params)
% Obsolete LFP analysis
% Replaced by automatic experiment analysis package

LFPfigs=[];

% Align all LFP traces to initial value = 0
for i=1:size(LFPdata,1)
    initVal=LFPdata(i,1);
    LFPdata(i,:)=LFPdata(i,:)-initVal;
end
% Align band-passed LFP data to mean, rather than initial, value
mVals=mean(bandPassedLFP,2);
for i=1:size(bandPassedLFP,2)
    bandPassedLFP(:,i)=bandPassedLFP(:,i)-mVals;
end

% Plan figure layout
nfigRows=3;
nfigCols=3;
% figParams.matpos=[0 0 0 0];
% figParams.figmargin=[0.05 0.05 0.05 0.05];
% figParams.matmargin=[0 0 0 0];
% figParams.cellmargin=[0 0 0 0];

% Get alpha band-passed LFP
alphaLFP=bandPassLFP(LFPdata,LFP_Fs,8,12,1);
mVals=mean(alphaLFP,2);
for i=1:size(alphaLFP,2)
    alphaLFP(:,i)=alphaLFP(:,i)-mVals;
end

% Plot average stimulus-triggered LFP for raw and band-passed LFP
f1=figure;
LFPfigs=[LFPfigs; f1];
axesmatrix(nfigRows,nfigCols,1);
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,mean(LFPdata,1),'Color','k');
title('Av. Stim.-Triggered LFP - Init. Values Aligned');
axesmatrix(nfigRows,nfigCols,4);
plot(0:params.totalTrialLength/(size(bandPassedLFP,2)-1):params.totalTrialLength,mean(bandPassedLFP,1),'Color','k');
title('Av. Stim.-Triggered LFP - Gamma Filtered (30-80 Hz)');
axesmatrix(nfigRows,nfigCols,7);
plot(0:params.totalTrialLength/(size(alphaLFP,2)-1):params.totalTrialLength,mean(alphaLFP,1),'Color','k');
title('Av. Stim.-Triggered LFP - Alpha Filtered (8-12 Hz)');

% Plot average stimulus-triggered LFP for different "brain states" -- 
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
usedRed=redTrials;
usedBlue=blueTrials;

% Plot red vs. blue trials with standard deviations
redTrials=redTrials';
blueTrials=blueTrials';
axesmatrix(nfigRows,nfigCols,2);
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,mean(LFPdata(redTrials,:),1),'Color','r');
hold on;
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,mean(LFPdata(redTrials,:),1)+std(LFPdata(redTrials,:),1),'Color','m');
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,mean(LFPdata(redTrials,:),1)-std(LFPdata(redTrials,:),1),'Color','m');
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,mean(LFPdata(blueTrials,:),1),'Color','b');
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,mean(LFPdata(blueTrials,:),1)+std(LFPdata(blueTrials,:),1),'Color','c');
plot(0:params.totalTrialLength/(size(LFPdata,2)-1):params.totalTrialLength,mean(LFPdata(blueTrials,:),1)-std(LFPdata(blueTrials,:),1),'Color','c');
title('Av. Stim.-Triggered LFP - Red vs. Blue Trials');

axesmatrix(nfigRows,nfigCols,5);
plot(0:params.totalTrialLength/(size(bandPassedLFP,2)-1):params.totalTrialLength,mean(bandPassedLFP(redTrials,:),1),'Color','r');
hold on;
plot(0:params.totalTrialLength/(size(bandPassedLFP,2)-1):params.totalTrialLength,mean(bandPassedLFP(redTrials,:),1)+std(bandPassedLFP(redTrials,:),1),'Color','m');
plot(0:params.totalTrialLength/(size(bandPassedLFP,2)-1):params.totalTrialLength,mean(bandPassedLFP(redTrials,:),1)-std(bandPassedLFP(redTrials,:),1),'Color','m');
plot(0:params.totalTrialLength/(size(bandPassedLFP,2)-1):params.totalTrialLength,mean(bandPassedLFP(blueTrials,:),1),'Color','b');
plot(0:params.totalTrialLength/(size(bandPassedLFP,2)-1):params.totalTrialLength,mean(bandPassedLFP(blueTrials,:),1)+std(bandPassedLFP(blueTrials,:),1),'Color','c');
plot(0:params.totalTrialLength/(size(bandPassedLFP,2)-1):params.totalTrialLength,mean(bandPassedLFP(blueTrials,:),1)-std(bandPassedLFP(blueTrials,:),1),'Color','c');
title('Av. Stim.-Triggered LFP - Gamma, Red v. Blue');

axesmatrix(nfigRows,nfigCols,8);
plot(0:params.totalTrialLength/(size(alphaLFP,2)-1):params.totalTrialLength,mean(alphaLFP(redTrials,:),1),'Color','r');
hold on;
plot(0:params.totalTrialLength/(size(alphaLFP,2)-1):params.totalTrialLength,mean(alphaLFP(redTrials,:),1)+std(alphaLFP(redTrials,:),1),'Color','m');
plot(0:params.totalTrialLength/(size(alphaLFP,2)-1):params.totalTrialLength,mean(alphaLFP(redTrials,:),1)-std(alphaLFP(redTrials,:),1),'Color','m');
plot(0:params.totalTrialLength/(size(alphaLFP,2)-1):params.totalTrialLength,mean(alphaLFP(blueTrials,:),1),'Color','b');
plot(0:params.totalTrialLength/(size(alphaLFP,2)-1):params.totalTrialLength,mean(alphaLFP(blueTrials,:),1)+std(alphaLFP(blueTrials,:),1),'Color','c');
plot(0:params.totalTrialLength/(size(alphaLFP,2)-1):params.totalTrialLength,mean(alphaLFP(blueTrials,:),1)-std(alphaLFP(blueTrials,:),1),'Color','c');
title('Av. Stim.-Triggered LFP - Alpha, Red v. Blue');

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
title('Red v. Blue, Black=Significance');

% For filtered LFP sweeps, get significance by comparing the INTEGRAL of
% rectified signal in time bin (sum of rectified data points in time bin)
axesmatrix(nfigRows,nfigCols,6);
plot(0:params.totalTrialLength/(size(bandPassedLFP,2)-1):params.totalTrialLength,mean(bandPassedLFP(redTrials,:),1),'Color','r');
hold on;
plot(0:params.totalTrialLength/(size(bandPassedLFP,2)-1):params.totalTrialLength,mean(bandPassedLFP(blueTrials,:),1),'Color','b');
minValForPlot=min(min([mean(bandPassedLFP(redTrials,:),1); mean(bandPassedLFP(blueTrials,:),1)],[],2),[],1);
binSize=0.04;
for i=binSize:binSize:params.totalTrialLength
    thisRedSeg=bandPassedLFP(redTrials,floor((i-binSize)/(params.totalTrialLength/(size(bandPassedLFP,2)-1)))+1:floor(i/(params.totalTrialLength/(size(bandPassedLFP,2)-1)))+1);
    thisRedSeg(thisRedSeg<0)=-thisRedSeg(thisRedSeg<0);
    thisRedSeg=sum(thisRedSeg,2);
    thisBlueSeg=bandPassedLFP(blueTrials,floor((i-binSize)/(params.totalTrialLength/(size(bandPassedLFP,2)-1)))+1:floor(i/(params.totalTrialLength/(size(bandPassedLFP,2)-1)))+1);
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
title('Gamma, Red v. Blue, Black=Significance');

axesmatrix(nfigRows,nfigCols,9);
plot(0:params.totalTrialLength/(size(alphaLFP,2)-1):params.totalTrialLength,mean(alphaLFP(redTrials,:),1),'Color','r');
hold on;
plot(0:params.totalTrialLength/(size(alphaLFP,2)-1):params.totalTrialLength,mean(alphaLFP(blueTrials,:),1),'Color','b');
minValForPlot=min(min([mean(alphaLFP(redTrials,:),1); mean(alphaLFP(blueTrials,:),1)],[],2),[],1);
binSize=0.15;
for i=binSize:binSize:params.totalTrialLength
    thisRedSeg=alphaLFP(redTrials,floor((i-binSize)/(params.totalTrialLength/(size(alphaLFP,2)-1)))+1:floor(i/(params.totalTrialLength/(size(alphaLFP,2)-1)))+1);
    thisRedSeg(thisRedSeg<0)=-thisRedSeg(thisRedSeg<0);
    thisRedSeg=sum(thisRedSeg,2);
    thisBlueSeg=alphaLFP(blueTrials,floor((i-binSize)/(params.totalTrialLength/(size(alphaLFP,2)-1)))+1:floor(i/(params.totalTrialLength/(size(alphaLFP,2)-1)))+1);
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
title('Alpha, Red v. Blue, Black=Significance');

