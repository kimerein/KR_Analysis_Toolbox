function calculateStimCondResponseDistributions(data,ledForSweeps,stimsForSweeps,LFP_Fs,baselineWindow)
% Calculate the different distributions of data for each stimulus/LED
% condition; ideal observer analysis of trial-by-trial information about
% stimulus/LED condition
% 
% PARAMETERS:
% data: matrix with samples as columns and trials as rows; this is data
% from which response distributions will be drawn
% Each row is a response, and specific response windows are defined below
% Average over each response window for each trial, and this produces a
% distribution of values for each stimulus/LED condition
% Plot the distributions for all stimulus/LED conditions
% Also, for two specific distributions, compare an ideal observer's ability
% to discriminate stimulus parameter on the basis of trial-by-trial
% observations
% ledForSweeps: vector with 1's for trials with LED, 0's for trials without
% LED; rows correspond to trial rows of data
% stimsForSweeps: vector with each stimulus condition for all trials; rows
% correspond to trial rows of data
% LFP_Fs: sampling rate of data
% baselineWindow: time window relative to trial beginning to use to
% calculate baseline

data=data(2:end-1,:);
ledForSweeps=ledForSweeps(2:end-1);
stimsForSweeps=stimsForSweeps(2:end-1);

% Need to collapse stimsForSweeps over phase
spatFreq1_inds=find(stimsForSweeps==1|stimsForSweeps==2|stimsForSweeps==3|stimsForSweeps==4);
spatFreq2_inds=find(stimsForSweeps==5|stimsForSweeps==6|stimsForSweeps==7|stimsForSweeps==8);
spatFreq3_inds=find(stimsForSweeps==9|stimsForSweeps==10|stimsForSweeps==11|stimsForSweeps==12);
spatFreq4_inds=find(stimsForSweeps==13|stimsForSweeps==14|stimsForSweeps==15|stimsForSweeps==16);
spatFreq5_inds=find(stimsForSweeps==17|stimsForSweeps==18|stimsForSweeps==19|stimsForSweeps==20);
stimsForSweeps(spatFreq1_inds)=1;
stimsForSweeps(spatFreq2_inds)=2;
stimsForSweeps(spatFreq3_inds)=3;
stimsForSweeps(spatFreq4_inds)=4;
stimsForSweeps(spatFreq5_inds)=5;

ONresponseWindow=[0.4 1];
% afterResponseWindow=[1.5 1.9];
afterResponseWindow=[1.6 1.9];

totalTrialLength=(1/LFP_Fs)*size(data,2);
times=0:totalTrialLength/(size(data,2)-1):totalTrialLength;
baselineInds=floor(baselineWindow(1)*LFP_Fs)+1:floor(baselineWindow(2)*LFP_Fs);

ONInds=floor(ONresponseWindow(1)*LFP_Fs)+1:floor(ONresponseWindow(2)*LFP_Fs);
afterInds=floor(afterResponseWindow(1)*LFP_Fs)+1:floor(afterResponseWindow(2)*LFP_Fs);

% Align all trials to baseline
data=data-mean(data(:,baselineInds),2)*ones(1,size(data,2));

allStims=sort(unique(stimsForSweeps));

LEDon_ONdistributions=cell(1,length(allStims));
LEDoff_ONdistributions=cell(1,length(allStims));
LEDon_afterdistributions=cell(1,length(allStims));
LEDoff_afterdistributions=cell(1,length(allStims));

for i=1:length(allStims)
    currSet=stimsForSweeps==allStims(i)&ledForSweeps>0;
    LEDon_ONdistributions{i}=mean(data(currSet,ONInds),2);
    currSet=stimsForSweeps==allStims(i)&ledForSweeps==0;
    LEDoff_ONdistributions{i}=mean(data(currSet,ONInds),2);
    currSet=stimsForSweeps==allStims(i)&ledForSweeps>0;
    LEDon_afterdistributions{i}=mean(data(currSet,afterInds),2);
    currSet=stimsForSweeps==allStims(i)&ledForSweeps==0;
    LEDoff_afterdistributions{i}=mean(data(currSet,afterInds),2);
end

temp=LEDon_ONdistributions{3};
temp=temp(temp>-0.02);
LEDon_ONdistributions{3}=temp;
temp=LEDoff_ONdistributions{3};
temp=temp(temp>-0.02);
LEDoff_ONdistributions{3}=temp;
temp=LEDon_afterdistributions{3};
temp=temp(temp>-0.02);
LEDon_afterdistributions{3}=temp;
temp=LEDoff_afterdistributions{3};
temp=temp(temp>-0.02);
LEDoff_afterdistributions{3}=temp;

% temp=LEDon_ONdistributions{3};
% temp=temp(temp>-0.5);
% LEDon_ONdistributions{3}=temp;
% temp=LEDoff_ONdistributions{3};
% temp=temp(temp>-0.5);
% LEDoff_ONdistributions{3}=temp;
% temp=LEDon_afterdistributions{3};
% temp=temp(temp>-0.5);
% LEDon_afterdistributions{3}=temp;
% temp=LEDoff_afterdistributions{3};
% temp=temp(temp>-0.5);
% LEDoff_afterdistributions{3}=temp;

% Plot ON + LED response distributions
figure;
for i=1:length(allStims)
    if i==2|i==3|i==4
        continue;
    end
    [n,xout] = histnorm(LEDon_ONdistributions{i},10);
    plot(xout,n);
    hold all;
end
legend;
% Plot ON + LED off response distributions
figure;
for i=1:length(allStims)
    if i==2|i==3|i==4
        continue;
    end
    [n,xout] = histnorm(LEDoff_ONdistributions{i},10);
    plot(xout,n);
    hold all;
end
legend;
% Plot after + LED response distributions
figure;
for i=1:length(allStims)
    [n,xout] = histnorm(LEDon_afterdistributions{i},10);
    plot(xout,n);
    hold all;
end
legend;
% Plot after + LED off response distributions
figure;
for i=1:length(allStims)
    [n,xout] = histnorm(LEDoff_afterdistributions{i},10);
    plot(xout,n);
    hold all;
end
legend;

disp(length(LEDoff_ONdistributions{1})+length(LEDoff_ONdistributions{5}));
disp(length(LEDoff_ONdistributions{1}));
disp(length(LEDoff_ONdistributions{5}));
nTotal=length(LEDoff_ONdistributions{1})+length(LEDoff_ONdistributions{5});
n1=length(LEDoff_ONdistributions{1});
n2=length(LEDoff_ONdistributions{5});

% Calc. p-val
dist1=LEDon_ONdistributions{1};
dist2=LEDon_ONdistributions{5};


tryThresh=-0.05:0.001:0.05;
bestRatio=0;
bestFractionCorrect=0;
bestFractionIncorrect=0;
bestThresh=0;
for i=1:length(tryThresh)
    count1=sum(dist1>tryThresh(i));
    count2=sum(dist2<=tryThresh(i));
    correctSum=count1+count2;
    incorrectSum=sum(dist1<=tryThresh(i))+sum(dist2>tryThresh(i));
    if correctSum/incorrectSum>bestRatio
        bestRatio=correctSum/incorrectSum;
        bestFractionCorrect=correctSum/(correctSum+incorrectSum);
        bestFractionIncorrect=incorrectSum/(correctSum+incorrectSum);
        bestThresh=tryThresh(i);
    end
end
disp(bestRatio);
disp('bestFractionCorrect');
disp(bestFractionCorrect);
disp(bestFractionIncorrect);
disp(bestThresh);

n1=length(dist1);
n2=length(dist2);
nTotal=n1+n2;
obs=[zeros(1,n1) ones(1,n2)];
obs=obs(randperm(nTotal));
fractionRight=zeros(1,10000);
for i=1:10000
    guesses=zeros(1,nTotal);
    for j=1:nTotal
        guess=rand(1,1);
        if guess>0.5
            guesses(j)=1;
        elseif guess<0.5
            guesses(j)=0;
        else
            j=j-1;
        end
    end
    fractionRight(i)=(sum(obs&logical(guesses))+sum(~obs&~logical(guesses)))/nTotal;
end
%disp(sum(fractionRight>0.5705)/10000);
disp(sum(fractionRight>0.6282));
disp(sum(fractionRight>0.6282)/10000);
 