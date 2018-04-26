function [x,av,lowBound,highBound]=plotBootstrappedPowerSpec(trialSpecgrams,Nsamples,sampling_rate,freq,freqBand,bootstrapError)

X_div=0.5; % in seconds

powerSpecs=zeros(length(trialSpecgrams.specgrams),size(trialSpecgrams.specgrams{1},2));
for i=1:length(trialSpecgrams.specgrams)
    temp=trialSpecgrams.specgrams{i};
    powerSpecs(i,:)=mean(temp(freq>=freqBand(1) & freq<=freqBand(2),:),1);
end

figure();
x=linspace(0,Nsamples/sampling_rate,size(trialSpecgrams.specgrams{1},2));
av=mean(powerSpecs,1);
plot(x,av,'Color','k');
hold on; 
if bootstrapError==0
    % Get 95% confidence intervals
    lowBound=zeros(1,size(powerSpecs,2));
    highBound=zeros(1,size(powerSpecs,2));
    for i=1:size(powerSpecs,2)
        currSet=powerSpecs(:,i);
        rankedCurrSet=sort(currSet);
        n95=floor(0.95*length(currSet));
        highBound(i)=rankedCurrSet(n95);
        lowBound(i)=rankedCurrSet(end-n95);
    end
    plot(x,lowBound,'Color',[0.5 0.5 0.5]);
    plot(x,highBound,'Color',[0.5 0.5 0.5]);
else
    nBootTrials=100;
    fractionToUse=0.8;
    nToUse=floor(size(powerSpecs,1)*fractionToUse);
    bootPowerSpecs=zeros(nBootTrials,size(powerSpecs,2));
    for i=1:nBootTrials
        avNow=mean(powerSpecs(randsample(size(powerSpecs,1),nToUse),:),1);
        bootPowerSpecs(i,:)=avNow;
    end
    lowBound=zeros(1,size(bootPowerSpecs,2));
    highBound=zeros(1,size(bootPowerSpecs,2));
    for i=1:size(bootPowerSpecs,2)
        currSet=bootPowerSpecs(:,i);
        rankedCurrSet=sort(currSet);
        n95=floor(0.95*length(currSet));
        highBound(i)=rankedCurrSet(n95);
        lowBound(i)=rankedCurrSet(end-n95);
    end
    plot(x,lowBound,'Color',[0.5 0.5 0.5]);
    plot(x,highBound,'Color',[0.5 0.5 0.5]);
end