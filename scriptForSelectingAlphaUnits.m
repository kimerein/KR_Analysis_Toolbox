function [t,f,y,outS,nCriterionUnits]=scriptForSelectingAlphaUnits(noTheta_trialAv_noLED,alphaThresh,dd)

forShort=1;

% spontWindowEnd=3.9;
% spontWindowStartAgain=10;
% spontWindowEnd=0.9;
% spontWindowStartAgain=4.5;
spontWindowEnd=1;
spontWindowStartAgain=5;
% spontWindowEnd=1.2;
% spontWindowStartAgain=2.5;

if ~isempty(dd)
    a=load([dd '\' 'UbyUnoTheta_trialAv_noLED']);
    noTheta_noLED=a.noTheta_trialAv_noLED;
    a=load([dd '\' 'UbyUnoTheta_trialAv_LED']);
    noTheta_LED=a.noTheta_trialAv_LED;
    for i=1:length(noTheta_noLED.low.S)
        FORALPHAnoTheta_trialAv_noLED.low.S{i}=(noTheta_noLED.low.S{i}+noTheta_LED.low.S{i})./2;
        FORALPHAnoTheta_trialAv_noLED.high.S{i}=(noTheta_noLED.high.S{i}+noTheta_LED.high.S{i})./2;
    end
else
    FORALPHAnoTheta_trialAv_noLED=noTheta_trialAv_noLED;
end

t=noTheta_trialAv_noLED.t;
f=noTheta_trialAv_noLED.f;
for i=1:(length(FORALPHAnoTheta_trialAv_noLED.low.S))
slow=FORALPHAnoTheta_trialAv_noLED.low.S{i};
shigh=FORALPHAnoTheta_trialAv_noLED.high.S{i};
% spontAlpha(i,:)=nanmean(slow(t<3.9 | t>10,:),1)'+nanmean(shigh(t<3.9 | t>10,:),1)';
% spontAlpha(i,:)=nanmean(slow(t<spontWindowEnd | t>spontWindowStartAgain,:),1)';
spontAlpha(i,:)=nanmean(slow(t<spontWindowEnd | t>spontWindowStartAgain,:)+shigh(t<spontWindowEnd | t>spontWindowStartAgain,:),1)';
end

alphaRatio=zeros(size(spontAlpha,1),1);
for i=1:size(spontAlpha,1)
    if forShort==1
%         alphaRatio(i)=nanmean(spontAlpha(i,f>=6 & f<=20),2)./nanmean(spontAlpha(i,f>20 & f<30),2);
        alphaRatio(i)=nanmean(spontAlpha(i,f>=6 & f<=20),2)./nanmean(spontAlpha(i,f<6 | f>20),2);
    else
        alphaRatio(i)=nanmean(spontAlpha(i,f>=6 & f<=20),2)./nanmean(spontAlpha(i,f<6 | f>20),2);
    end
end
figure(); plot(f,spontAlpha(alphaRatio>alphaThresh,:)','Color','r');
hold on;
if any(alphaRatio<=alphaThresh)
    plot(f,spontAlpha(alphaRatio<=alphaThresh,:)','Color','k');
end
figure(); hist(alphaRatio,30);

integrates=find(alphaRatio>alphaThresh);
runningSum=zeros(size(noTheta_trialAv_noLED.low.S{1}));
for i=1:length(integrates)
% runningSum=runningSum+noTheta_trialAv_noLED.low.S{integrates(i)}+noTheta_trialAv_noLED.high.S{integrates(i)};
if isnan(noTheta_trialAv_noLED.low.S{integrates(i)})
    continue
else
    runningSum=runningSum+noTheta_trialAv_noLED.low.S{integrates(i)};
end
end
figure(); 
% plot(t,nanmean(runningSum(:,f>11.5 & f<=12.5),2));
plot(t,nanmean(runningSum(:,f>11 & f<=20),2));
% y=nanmean(runningSum(:,f>11.5 & f<=12.5),2);
y=nanmean(runningSum(:,f>11 & f<=20),2);
nCriterionUnits=length(integrates);
figure(); imagesc(t,f,runningSum');
outS=runningSum;


end