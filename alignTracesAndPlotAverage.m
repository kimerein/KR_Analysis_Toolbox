function [amps,alignedResps,normAmps]=alignTracesAndPlotAverage(x,traces,scaleVals,sineFreq)

% pulseTimes=[1 1.221 1.442 1.663 1.884 2.106 2.327 2.548 2.769 2.99]; % 10x = 4.5 Hz
% alignBefore=0.1; % in s
alignBefore=0.05; % in s
earliestOnset=0.02; % in s after pulse
baseBefore=[0.025 0.027];
% takeAfter=0.1; % in s
takeAfter=0.05; % in s
firstTroughWindow=0.05; % in s
band=[1 10000];
pulseFreq=1;
% sineFreq=10;
timeOffset=1;
isSpikes=1;

if isSpikes==1
    baseBefore=[0 0.027];
end

if pulseFreq==1
    StartPt=1;
    EndPt=3;
    V=[0.05 0.06:0.001:0.1 4.9 4.91 4.92 4.93 4.94 4.95 4.96 4.97 4.98];
    I=[1.2 7.4 33.8 59.6 84.6 108.6 131.5 153.7 175 196 217 237 256 275 293 311 329 346 363 380 398 414 431 447 464 480 496 511 528 545 561 577 590 606 622 639 650 665 682 697 704 718 728 745 759 772 780 787 810 824 838];
    LEDrange=I(end)-I(1);
    useLEDrange=LEDrange;
    offse=(LEDrange/2)+I(1);
    a=useLEDrange/2;
    herex=linspace(0,EndPt-StartPt,(EndPt-StartPt)*1000);
    y=a.*sin(2*pi*sineFreq.*herex)+offse;
    y_volt=interp1(I,V,y);
    onTimes=y_volt>0.2;
    pulseTimes=[];
    inPulse=0;
    for i=1:length(onTimes)
        if onTimes(i)==1
            if inPulse==1
            else
                pulseTimes=[pulseTimes herex(i)+timeOffset];
                inPulse=1;
            end
        else
            if inPulse==1
                inPulse=0;
            else
            end
        end
    end
%     disp(pulseTimes);    
else
    StartPt=1;
    EndPt=3;
    % StartPt=1;
    % EndPt=11;
    indPulseWidth=0.01;
    nPulses=10;
    totalPulsesTime=indPulseWidth*nPulses;
    remainingTime=(EndPt-StartPt)-totalPulsesTime;
    indWidthOfIPI=remainingTime/(nPulses-1);
    pulseTimes=zeros(1,nPulses);
    for i=1:nPulses
        pulseTimes(i)=StartPt+(i-1)*(indPulseWidth+indWidthOfIPI);
    end
%     disp(pulseTimes);
end

if isSpikes~=1
    for i=1:size(traces,1)
        traces(i,:)=bandPassLFP(traces(i,:),1/(x(2)-x(1)),band(1),band(2),0);
    end
end

alignedResps=cell(1,length(pulseTimes));
amps=zeros(size(traces,1),length(pulseTimes));
takeXs=cell(1,length(pulseTimes));
for i=1:length(pulseTimes)
    takeX=x(x>=pulseTimes(i)-alignBefore & x<=pulseTimes(i)+takeAfter);
    takeY=traces(:,x>=pulseTimes(i)-alignBefore & x<=pulseTimes(i)+takeAfter);
    takeXs{i}=takeX;
    % Align by before
    for j=1:size(takeY,1)
        currBase=mean(takeY(j,takeX>=pulseTimes(i)+baseBefore(1) & takeX<=pulseTimes(i)+baseBefore(2)));
        takeY(j,:)=takeY(j,:)-currBase;
    end
    subTakeY=takeY(:,takeX>=pulseTimes(i)+earliestOnset & takeX<=pulseTimes(i)+firstTroughWindow);
    if isSpikes==1
        amps(:,i)=max(subTakeY,[],2);
    else
        amps(:,i)=-min(subTakeY,[],2);
    end
    alignedResps{i}=takeY;
end

figure();
for i=1:length(pulseTimes)
    if isSpikes==1
        line([pulseTimes(i) pulseTimes(i)],[min(alignedResps{i}) max(alignedResps{i})],'Color','k');
    else
        line([pulseTimes(i) pulseTimes(i)],[-2 0.5],'Color','k');
    end
    hold on;
    plot(takeXs{i},alignedResps{i});
end

normAmps=zeros(size(amps));
for i=1:size(amps,1)
    if isempty(scaleVals)
        resp1Amp=amps(i,1);
    else
        resp1Amp=scaleVals(i);
    end
    scale=1/resp1Amp;
    normAmps(i,:)=amps(i,:).*scale;
end
% figure(); 
% errorbar(1:length(pulseTimes),mean(normAmps,1),std(normAmps,[],1)./sqrt(size(normAmps,1)));