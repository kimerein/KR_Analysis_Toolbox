function makeUnitResponseHeatmap(spikes,fileInd,includeStimVals,binsize,duration)

ledOffVal=[];
ledOnVal=[];
norm=1;
chopSpikesAccordingToPhase=1;
stimDuration=[1 3]; % in s
onsetPeakWindow=0.2; % in s
offsetPeakWindow=0.2; % in s
frequencyWindow=0.5; % in s
freqRangeForPeak=[3 10]; % in Hz

spikes=filtspikes(spikes,0,'fileInd',fileInd);
a=unique(spikes.assigns);
% s=unique(spikes.stimcond(~isnan(spikes.stimcond)));
s=includeStimVals;

nTrials=zeros(length(s),2);
theseTrials=cell(length(s),2);
for i=1:length(s)
    subS=filtspikes(spikes,0,'stimcond',s(i),'fileInd',fileInd);
    if ~isempty(ledOffVal)
        subS1=filtspikes(subS,0,'led',ledOffVal);
    else
        subS1=subS;
    end
    nTrials(i,1)=length(unique(subS1.trials));
    theseTrials{i,1}=unique(subS1.trials);
    if ~isempty(ledOnVal)
        subS1=filtspikes(subS,0,'led',ledOnVal);
    else
        subS1=subS;
    end
    nTrials(i,2)=length(unique(subS1.trials));
    theseTrials{i,2}=unique(subS1.trials);
end    
if ~isempty(ledOffVal)
    subS=filtspikes(spikes,0,'led',ledOffVal);
else
    subS=spikes;
end
temp{1}=unique(subS.trials);
if ~isempty(ledOnVal)
    subS=filtspikes(spikes,0,'led',ledOnVal);
else
    subS=spikes;
end
temp{2}=unique(subS.trials);
[x,y1,y2]=makePSTH(spikes,fileInd,ledOffVal,ledOnVal,binsize,sum(nTrials,1),duration,temp);
figure();
plot(x,y1,'Color','k');
hold on;
plot(x,y2,'Color','r');
allSpikesx=x;

a=unique(spikes.assigns);
ledOff=nan(length(a),length(s),length(y1));
ledOn=nan(length(a),length(s),length(y1));
ledBoth=nan(length(a),length(s),length(y1));
allSpikes=nan(length(s),length(y1));
for i=1:length(s)
    temp{1}=[theseTrials{i,1} theseTrials{i,2}];
    temp{2}=theseTrials{i,2};
    [x,y1,y2]=makePSTH(filtspikes(spikes,0,'stimcond',s(i)),fileInd,ledOffVal,[ledOffVal ledOnVal],binsize,[nTrials(i,2) sum(nTrials(i,:),2)],duration,temp);
    allSpikes(i,:)=y1;
end
for i=1:length(a)
    disp(i);
    for j=1:length(s)
        temp{1}=theseTrials{j,1};
        temp{2}=theseTrials{j,2};
        [x,y1,y2]=makePSTH(filtspikes(spikes,0,'assigns',a(i),'stimcond',s(j)),fileInd,ledOffVal,ledOnVal,binsize,nTrials(j,:),duration,temp);
        if isempty(y1) || isempty(y2)
            continue
        end
        ledOff(i,j,:)=y1;
        ledOn(i,j,:)=y2;
        temp{1}=[theseTrials{j,1} theseTrials{j,2}];
        temp{2}=[theseTrials{j,1} theseTrials{j,2}];
        [x,y1,y2]=makePSTH(filtspikes(spikes,0,'assigns',a(i),'stimcond',s(j)),fileInd,[ledOffVal ledOnVal],[ledOffVal ledOnVal],binsize,[sum(nTrials(j,:),2) sum(nTrials(j,:),2)],duration,temp);
        ledBoth(i,j,:)=y1;
    end
end
PSTHs.ledOff=ledOff;
PSTHs.ledOn=ledOn;
PSTHs.ledBoth=ledBoth;

sepPSTHbyStim=cell(1,length(s));
for i=1:length(s)
    currStim=reshape(ledOff(:,i,:),length(a),size(ledOff,3));
    sepPSTHbyStim{i}=currStim;
    idx1=kmeans(currStim,2);
    [s1,i1]=sort(idx1);
    currStim=currStim(i1,:);
    if norm==1
        m1=max(currStim,[],2);
        currStim(m1~=0,:)=currStim(m1~=0,:)./repmat(m1(m1~=0),1,size(currStim,2));
    end
    figure();
    imagesc(currStim);
    title(['currStim' num2str(i)]);
end

togetherSepPSTH=[];
topRow=[];
for i=1:length(s)
    topRow=[topRow ones(1,size(sepPSTHbyStim{i},2)).*i];
    togetherSepPSTH=[togetherSepPSTH sepPSTHbyStim{i}];
end
idx2=kmeans(togetherSepPSTH,2);
[s2,i2]=sort(idx2);
togetherSepPSTH=togetherSepPSTH(i2,:);
if norm==1
    m1=max(togetherSepPSTH,[],2);
    togetherSepPSTH(m1~=0,:)=togetherSepPSTH(m1~=0,:)./repmat(m1(m1~=0),1,size(togetherSepPSTH,2));
end
figure();
imagesc([topRow./max(topRow); togetherSepPSTH]);

if chopSpikesAccordingToPhase==1
    % Find phase triggered on alpha onset at stim. onset
    % Find phase triggered on alpha at stim. offset -- inverse time
    onsetWindow=[stimDuration(1) stimDuration(1)+onsetPeakWindow];
    pre_offsetWindow=[stimDuration(2)-offsetPeakWindow stimDuration(2)];
    peri_offsetWindow=[stimDuration(2)-offsetPeakWindow/2 stimDuration(2)+offsetPeakWindow/2];
    post_offsetWindow=[stimDuration(2) stimDuration(2)+offsetPeakWindow];
    % Get spectrograms for different stimuli
    params.Fs=1/(allSpikesx(2)-allSpikesx(1));
    params.tapers=[0.5 allSpikesx(end)-allSpikesx(1) 0];
    allSp=cell(1,length(s));
    allf=cell(1,length(s));
    allt=cell(1,length(s));
    for i=1:length(s)
        [Sp,t,f]=mtspecgrampb(allSpikes(i,:),[0.5 0.1],params);
        allSp{i}=Sp;
        allf{i}=f;
        allt{i}=t;
    end    
    % Find peak
    % Find frequency around peak
    onsetPeaks=zeros(1,length(s));
    pre_offsetPeaks=zeros(1,length(s));
    peri_offsetPeaks=zeros(1,length(s));
    post_offsetPeaks=zeros(1,length(s));
    onsetFreqs=zeros(1,length(s));
    pre_offsetFreqs=zeros(1,length(s));
    peri_offsetFreqs=zeros(1,length(s));
    post_offsetFreqs=zeros(1,length(s));
    for i=1:length(s)
        currWindow=onsetWindow;
        [onsetPeaks(i),onsetFreqs(i)]=getPeakAndFreq(currWindow,allSpikesx,allSpikes(i,:),frequencyWindow,allSp{i},allf{i},allt{i},freqRangeForPeak);
        currWindow=pre_offsetWindow;
        [pre_offsetPeaks(i),pre_offsetFreqs(i)]=getPeakAndFreq(currWindow,allSpikesx,allSpikes(i,:),frequencyWindow,allSp{i},allf{i},allt{i},freqRangeForPeak);
        currWindow=peri_offsetWindow;
        [peri_offsetPeaks(i),peri_offsetFreqs(i)]=getPeakAndFreq(currWindow,allSpikesx,allSpikes(i,:),frequencyWindow,allSp{i},allf{i},allt{i},freqRangeForPeak);
        currWindow=post_offsetWindow;
        [post_offsetPeaks(i),post_offsetFreqs(i)]=getPeakAndFreq(currWindow,allSpikesx,allSpikes(i,:),frequencyWindow,allSp{i},allf{i},allt{i},freqRangeForPeak);
    end    
    timebinsSet=cell(1,4);
    for i=1:length(s)
        % Get different phase window sets
        timebinsSet{1}=getPhaseWindowSets(onsetFreqs(i),onsetPeaks(i),duration);
        timebinsSet{2}=getPhaseWindowSets(pre_offsetFreqs(i),pre_offsetPeaks(i),duration);
        timebinsSet{3}=getPhaseWindowSets(peri_offsetFreqs(i),peri_offsetPeaks(i),duration);
        timebinsSet{4}=getPhaseWindowSets(post_offsetFreqs(i),post_offsetPeaks(i),duration);
        % Label spikes by different phase window sets
        spikes=chopSpikesByPhase(spikes,timebinsSet,s(i));
    end
    
    % Plot tuning of different phase sets
    for i=1:length(timebinsSet)
       figure();
       tuning=zeros(1,length(s));
       tuningDuringStim=zeros(1,length(s));
       tuningAfterStim=zeros(1,length(s));
       for j=1:length(s)
           [px,py1]=makePSTHForMUA(filtspikes(spikes,0,'stimcond',s(j),['set' num2str(i)],i),[],ledOffVal,ledOnVal,50,duration);
           plot(px,py1);
           tuning(j)=mean(py1);
           tuningDuringStim(j)=mean(py1(px>=stimDuration(1) & px<=stimDuration(2)));
           tuningDuringStim(j)=mean(py1(px>=3 & px<=4));
           tuningAfterStim(j)=mean(py1(px>=stimDuration(2)));
           hold all;
       end
       title(['set' num2str(i)]);
       if i==1 || i==2 || i==3 || i==4
        figure();
        plot(s,tuningAfterStim);
        title(['tuning set after timebins' num2str(i)]);
        figure();
        plot(s,tuningDuringStim);
        title(['tuning set during timebins' num2str(i)]);
       end
    end
    
    % Plot different phase sets for each stim.
    for j=1:length(s)
       figure();
       for i=1:length(timebinsSet)
           if ismember(i,[2 3]) % exclude which offset time bins
               continue
           end
           [px,py1]=makePSTHForMUA(filtspikes(spikes,0,'stimcond',s(j),['set' num2str(i)],i),[],ledOffVal,ledOnVal,50,duration);
           if i==1
               c='k';
           else
               c='r';
           end
           plot(px,py1,'Color',c);
           hold on;
       end
       title(['stimcond' num2str(j)]);
    end
    
end

end

function timebinsSet=getPhaseWindowSets(onsetFreqs,onsetPeaks,duration)
    currPeriod=1/onsetFreqs;
    starts=[fliplr(onsetPeaks-currPeriod/4:-currPeriod:0) onsetPeaks-currPeriod/4+currPeriod:currPeriod:duration];
    ends=[fliplr(onsetPeaks+currPeriod/4:-currPeriod:0) onsetPeaks+currPeriod/4+currPeriod:currPeriod:duration];
    ends=ends(ends>=starts(1));
    starts=starts(starts<=ends(end));
    timebinsSet=[starts' ends'];   
end

function [peakTime,peakFreq]=getPeakAndFreq(currWindow,x,allSpikes,frequencyWindow,sp,f,t,freqRange)
% Get peak
useInds=x>=currWindow(1) & x<=currWindow(2);
[~,peakInd]=max(allSpikes(useInds));
peakTime=x(find(useInds,1,'first')+(peakInd-1));
% Get freq
useInds=t>=mean(currWindow)-frequencyWindow/2 & t<=mean(currWindow)+frequencyWindow/2;
useFreqInds=f>=freqRange(1) & f<=freqRange(2);
[~,peakFreqInd]=max(mean(sp(useInds,useFreqInds),1),[],2);    
peakFreq=f(find(useFreqInds,1,'first')+(peakFreqInd-1));
end

function spikes=chopSpikesByPhase(spikes,timebinsSet,useStimcond)
for j=1:length(timebinsSet)
    if ~isfield(spikes,['set' num2str(j)])
        spikes.(['set' num2str(j)])=zeros(size(spikes.spiketimes));
    end
    timebinsSet1=timebinsSet{j};
    for i=1:size(timebinsSet1,1)
        currsub=spikes.(['set' num2str(j)]);
        currsub(spikes.stimcond==useStimcond & spikes.spiketimes>=timebinsSet1(i,1) & spikes.spiketimes<=timebinsSet1(i,2))=j;
        spikes.(['set' num2str(j)])=currsub;
    end
end
end

function [xpoints1,ypoints1,ypoints2]=makePSTH(spikes,useTheseFileInds,noLedValues,ledValue,binsize,nTrials,duration,theseTrials)

useBlackTrials=[];
useRedTrials=[];
% binsize=5; % ms

stimconds=[1:128];

nunits=1;

spikes=filtspikes(spikes,0,'stimcond',stimconds); 

if isempty(ledValue)
    allUnits_withLED=filtspikes(spikes,0,'fileInd',useTheseFileInds);
    allUnits_noLED=filtspikes(spikes,0,'fileInd',useTheseFileInds);
else
    if ~isempty(useTheseFileInds)
        temp=[];
        temp1=[];
        for i=1:length(ledValue)
            spikes=makeTempField(spikes,'led',ledValue(i));
            temp(i,:)=spikes.temp;
            temp1(i,:)=spikes.sweeps.temp;
            spikes.temp=[];
            spikes.sweeps.temp=[];
        end
        spikes.temp=sum(temp,1)>=1;
        spikes.sweeps.temp=sum(temp1,1)>=1;
        allUnits_withLED=filtspikes(spikes,0,'temp',1,'fileInd',useTheseFileInds);
        temp=[];
        temp1=[];
        for i=1:length(noLedValues)
            spikes=makeTempField(spikes,'led',noLedValues(i));
            temp(i,:)=spikes.temp;
            temp1(i,:)=spikes.sweeps.temp;
            spikes.temp=[];
            spikes.sweeps.temp=[];
        end
        spikes.temp=sum(temp,1)>=1;
        spikes.sweeps.temp=sum(temp1,1)>=1;
        allUnits_noLED=filtspikes(spikes,0,'temp',1,'fileInd',useTheseFileInds);
    else
        temp=[];
        temp1=[];
        for i=1:length(ledValue)
            spikes=makeTempField(spikes,'led',ledValue(i));
            temp(i,:)=spikes.temp;
            temp1(i,:)=spikes.sweeps.temp;
            spikes.temp=[];
            spikes.sweeps.temp=[];
        end
        spikes.temp=sum(temp,1)>=1;
        spikes.sweeps.temp=sum(temp1,1)>=1;
        allUnits_withLED=filtspikes(spikes,0,'temp',1);
        temp=[];
        temp1=[];
        for i=1:length(noLedValues)
            spikes=makeTempField(spikes,'led',noLedValues(i));
            temp(i,:)=spikes.temp;
            temp1(i,:)=spikes.sweeps.temp;
            spikes.temp=[];
            spikes.sweeps.temp=[];
        end
        spikes.temp=sum(temp,1)>=1;
        spikes.sweeps.temp=sum(temp1,1)>=1;
        allUnits_noLED=filtspikes(spikes,0,'temp',1);
    end
end

allSpiketimes_withLED=[];
allSpiketimes_noLED=[];

% if ~isempty(useBlackTrials)
%     allUnits_noLED=filtspikes(allUnits_noLED,0,'trials',sort(useBlackTrials));
% end
% if ~isempty(useRedTrials)
%     allUnits_withLED=filtspikes(allUnits_withLED,0,'trials',sort(useRedTrials));
% end
%     
[n1,centers1,edges1,xpoints1,ypoints1,stds1]=psth_wStd(allUnits_noLED,binsize,0,duration,nTrials(1),theseTrials{1});

[n2,centers2,edges2,xpoints2,ypoints2,stds2]=psth_wStd(allUnits_withLED,binsize,0,duration,nTrials(2),theseTrials{2});

ypoints1=ypoints1/nunits;
ypoints2=ypoints2/nunits;
end

function [xpoints1,ypoints1,ypoints2]=makePSTHForMUA(spikes,useTheseFileInds,noLedValues,ledValue,binsize,duration)

useBlackTrials=[];
useRedTrials=[];
% binsize=5; % ms

stimconds=[1:128];

nunits=1;

spikes=filtspikes(spikes,0,'stimcond',stimconds); 

if isempty(ledValue)
    allUnits_withLED=filtspikes(spikes,0,'fileInd',useTheseFileInds);
    allUnits_noLED=filtspikes(spikes,0,'fileInd',useTheseFileInds);
else
    if ~isempty(useTheseFileInds)
        temp=[];
        temp1=[];
        for i=1:length(ledValue)
            spikes=makeTempField(spikes,'led',ledValue(i));
            temp(i,:)=spikes.temp;
            temp1(i,:)=spikes.sweeps.temp;
            spikes.temp=[];
            spikes.sweeps.temp=[];
        end
        spikes.temp=sum(temp,1)>=1;
        spikes.sweeps.temp=sum(temp1,1)>=1;
        allUnits_withLED=filtspikes(spikes,0,'temp',1,'fileInd',useTheseFileInds);
        temp=[];
        temp1=[];
        for i=1:length(noLedValues)
            spikes=makeTempField(spikes,'led',noLedValues(i));
            temp(i,:)=spikes.temp;
            temp1(i,:)=spikes.sweeps.temp;
            spikes.temp=[];
            spikes.sweeps.temp=[];
        end
        spikes.temp=sum(temp,1)>=1;
        spikes.sweeps.temp=sum(temp1,1)>=1;
        allUnits_noLED=filtspikes(spikes,0,'temp',1,'fileInd',useTheseFileInds);
    else
        temp=[];
        temp1=[];
        for i=1:length(ledValue)
            spikes=makeTempField(spikes,'led',ledValue(i));
            temp(i,:)=spikes.temp;
            temp1(i,:)=spikes.sweeps.temp;
            spikes.temp=[];
            spikes.sweeps.temp=[];
        end
        spikes.temp=sum(temp,1)>=1;
        spikes.sweeps.temp=sum(temp1,1)>=1;
        allUnits_withLED=filtspikes(spikes,0,'temp',1);
        temp=[];
        temp1=[];
        for i=1:length(noLedValues)
            spikes=makeTempField(spikes,'led',noLedValues(i));
            temp(i,:)=spikes.temp;
            temp1(i,:)=spikes.sweeps.temp;
            spikes.temp=[];
            spikes.sweeps.temp=[];
        end
        spikes.temp=sum(temp,1)>=1;
        spikes.sweeps.temp=sum(temp1,1)>=1;
        allUnits_noLED=filtspikes(spikes,0,'temp',1);
    end
end

allSpiketimes_withLED=[];
allSpiketimes_noLED=[];

% if ~isempty(useBlackTrials)
%     allUnits_noLED=filtspikes(allUnits_noLED,0,'trials',sort(useBlackTrials));
% end
% if ~isempty(useRedTrials)
%     allUnits_withLED=filtspikes(allUnits_withLED,0,'trials',sort(useRedTrials));
% end
%     
allt=unique(allUnits_noLED.trials);
[n1,centers1,edges1,xpoints1,ypoints1,stds1]=psth_wStd(allUnits_noLED,binsize,0,duration,length(allt),allt);
allt=unique(allUnits_withLED.trials);
[n2,centers2,edges2,xpoints2,ypoints2,stds2]=psth_wStd(allUnits_withLED,binsize,0,duration,length(allt),allt);

ypoints1=ypoints1/nunits;
ypoints2=ypoints2/nunits;
end

function [varargout]=psth_wStd(spikes,binsize,bsmooth,duration,nTrials,theseTrials)

if nargin < 2
    binsize = 50; 
end
if nargin < 3
    bsmooth = 1;
end
% Set duration and number of trials
if ~isempty(nTrials)
    numtrials=nTrials;
elseif isfield(spikes,'sweeps')
    a=unique(spikes.trials);
    if length(spikes.sweeps.trials)==4*length(a)
        numtrials=length(unique(spikes.trials));
    else
        numtrials = length(spikes.sweeps.trials);
    end
else
    numtrials=length(unique(spikes.trials));
end
% Set spiketimes
spiketimes = spikes.spiketimes;
% Convert binsize from ms to s
binsize = binsize/1000;
% Get counts
edges = 0:binsize:duration;
n = histc(spiketimes,edges);
n = n/numtrials/binsize;

nsForStdev=zeros(numtrials,size(n,2));
if ~isempty(theseTrials)
    allTrials=theseTrials;
else
    allTrials=unique(spikes.trials);
end
if length(allTrials)~=numtrials
    if ~isempty(theseTrials)
        allTrials=theseTrials;
    elseif length(spikes.sweeps.trials)==numtrials
        allTrials=spikes.sweeps.trials;
    else
        disp('Needed to fill in trials -- be sure you are using contiguous daq files');
        allTrials=min(unique(spikes.trials)):max(unique(spikes.trials));
    end
end      
for i=1:length(allTrials)
    cspikes=filtspikes(spikes,0,'trials',allTrials(i));
    if isempty(cspikes.spiketimes)
        continue
    end
    nsForStdev(i,:)=histc(cspikes.spiketimes,edges);
end
nsForStdev=nsForStdev/binsize;
if all(isnan(n))
    n = 0;
end
% Compute center of bins
centers = edges + diff(edges(1:2))/2;
% Last point of n contains values falling on edge(end) -- usually zero
if bsmooth
    xpoints=centers(1:end-1);
    ypoints=smooth(n(1:end-1),3);
else
    xpoints=centers(1:end-1);
    ypoints=n(1:end-1);
end

varargout{1} = n;
varargout{2} = centers;
varargout{3} = edges;
varargout{4} = xpoints;
varargout{5} = ypoints;
varargout{6} = 1*std(nsForStdev(:,1:end-1),0,1);
end