function [PSTHs,preffed,nonpreffed,x]=testAlpha(spikes,fileInd,includeStimVals,binsize,duration)

ledOffVal=0;
ledOnVal=5;
stimWindow=[1 4];
% avEvWindow=[2.5 3.7];
avEvWindow=[1 4];
avBaseline=[0 0.9];
useAllSpikesPref=0;

s=includeStimVals;
spikes=filtspikes(spikes,0,'fileInd',fileInd);

nTrials=zeros(length(s),2);
theseTrials=cell(length(s),2);
for i=1:length(s)
    subS=filtspikes(spikes,0,'stimcond',s(i),'fileInd',fileInd);
    subS1=filtspikes(subS,0,'led',ledOffVal);
    nTrials(i,1)=length(unique(subS1.trials));
    theseTrials{i,1}=unique(subS1.trials);
    subS1=filtspikes(subS,0,'led',ledOnVal);
    nTrials(i,2)=length(unique(subS1.trials));
    theseTrials{i,2}=unique(subS1.trials);
end    

subS=filtspikes(spikes,0,'led',ledOffVal);
temp{1}=unique(subS.trials);
subS=filtspikes(spikes,0,'led',ledOnVal);
temp{2}=unique(subS.trials);
[x,y1,y2]=makePSTH(spikes,fileInd,ledOffVal,ledOnVal,binsize,sum(nTrials,1),duration,temp);
allSpikes.x=x;
allSpikes.y1=y1;
allSpikes.y2=y2;
figure();
plot(x,y1,'Color','k');
hold on;
plot(x,y2,'Color','r');

a=unique(spikes.assigns);
% s=unique(spikes.stimcond);
ledOff=nan(length(a),length(s),length(y1));
ledOn=nan(length(a),length(s),length(y1));
ledBoth=nan(length(a),length(s),length(y1));
allSpikes=nan(length(s),length(y1));
for i=1:length(s)
    temp{1}=[theseTrials{i,1} theseTrials{i,2}];
    temp{2}=theseTrials{i,2};
    [x,y1,y2]=makePSTH(filtspikes(spikes,0,'stimcond',s(i)),fileInd,[ledOffVal ledOnVal],ledOnVal,binsize,[sum(nTrials(i,:),2) nTrials(i,2)],duration,temp);
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

preffedStim=nan(length(a),1);
preffedVals=nan(length(a),1);
nonpreffedStim=nan(length(a),1);
nonpreffedVals=nan(length(a),1);
xtimes=linspace(0,5,size(ledOff,3));
if useAllSpikesPref==1
    for i=1:length(s)
        evResponse(i)=mean(allSpikes(i,xtimes>=avEvWindow(1) & xtimes<=avEvWindow(2)),2)-mean(allSpikes(i,xtimes>=avBaseline(1) & xtimes<=avBaseline(2)),2);
    end
    [allprefVal,allprefStim]=nanmax(evResponse);
    [allnonprefVal,allnonprefStim]=nanmin(evResponse);
    preffed.stim=ones(size(preffedStim)).*allprefStim;
    preffed.vals=ones(size(preffedVals)).*allprefVal;
    nonpreffed.stim=ones(size(nonpreffedStim)).*allnonprefStim;
    nonpreffed.vals=ones(size(nonpreffedVals)).*allnonprefVal;
else
    for i=1:length(a)
        evResponse=zeros(1,length(s));
        for j=1:length(s)
%             evResponse(j)=mean(ledBoth(i,j,xtimes>=avEvWindow(1) & xtimes<=avEvWindow(2)),3);
%             evResponse(j)=mean(ledOff(i,j,xtimes>=avEvWindow(1) & xtimes<=avEvWindow(2)),3)-mean(ledOff(i,j,xtimes>=avBaseline(1) & xtimes<=avBaseline(2)),3);
              evResponse(j)=mean(ledOff(i,j,xtimes>=avEvWindow(1) & xtimes<=avEvWindow(2)),3);
%             evResponse(j)=mean(ledBoth(i,j,xtimes>=avEvWindow(1) & xtimes<=avEvWindow(2)),3)-mean(ledBoth(i,j,xtimes>=avBaseline(1) & xtimes<=avBaseline(2)),3);
        end
        [preffedVals(i),preffedStim(i)]=nanmax(evResponse);
        [nonpreffedVals(i),nonpreffedStim(i)]=nanmin(evResponse);
    end
    preffed.stim=preffedStim;
    preffed.vals=preffedVals;
    nonpreffed.stim=nonpreffedStim;
    nonpreffed.vals=nonpreffedVals;
end

% for i=1:size(ledBoth,1)
%     figure(); 
%     for j=1:size(ledBoth,2)
%         switch j
%             case 1
%                 c='k';
%             case 2
%                 c='r';
%             case 3
%                 c='b';
%             case 4
%                 c='g';
%             case 5
%                 c='y';
%             case 6
%                 c='c';
%             case 7
%                 c=[0.5 0.5 0.5];
%             case 8
%                 c='m';
%         end      
%         plot(x,reshape(ledBoth(i,j,:),1,length(ledBoth(i,j,:))),'Color',c);
%         hold on;
%     end
% end
end

function [xpoints1,ypoints1,ypoints2]=makePSTH(spikes,useTheseFileInds,noLedValues,ledValue,binsize,nTrials,duration,theseTrials)

useBlackTrials=[];
useRedTrials=[];
% binsize=5; % ms

stimconds=[1:128];

nunits=1;

spikes=filtspikes(spikes,0,'stimcond',stimconds); 

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