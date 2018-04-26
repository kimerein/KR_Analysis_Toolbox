function response=getSignificantPreferredResponse(allSpikes)

useTheseStimcond=1:8;
stimResponseWindow=[1 4];
nBest=1;
nWorst=1;
stimResponseLED=0;

% Get significance of stimulus responses for all units
allAssigns=unique(allSpikes.assigns);

for i=1:length(allAssigns)
    unitByUnitTrials{i}=allSpikes.sweeps.trials;
    unitByUnitStimcond{i}=allSpikes.sweeps.stimcond;
    unitByUnitLED{i}=allSpikes.sweeps.led;
end

[prefR,prefStimcond,nonprefR,nonprefStimcond,sigStim]=getSignificanceOfPreferredResponse(allSpikes,allAssigns,useTheseStimcond,stimResponseWindow,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED,nBest,nWorst,stimResponseLED);
response.prefR=prefR;
response.prefStimcond=prefStimcond;
response.nonprefR=nonprefR;
response.nonprefStimcond=nonprefStimcond;
response.sigStim=sigStim;

end

function FRsAcrossStims=getFRsAcrossStims(spikes,useTheseAssigns,stimcondsToUse,ledCond1,window,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED)

FRsAcrossStims=zeros(length(useTheseAssigns),length(stimcondsToUse));
for i=1:length(stimcondsToUse)
    FRsAcrossStims(:,i)=getType1WindowFR(spikes,useTheseAssigns,stimcondsToUse(i),ledCond1,window,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED);
end
end

function meanFRsForUnits=getType1WindowFR(spikes,useTheseAssigns,type1Stimcond,ledCond1,window,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED)

meanFRsForUnits=zeros(length(useTheseAssigns),1);
for i=1:length(useTheseAssigns)
    useSpikes=filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',type1Stimcond);
    currTrialCon=unitByUnitTrials{i};
    currStimCon=unitByUnitStimcond{i};
    currLEDCon=unitByUnitLED{i};
    unitByUnitConsensus=currTrialCon(ismember(currStimCon,type1Stimcond) & ismember(currLEDCon,ledCond1));
    if isempty(unitByUnitConsensus)
        disp('unitByUnitConsensus is empty');
    end
    m=calcMeanAndStdForUnit(filtspikes(useSpikes,0,'led',ledCond1),window,unitByUnitConsensus);
    meanFRsForUnits(i)=m;
end
end

function [prefR,prefStimcond,nonprefR,nonprefStimcond,sigStim]=getSignificanceOfPreferredResponse(spikes,useTheseAssigns,useTheseStimcond,onsetResponseWindow,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED,nBest,nWorst,useLED)

prefR=zeros(length(useTheseAssigns),1);
prefStimcond=cell(length(useTheseAssigns),1);
nonprefR=zeros(length(useTheseAssigns),1);
nonprefStimcond=cell(length(useTheseAssigns),1);
sigStim=zeros(length(useTheseAssigns),1);
for i=1:length(useTheseAssigns)
    stimcondR=zeros(length(useTheseStimcond),1);
    for j=1:length(useTheseStimcond)
        currTrialCon=unitByUnitTrials{i};
        currStimCon=unitByUnitStimcond{i};
        currLEDCon=unitByUnitLED{i};
        unitByUnitConsensus=currTrialCon(ismember(currStimCon,useTheseStimcond(j)) & ismember(currLEDCon,useLED));
        if isempty(unitByUnitConsensus)
            disp('unitByUnitConsensus is empty');
        end
        [stimcondR(j),temp,d{j}]=calcMeanAndStdForUnit(filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',useTheseStimcond(j)),onsetResponseWindow,unitByUnitConsensus);
    end
    [sortVals,sortInd]=sort(stimcondR,'descend');
    prefR(i)=mean(sortVals(1:nBest));
    ind=sortInd(1:nBest);
    prefStimcond{i}=useTheseStimcond(ind);
    maxprefStimcond=useTheseStimcond(ind(1));
    [sortVals,sortInd]=sort(stimcondR,'ascend');
    nonprefR(i)=mean(sortVals(1:nWorst));
    ind=sortInd(1:nWorst);
    nonprefStimcond{i}=useTheseStimcond(ind);
    minnonprefStimcond=useTheseStimcond(ind(1));
    sigStim(i)=mattest(d{maxprefStimcond}',d{minnonprefStimcond}');
end
disp('Done with calculating preferred and non-preferred responses');
end

function [pref,nonpref]=getPrefVsNonPrefTimeCourse(spikes,useTheseAssigns,prefStimcond,nonprefStimcond,ledOnCond,bin,trialDuration,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED)

for i=1:length(useTheseAssigns)
    if i==1
        useSpikes=filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',prefStimcond{i});
        currTrialCon=unitByUnitTrials{i};
        currStimCon=unitByUnitStimcond{i};
        currLEDCon=unitByUnitLED{i};
        unitByUnitConsensus=currTrialCon(ismember(currStimCon,prefStimcond{i}) & ismember(currLEDCon,ledOnCond));
        if isempty(unitByUnitConsensus)
            disp('ack');
        end
        [~,~,~,xpoints1,ypoints1]=psth_wStd(filtspikes(useSpikes,0,'led',ledOnCond),bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
        prefPSTH_xpoints=xpoints1;
        prefPSTH_ypoints=zeros(length(useTheseAssigns),length(ypoints1));
        prefPSTH_ypoints(i,:)=ypoints1;
        useSpikes=filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',nonprefStimcond{i});
        unitByUnitConsensus=currTrialCon(ismember(currStimCon,nonprefStimcond{i}) & ismember(currLEDCon,ledOnCond));
        if isempty(unitByUnitConsensus)
            disp('ack');
        end
        [~,~,~,xpoints1,ypoints1]=psth_wStd(filtspikes(useSpikes,0,'led',ledOnCond),bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
        nonprefPSTH_xpoints=xpoints1;
        nonprefPSTH_ypoints=zeros(length(useTheseAssigns),length(ypoints1));
        nonprefPSTH_ypoints(i,:)=ypoints1;
    else
        useSpikes=filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',prefStimcond{i});
        currTrialCon=unitByUnitTrials{i};
        currStimCon=unitByUnitStimcond{i};
        currLEDCon=unitByUnitLED{i};
        unitByUnitConsensus=currTrialCon(ismember(currStimCon,prefStimcond{i}) & ismember(currLEDCon,ledOnCond));
        if isempty(unitByUnitConsensus)
            disp('ack');
        end
        [~,~,~,xpoints1,ypoints1]=psth_wStd(filtspikes(useSpikes,0,'led',ledOnCond),bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
        prefPSTH_ypoints(i,:)=ypoints1;
        useSpikes=filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',nonprefStimcond{i});
        unitByUnitConsensus=currTrialCon(ismember(currStimCon,nonprefStimcond{i}) & ismember(currLEDCon,ledOnCond));
        
        [~,~,~,xpoints1,ypoints1]=psth_wStd(filtspikes(useSpikes,0,'led',ledOnCond),bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
        nonprefPSTH_ypoints(i,:)=ypoints1;
    end
end
pref.xpoints=prefPSTH_xpoints;
pref.ypoints=prefPSTH_ypoints;
nonpref.xpoints=nonprefPSTH_xpoints;
nonpref.ypoints=nonprefPSTH_ypoints;
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

function [varargout]=calcMeanAndStdForUnit(spikes,window,theseTrials)

binsize=1; % in ms
if isempty(theseTrials)
    disp('theseTrials should not be empty in calcMeanAndStdForUnit');
end
numtrials=length(theseTrials);
% Set spiketimes
spiketimes = spikes.spiketimes;
% Convert binsize from ms to s
binsize = binsize/1000;
% Get counts
edges=window(1):binsize:window(2);
n = histc(spiketimes,edges);
n = n/numtrials/binsize;
m = mean(n);

nsForStdev=zeros(numtrials,size(n,2));
allTrials=theseTrials;
for i=1:length(allTrials)
    cspikes=filtpartialspikes(spikes,0,'trials',allTrials(i));
    if isempty(cspikes.spiketimes)
        continue
    end
    nsForStdev(i,:)=histc(cspikes.spiketimes,edges);
end
nsForStdev=nsForStdev/binsize;
s=std(mean(nsForStdev,2),0,1);
if all(isnan(n))
    n = 0;
end
varargout{1} = mean(mean(nsForStdev,2));
varargout{2} = s;
varargout{3} = mean(nsForStdev,2);
end