function [saveStep2]=bootstrapLedConditions(allSpikes,useTheseAssigns,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED,type1Stimcond,type2Stimcond,ledCond1,ledCond2,PSTHbin,stimResponseWindow,unitsFractionForBootstrap,trialDuration,saveStep2)

nBootstrapTrials=250;
saveSteps=1;
ledOnset=1.3; % in seconds from the start of the trial
normalizePrecedingWindow=1;
precedingWindow=[1.29 1.3]; % s before thalamic silencing
zeroWindow=[1.37 1.38]; % For normalization, when traces should be at zero

if stimResponseWindow(2)>ledOnset
    disp('stimResponseWindow should end before ledOnset');
    return
end

% Take unitsFractionForBootstrap fraction of total units to include in each
% trial of bootstrapping - randomly select unitsFractionForBootstrap for
% each trial
nUnitsPerTrial=floor(unitsFractionForBootstrap*length(useTheseAssigns));
% Get shut-off time course for all units
% Get time course library for all units
if isempty(saveStep2)
    [pref,nonpref]=getType1VsType2TimeCourse(allSpikes,useTheseAssigns,type1Stimcond,type2Stimcond,ledCond1,ledCond2,PSTHbin,trialDuration,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED);
    saveStep2.pref=pref;
    saveStep2.nonpref=nonpref;
    if saveSteps==1
        save('Z:\Updating Data\Analysis Computer\Changing Tau - Anesth\Analysis Steps\saveStep2.mat','saveStep2');
    end
else
    pref=saveStep2.pref;
    nonpref=saveStep2.nonpref;
end
if size(pref.ypoints,1)~=length(useTheseAssigns)
    disp('saveStep2.pref number of rows is not the same as length of current assigns list');
    return
end

all_xpoints=pref.xpoints;
running_pref_ypoints=zeros(nBootstrapTrials,length(all_xpoints));
running_nonpref_ypoints=zeros(nBootstrapTrials,length(all_xpoints));
for i=1:nBootstrapTrials
    % Get random set of units for this bootstrap trial - sample without
    % replacement
    r=randperm(length(useTheseAssigns));
    useTheseIndsForTrial=r(1:nUnitsPerTrial);
%     useTheseAssignsForTrial=allAssigns(useTheseIndsForTrial);
%     useTheseUnitByUnitTrials=unitByUnitTrials(useTheseIndsForTrial);
    % For each trial, get preferred vs. non-preferred time course of
    % cortical shut-off using only useTheseAssignsForTrial
    disp(i);
%     [pr,npr]=getPrefVsNonPrefTimeCourse(filtspikes(allSpikes,0,'assigns',useTheseAssignsForTrial),useTheseAssignsForTrial,prefStimcond(useTheseIndsForTrial),nonprefStimcond(useTheseIndsForTrial),ledOnCond,PSTHbin,trialDuration,useTheseUnitByUnitTrials);
    prefTrace=mean(pref.ypoints(useTheseIndsForTrial,:),1);
    nonprefTrace=mean(nonpref.ypoints(useTheseIndsForTrial,:),1);
    if normalizePrecedingWindow==1
        % Normalize both pref and nonpref traces to 1 during precedingWindow
        % before LED onset and normalize them to 0 during zeroWindow
        useXinds1=all_xpoints>precedingWindow(1) & all_xpoints<precedingWindow(2);
        useXinds0=all_xpoints>zeroWindow(1) & all_xpoints<zeroWindow(2);
        base=mean(prefTrace(useXinds0));
        prefTrace=prefTrace-base;
        top=mean(prefTrace(useXinds1));
        prefTrace=prefTrace/top;
        base=mean(nonprefTrace(useXinds0));
        nonprefTrace=nonprefTrace-base;
        top=mean(nonprefTrace(useXinds1));
        nonprefTrace=nonprefTrace/top;
    end
    running_pref_ypoints(i,:)=prefTrace;
    running_nonpref_ypoints(i,:)=nonprefTrace;
end

% Plot time course figure with error bars from bootstrapped standard
% deviation
figure(); 
hax=axes();
hl=plot(hax,all_xpoints,mean(running_pref_ypoints,1),'Color','k');
addErrBar(all_xpoints,mean(running_pref_ypoints,1),std(running_pref_ypoints,0,1),'y',hax,hl)
% errorbar(all_xpoints,mean(running_pref_ypoints,1),std(running_pref_ypoints,0,1),'Color','k');
hold on; 
hl=plot(hax,all_xpoints,mean(running_nonpref_ypoints,1),'Color','blue');
addErrBar(all_xpoints,mean(running_nonpref_ypoints,1),std(running_nonpref_ypoints,0,1),'y',hax,hl);
%errorbar(all_xpoints,mean(running_nonpref_ypoints,1),std(running_nonpref_ypoints,0,1),'Color','blue');

end


function [pref,nonpref]=getType1VsType2TimeCourse(spikes,useTheseAssigns,type1Stimcond,type2Stimcond,ledCond1,ledCond2,bin,trialDuration,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED)

for i=1:length(useTheseAssigns)
    if i==1
        useSpikes=filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',type1Stimcond);
        currTrialCon=unitByUnitTrials{i};
        currStimCon=unitByUnitStimcond{i};
        currLEDCon=unitByUnitLED{i};
        unitByUnitConsensus=currTrialCon(ismember(currStimCon,type1Stimcond) & ismember(currLEDCon,ledCond1));
        if isempty(unitByUnitConsensus)
            disp('ack');
        end
        [~,~,~,xpoints1,ypoints1]=psth_wStd(filtspikes(useSpikes,0,'led',ledCond1),bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
        prefPSTH_xpoints=xpoints1;
        prefPSTH_ypoints=zeros(length(useTheseAssigns),length(ypoints1));
        prefPSTH_ypoints(i,:)=ypoints1;
        useSpikes=filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',type2Stimcond);
        unitByUnitConsensus=currTrialCon(ismember(currStimCon,type2Stimcond) & ismember(currLEDCon,ledCond2));
        if isempty(unitByUnitConsensus)
            disp('ack');
        end
        [~,~,~,xpoints1,ypoints1]=psth_wStd(filtspikes(useSpikes,0,'led',ledCond2),bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
        nonprefPSTH_xpoints=xpoints1;
        nonprefPSTH_ypoints=zeros(length(useTheseAssigns),length(ypoints1));
        nonprefPSTH_ypoints(i,:)=ypoints1;
    else
        useSpikes=filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',type1Stimcond);
        currTrialCon=unitByUnitTrials{i};
        currStimCon=unitByUnitStimcond{i};
        currLEDCon=unitByUnitLED{i};
        unitByUnitConsensus=currTrialCon(ismember(currStimCon,type1Stimcond) & ismember(currLEDCon,ledCond1));
        if isempty(unitByUnitConsensus)
            disp('ack');
        end
        [~,~,~,xpoints1,ypoints1]=psth_wStd(filtspikes(useSpikes,0,'led',ledCond1),bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
        prefPSTH_ypoints(i,:)=ypoints1;
        useSpikes=filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',type2Stimcond);
        unitByUnitConsensus=currTrialCon(ismember(currStimCon,type2Stimcond) & ismember(currLEDCon,ledCond2));
        if isempty(unitByUnitConsensus)
            disp('ack');
        end
        [~,~,~,xpoints1,ypoints1]=psth_wStd(filtspikes(useSpikes,0,'led',ledCond2),bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
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
varargout{1} = m;
varargout{2} = s;
varargout{3} = mean(nsForStdev,2);
end