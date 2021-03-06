function [saveStep1,saveStep2,FRsAcrossStims1,FRsAcrossStims2]=compareUnitAverageRankedTuning(allSpikes,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED,useTheseStimcond,useOnlySigPrefUnits,ledCond1,ledCond2,stimResponseWindow,stimResponseLED,saveStep1,saveStep2)

stimcondsToUse=useTheseStimcond;
considerWindow1=[1.325 1.4];
considerWindow2=[1.325 1.4];
% considerWindow1=[1 4];
% considerWindow2=[1 4];
% considerWindow=[1.125 1.2];
% considerWindow=[1 4];
nBest=1;
nWorst=1;
significanceThresh=0.0001;
saveSteps=1;
rankFRs=1;
considerPrefDirAndOrientSeparately=1;
alignToPreferred=1;
alignToTuning=0;
rankDirOrient=0;

% Get significance of stimulus responses for all units
allAssigns=unique(allSpikes.assigns);

if isempty(saveStep2) && isempty(saveStep1)
    [prefR,prefStimcond,nonprefR,nonprefStimcond,sigStim]=getSignificanceOfPreferredResponse(allSpikes,allAssigns,useTheseStimcond,stimResponseWindow,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED,nBest,nWorst,stimResponseLED);
    saveStep1.prefR=prefR;
    saveStep1.prefStimcond=prefStimcond;
    saveStep1.nonprefR=nonprefR;
    saveStep1.nonprefStimcond=nonprefStimcond;
    saveStep1.sigStim=sigStim;
    if saveSteps==1
        save('Z:\Updating Data\Analysis Computer\Changing Tau - Anesth\Information\saveStep1.mat','saveStep1');
    end
elseif ~isempty(saveStep1)
    prefR=saveStep1.prefR;
    prefStimcond=saveStep1.prefStimcond;
    nonprefR=saveStep1.nonprefR;
    nonprefStimcond=saveStep1.nonprefStimcond;
    sigStim=saveStep1.sigStim;
end

if isempty(saveStep1) && ~isempty(saveStep2)
    disp('Need the data from both saveStep1 and saveStep2 to skip ahead in code');
    return
end

% If useOnlySigForStimUnits==1, use only the units with a significant
% preferred stimulus, that is, response to preferred stimulus is
% SIGNIFICANTLY greater than response to non-preferred stimulus
% Preferred and non-preferred orientations are not forced to be 90 deg.
% apart
if useOnlySigPrefUnits==1
    if length(sigStim)~=length(allAssigns)
        disp('sigStim and allAssigns lengths do not match');
        return
    else
        allSpikes=filtspikes(allSpikes,0,'assigns',allAssigns(sigStim<significanceThresh));
        unitByUnitTrials=unitByUnitTrials(sigStim<significanceThresh);
        unitByUnitStimcond=unitByUnitStimcond(sigStim<significanceThresh);
        unitByUnitLED=unitByUnitLED(sigStim<significanceThresh);
        allAssigns=unique(allSpikes.assigns);
        prefR=prefR(sigStim<significanceThresh);
        prefStimcond=prefStimcond(sigStim<significanceThresh);
        nonprefR=nonprefR(sigStim<significanceThresh);
        nonprefStimcond=nonprefStimcond(sigStim<significanceThresh);
    end
end
if length(allAssigns)~=length(unitByUnitTrials)
    disp('units count is off');
    return
end
disp('n units');
disp(length(allAssigns));

% Get FR library for all units
if isempty(saveStep2)
    FRsAcrossStims1=getFRsAcrossStims(allSpikes,allAssigns,stimcondsToUse,ledCond1,considerWindow1,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED);
    FRsAcrossStims2=getFRsAcrossStims(allSpikes,allAssigns,stimcondsToUse,ledCond2,considerWindow2,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED);
    saveStep2.ledCond1=FRsAcrossStims1;
    saveStep2.ledCond2=FRsAcrossStims2;
    if saveSteps==1
        save('Z:\Updating Data\Analysis Computer\Changing Tau - Anesth\Information\saveStep2.mat','saveStep2');
    end
else
    FRsAcrossStims1=saveStep2.ledCond1;
    FRsAcrossStims2=saveStep2.ledCond2;
end
if size(saveStep2.ledCond1,1)~=length(allAssigns)
    disp('saveStep2.ledCond1 number of rows is not the same as length of current assigns list');
    return
end

if considerPrefDirAndOrientSeparately==1
    % Align tuning curves to max stimulus condition for first data set
    if alignToTuning==1
        % Align tuning curves to preferred or non-preferred orientation
        for i=1:size(FRsAcrossStims1,1)
            if alignToPreferred==1
                currAlignStim=prefStimcond(i);
            else
                currAlignStim=nonprefStimcond(i);
            end
            currAlignInd=find(stimcondsToUse==currAlignStim{1});
            temp=FRsAcrossStims1(i,:);
            FRsAcrossStims1(i,1:length(temp(currAlignInd:end)))=temp(currAlignInd:end);
            if currAlignInd>1
                FRsAcrossStims1(i,length(temp(currAlignInd:end))+1:end)=temp(1:currAlignInd-1);
            end
            temp=FRsAcrossStims2(i,:);
            FRsAcrossStims2(i,1:length(temp(currAlignInd:end)))=temp(currAlignInd:end);
            if currAlignInd>1
                FRsAcrossStims2(i,length(temp(currAlignInd:end))+1:end)=temp(1:currAlignInd-1);
            end
        end
    else
        for i=1:size(FRsAcrossStims1,1)
            if alignToPreferred==1
                [ma,in]=max(FRsAcrossStims1(i,:),[],2);
            else
                [ma,in]=min(FRsAcrossStims1(i,:),[],2);
            end
            currAlignInd=in;
            temp=FRsAcrossStims1(i,:);
            FRsAcrossStims1(i,1:length(temp(currAlignInd:end)))=temp(currAlignInd:end);
            if currAlignInd>1
                FRsAcrossStims1(i,length(temp(currAlignInd:end))+1:end)=temp(1:currAlignInd-1);
            end
            temp=FRsAcrossStims2(i,:);
            FRsAcrossStims2(i,1:length(temp(currAlignInd:end)))=temp(currAlignInd:end);
            if currAlignInd>1
                FRsAcrossStims2(i,length(temp(currAlignInd:end))+1:end)=temp(1:currAlignInd-1);
            end
        end
    end
    % So preferred (if alignToPreferred==1) or non-preferred stimulus
    % condition should now be first
    % Consider preferred direction and preferred orientation for each
    % unit separately
    % Preferred direction is the preferred stimulus and the two neighboring
    % stimuli
    % Preferred orientation collapses similarly oriented stimuli (across
    % different directions)
    dir1=zeros(size(FRsAcrossStims1,1),2);
    orient1=zeros(size(FRsAcrossStims1,1),4);
    dir2=zeros(size(FRsAcrossStims2,1),2);
    orient2=zeros(size(FRsAcrossStims2,1),4);
    for i=1:size(FRsAcrossStims1,1)
        dir1(i,1)=mean(FRsAcrossStims1(i,[1 2 8]));
        dir1(i,2)=mean(FRsAcrossStims1(i,[4 5 6]));
        dir2(i,1)=mean(FRsAcrossStims2(i,[1 2 8]));
        dir2(i,2)=mean(FRsAcrossStims2(i,[4 5 6]));
        orient1(i,1)=mean(FRsAcrossStims1(i,[1 5]));
        orient1(i,2)=mean(FRsAcrossStims1(i,[2 6]));
        orient1(i,3)=mean(FRsAcrossStims1(i,[3 7]));
        orient1(i,4)=mean(FRsAcrossStims1(i,[4 8]));
        orient2(i,1)=mean(FRsAcrossStims2(i,[1 5]));
        orient2(i,2)=mean(FRsAcrossStims2(i,[2 6]));
        orient2(i,3)=mean(FRsAcrossStims2(i,[3 7]));
        orient2(i,4)=mean(FRsAcrossStims2(i,[4 8]));
    end
end

if rankDirOrient==1
%     for i=1:size(dir1,1)
%         [~,in]=sort(dir1(i,:));
%         dir1(i,in)=1:length(in);
%     end
%     for i=1:size(dir2,1)
%         [~,in]=sort(dir2(i,:));
%         dir2(i,in)=1:length(in);
%     end
    for i=1:size(orient1,1)
        [~,in]=sort(orient1(i,:));
        orient1(i,in)=1:length(in);
    end
    for i=1:size(orient2,1)
        [~,in]=sort(orient2(i,:));
        orient2(i,in)=1:length(in);
    end
else
    mins1=min(orient1,[],2);
    useMins1=zeros(size(orient1,1),size(orient1,2));
    for i=1:size(orient1,2)
        useMins1(:,i)=mins1;
    end
    orient1=orient1-useMins1;
    maxs1=max(orient1,[],2);
    for i=1:size(orient1,1)
        if maxs1(i)==0
            orient1(i,:)=ones(size(orient1(i,:)));
        else
            currScale=1/maxs1(i);
            orient1(i,:)=orient1(i,:)*currScale;
        end
    end
    mins1=min(orient2,[],2);
    useMins1=zeros(size(orient2,1),size(orient2,2));
    for i=1:size(orient2,2)
        useMins1(:,i)=mins1;
    end
    orient2=orient2-useMins1;
    maxs1=max(orient2,[],2);
    for i=1:size(orient2,1)
        if maxs1(i)==0
            orient2(i,:)=ones(size(orient2(i,:)));
        else
            currScale=1/maxs1(i);
            orient2(i,:)=orient2(i,:)*currScale;
        end
    end
end

if rankFRs==1
    for i=1:size(FRsAcrossStims1,1)
        [~,in]=sort(FRsAcrossStims1(i,:));
        FRsAcrossStims1(i,in)=1:length(in);
    end
    for i=1:size(FRsAcrossStims2,1)
        [~,in]=sort(FRsAcrossStims2(i,:));
        FRsAcrossStims2(i,in)=1:length(in);
    end
end    
    
if considerPrefDirAndOrientSeparately==1
    figure(); 
    a=-0.4;
    b=0.4;
    for i=1:size(dir1,1)
        r1 = a + (b-a).*rand(length(dir1(i,:)),1);
        r2 = a + (b-a).*rand(length(dir2(i,:)),1);
        plot(dir1(i,:)+r1',dir2(i,:)+r2','Color','b');
        hold on;
    end
    figure(); 
    for i=1:size(orient1,1)
%         r1 = a + (b-a).*rand(length(orient1(i,:)),1);
%         r2 = a + (b-a).*rand(length(orient2(i,:)),1);
%         scatter(orient1(i,:)+r1',orient2(i,:)+r2',[],'k');
        scatter(orient1(i,:),orient2(i,:),[],'k');
        hold all; 
%         [sFRs,in]=sort(orient1(i,:));
%         plot(sFRs+r1(in)',orient2(i,in)+r2(in)');
        [sFRs,in]=sort(orient1(i,:));
        plot(sFRs,orient2(i,in));
    end
end

% Plot ranked stimcond relationships
figure();
concatFRs1=[];
concatFRs2=[];
for i=1:size(FRsAcrossStims1,1)
    concatFRs1=[concatFRs1 FRsAcrossStims1(i,:)];
    concatFRs2=[concatFRs2 FRsAcrossStims2(i,:)];
end
% concatFRs1(concatFRs1==1)=9;
% concatFRs1(concatFRs1==2)=10;
% concatFRs1(concatFRs1==3)=11;
% concatFRs1(concatFRs1==4)=12;
a=-0.4;
b=0.4;
r1 = a + (b-a).*rand(length(concatFRs1),1);
r2 = a + (b-a).*rand(length(concatFRs2),1);
scatter(concatFRs1+r1',concatFRs2+r2',[],'k');
[r,p,rlo,rup]=corrcoef(concatFRs1,concatFRs2);
disp('corrcoef');
disp(r);
disp('pval');
disp(p);

figure(); 
for i=1:size(FRsAcrossStims1,1)
    a=-0.1;
    b=0.1;
    r1 = a + (b-a).*rand(length(FRsAcrossStims1(i,:)),1);
    r2 = a + (b-a).*rand(length(FRsAcrossStims2(i,:)),1);
    scatter(FRsAcrossStims1(i,:)+r1',FRsAcrossStims2(i,:)+r2',[],'k');
    hold all;
    [sFRs,in]=sort(FRsAcrossStims1(i,:));
    plot(sFRs+r1(in)',FRsAcrossStims2(i,in)+r2(in)');
end

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