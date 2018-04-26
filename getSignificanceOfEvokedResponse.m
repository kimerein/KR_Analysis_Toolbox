function [isSigStim,prefR,prefStimcond,nonprefR,nonprefStimcond,prefSigStim]=getSignificanceOfEvokedResponse(spikes,useTheseAssigns,stim_useTheseStimcond,blank_useTheseStimcond,useLEDcond,onsetResponseWindow)

prefR=zeros(length(useTheseAssigns),1);
prefStimcond=zeros(length(useTheseAssigns),1);
nonprefR=zeros(length(useTheseAssigns),1);
nonprefStimcond=zeros(length(useTheseAssigns),1);
sigStim=zeros(length(useTheseAssigns),1);
[currTrialCon,ind]=unique(spikes.trials);
currStimCon=spikes.stimcond(ind);
currLedCon=spikes.led(ind);
prefSigStim=zeros(length(useTheseAssigns),1);
for i=1:length(useTheseAssigns)
    unitByUnitConsensus=currTrialCon(ismember(currStimCon,blank_useTheseStimcond) & ismember(currLedCon,useLEDcond));
    if isempty(unitByUnitConsensus)
        disp('ack');
    end
    [blankR,temp,d_blank]=calcMeanAndStdForUnit(filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',blank_useTheseStimcond),onsetResponseWindow,unitByUnitConsensus);
    stimcondR=zeros(length(stim_useTheseStimcond),1);
    sigStim(i)=0;
    for j=1:length(stim_useTheseStimcond)
        unitByUnitConsensus=currTrialCon(ismember(currStimCon,stim_useTheseStimcond(j)) & ismember(currLedCon,useLEDcond));
        if isempty(unitByUnitConsensus)
            disp('ack');
        end
        [stimcondR(j),temp,d{j}]=calcMeanAndStdForUnit(filtspikes(spikes,0,'assigns',useTheseAssigns(i),'stimcond',stim_useTheseStimcond(j)),onsetResponseWindow,unitByUnitConsensus);
        sigVal=mattest(d{j}',d_blank');
        if sigVal<0.05
            sigStim(i)=1;
        end
    end
    [prefR(i),ind]=max(stimcondR);
    prefStimcond(i)=stim_useTheseStimcond(ind);
    [nonprefR(i),ind]=min(stimcondR);
    nonprefStimcond(i)=stim_useTheseStimcond(ind);
    prefSigStim(i)=mattest(d{prefStimcond(i)}',d{nonprefStimcond(i)}');
end
isSigStim=sigStim;
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