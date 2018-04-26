function [fx,spikes_ledNum,spikes_ledDen,trials_nsNum,trials_nsDen,assignsToUse]=getUnitByUnitEffects(spikes,trialCount_spikes,assignsToUse,ledCondNum,ledCondDen,stimCondNum,stimCondDen,window)

% Calculates, for each unit in assignsToUse, the value
% mean(spikes for that unit with ledCondNum during window)/
% mean(spikes for that unit with ledCondDen during window)
% 
% Returns this value for each unit in a vector of length(assignsToUse)
% where each index into fx corresponds to the indexed unit in assignsToUse

trials_nsNum=[];
trials_nsDen=[];
for i=1:length(assignsToUse)
%     disp(assignsToUse(i));
%     sp=filtspikes(spikes,0,'stimcond',stimCondNum,'led',ledCondNum);
    sp1=filtspikes(trialCount_spikes,0,'stimcond',stimCondNum);
    sp1=prep_filtspikes(sp1,'led',ledCondNum,[]);
    trialsAll=unique(sp1.trials);
    sp=filtspikes(spikes,0,'stimcond',stimCondNum);
    sp=prep_filtspikes(sp,'led',ledCondNum,[]);
    sp=filtspikes(sp,0,'assigns',assignsToUse(i));
%     sp=filtspikes(spikes,0,'assigns',assignsToUse(i),'led',ledCondNum);
%     sp=filtspikes(sp,0,'stimcond',stimCondNum);
    [m_ledCondNum,~,nsNum]=calcMeanAndStdForUnit(sp,window,trialsAll);
%     sp=filtspikes(spikes,0,'stimcond',stimCondDen,'led',ledCondDen);
    sp1=filtspikes(trialCount_spikes,0,'stimcond',stimCondDen);
    sp1=prep_filtspikes(sp1,'led',ledCondDen,[]);
    trialsAll=unique(sp1.trials);
    sp=filtspikes(spikes,0,'stimcond',stimCondDen);
    sp=prep_filtspikes(sp,'led',ledCondDen,[]);
    sp=filtspikes(sp,0,'assigns',assignsToUse(i));
%     sp=filtspikes(spikes,0,'assigns',assignsToUse(i),'led',ledCondDen);
%     sp=filtspikes(sp,0,'stimcond',stimCondDen);
    [m_ledCondDen,~,nsDen]=calcMeanAndStdForUnit(sp,window,trialsAll);
    if isnan(m_ledCondNum)
        disp('isnan');
    elseif isnan(m_ledCondDen)
        disp('isnan2');
    end
%     trials_nsNum=[trials_nsNum; nsNum];
%     trials_nsDen=[trials_nsDen; nsDen];
%     trials_nsNum=[trials_nsNum; nsNum(2:2:end)];
%     trials_nsDen=[trials_nsDen; nsDen(2:2:end)];
    trials_nsNum=[trials_nsNum nsNum];
    trials_nsDen=[trials_nsDen nsDen];
    if m_ledCondDen==0
        fx(i)=nan;
        spikes_ledNum(i)=m_ledCondNum;
        spikes_ledDen(i)=m_ledCondDen;
    else
        fx(i)=m_ledCondNum/m_ledCondDen;
        spikes_ledNum(i)=m_ledCondNum;
        spikes_ledDen(i)=m_ledCondDen;
    end
end
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