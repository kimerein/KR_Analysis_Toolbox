function getNewUnitsFx_forFreqs(bestSpikes,spikes,params,saveDirName)

window1=[0 0.999];
window2=[1.3 3];
leds=params.useLEDcond{1};
baseUnitsFx.fr=cell(1,length(leds));
unitsFx.fr=cell(1,length(leds));
for i=1:length(leds)-1
    disp(i);
    currLED1=leds(i);
    currLED2=leds(i+1);
    [led1Spikes,led2Spikes,led1Spikes2,led2Spikes2,ass]=getUnitByUnitEffects(bestSpikes,spikes,unique(bestSpikes.assigns),currLED1,currLED2,params.useStimcond,params.useStimcond,window1,window2);
    baseUnitsFx.fr{i}=led1Spikes;
    baseUnitsFx.fr{i+1}=led2Spikes;
    baseUnitsFx.assigns=ass;
    unitsFx.fr{i}=led1Spikes2;
    unitsFx.fr{i+1}=led2Spikes2;
    unitsFx.assigns=ass;
end
save([saveDirName '\baseUnitsFx.mat'],'baseUnitsFx');
save([saveDirName '\unitsFx.mat'],'unitsFx');
end

function [spikes_ledNum,spikes_ledDen,spikes_ledNum2,spikes_ledDen2,assignsToUse]=getUnitByUnitEffects(spikes,trialCount_spikes,assignsToUse,ledCondNum,ledCondDen,stimCondNum,stimCondDen,window1,window2)

% Calculates, for each unit in assignsToUse, the value
% mean(spikes for that unit with ledCondNum during window)/
% mean(spikes for that unit with ledCondDen during window)
% 
% Returns this value for each unit in a vector of length(assignsToUse)
% where each index into fx corresponds to the indexed unit in assignsToUse

% trials_nsNum=[];
% trials_nsDen=[];
spikes_ledNum=zeros(1,length(assignsToUse));
spikes_ledDen=zeros(1,length(assignsToUse));
spikes_ledNum2=zeros(1,length(assignsToUse));
spikes_ledDen2=zeros(1,length(assignsToUse));
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
    [m_ledCondNum,m_ledCondNum2]=calcMeanAndStdForUnit(sp,window1,window2,trialsAll);
%     sp=filtspikes(spikes,0,'stimcond',stimCondDen,'led',ledCondDen);
    sp1=filtspikes(trialCount_spikes,0,'stimcond',stimCondDen);
    sp1=prep_filtspikes(sp1,'led',ledCondDen,[]);
    trialsAll=unique(sp1.trials);
    sp=filtspikes(spikes,0,'stimcond',stimCondDen);
    sp=prep_filtspikes(sp,'led',ledCondDen,[]);
    sp=filtspikes(sp,0,'assigns',assignsToUse(i));
%     sp=filtspikes(spikes,0,'assigns',assignsToUse(i),'led',ledCondDen);
%     sp=filtspikes(sp,0,'stimcond',stimCondDen);
    [m_ledCondDen,m_ledCondDen2]=calcMeanAndStdForUnit(sp,window1,window2,trialsAll);
    if isnan(m_ledCondNum)
        disp('isnan');
    elseif isnan(m_ledCondDen)
        disp('isnan2');
    end
% %     trials_nsNum=[trials_nsNum; nsNum];
% %     trials_nsDen=[trials_nsDen; nsDen];
% %     trials_nsNum=[trials_nsNum; nsNum(2:2:end)];
% %     trials_nsDen=[trials_nsDen; nsDen(2:2:end)];
%     trials_nsNum=[trials_nsNum nsNum];
%     trials_nsDen=[trials_nsDen nsDen];
    if m_ledCondDen==0
        spikes_ledNum(i)=m_ledCondNum;
        spikes_ledDen(i)=m_ledCondDen;
        spikes_ledNum2(i)=m_ledCondNum2;
        spikes_ledDen2(i)=m_ledCondDen2;
    else
        spikes_ledNum(i)=m_ledCondNum;
        spikes_ledDen(i)=m_ledCondDen;
        spikes_ledNum2(i)=m_ledCondNum2;
        spikes_ledDen2(i)=m_ledCondDen2;
    end
end
end


function [varargout]=calcMeanAndStdForUnit(spikes,window1,window2,theseTrials)

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
edges=window1(1):binsize:window1(2);
n = histc(spiketimes,edges);
n = n/numtrials/binsize;
m = mean(n);

% nsForStdev=zeros(numtrials,size(n,2));
% allTrials=theseTrials;
% for i=1:length(allTrials)
%     cspikes=filtpartialspikes(spikes,0,'trials',allTrials(i));
%     if isempty(cspikes.spiketimes)
%         continue
%     end
%     nsForStdev(i,:)=histc(cspikes.spiketimes,edges);
% end
% nsForStdev=nsForStdev/binsize;
% s=std(mean(nsForStdev,2),0,1);
% if all(isnan(n))
%     n = 0;
% end
% varargout{1} = mean(mean(nsForStdev,2));
varargout{1} = m;

edges=window2(1):binsize:window2(2);
n = histc(spiketimes,edges);
n = n/numtrials/binsize;
m = mean(n);

varargout{2} = m;
% varargout{3} = mean(nsForStdev,2);
end