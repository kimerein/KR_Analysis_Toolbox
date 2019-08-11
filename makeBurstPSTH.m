function [fractionpsth,allspikespsth,burstpsth,unitISIplots,spikes]=makeBurstPSTH(spikes,psth,noThetaTrials)

% trialDuration=14.5;
trialDuration=4;

spikes.trials=spikes.trials-min(spikes.trials)+1;
spikes.sweeps.trials=spikes.sweeps.trials-min(spikes.sweeps.trials)+1;
tri=psth.unitTrials{1};
spikes=filtspikes(spikes,0,'trials',tri);

u=unique(spikes.assigns);
if length(psth.psths)>length(u)
    psth.psths=psth.psths(1:length(u));
    psth.unitTrials=psth.unitTrials(1:length(u));
    psth.unitStimcond=psth.unitStimcond(1:length(u));
    psth.unitLED=psth.unitLED(1:length(u));
end

% timeBins=0:0.01:trialDuration;
timeBins=psth.t-((psth.t(2)-psth.t(1))/2);
useStimcond={[1:100000]};
% useLED=[0 0.05 0.1515 0.2525 0.7575 5.00 5.05];
useLED=unique(psth.unitLED{1});
useAssigns=unique(spikes.assigns);
bin=10; % in ms
trialDuration=4;
spikes.isFirstSpikeInBurst=zeros(size(spikes.trials));
spikes.isFollowingSpikeInBurst=zeros(size(spikes.trials));
spikes.isInBurst=zeros(size(spikes.trials));
unitISIplots=cell(length(psth.psths),2);
for i=1:length(psth.psths)
    curru=u(i);
    subspikes=filtspikes(spikes,0,'trials',tri,'assigns',curru);
    [firstSpikes,followSpikes,preISI,postISI,isThirdSpike]=sub_findBursts(subspikes.spiketimes,trialDuration);
    unitISIplots{i,1}=preISI;
    unitISIplots{i,2}=postISI;
    subspikes.isFirstSpikeInBurst=firstSpikes;
    subspikes.isFollowingSpikeInBurst=followSpikes;
    subspikes.isInBurst=firstSpikes | followSpikes;
    spikes.isFirstSpikeInBurst(ismember(spikes.trials,tri) & ismember(spikes.assigns,curru))=firstSpikes;
    spikes.isFollowingSpikeInBurst(ismember(spikes.trials,tri) & ismember(spikes.assigns,curru))=followSpikes;
    spikes.isThirdSpike(ismember(spikes.trials,tri) & ismember(spikes.assigns,curru))=isThirdSpike;
    spikes.isInBurst(ismember(spikes.trials,tri) & ismember(spikes.assigns,curru))=firstSpikes | followSpikes;
end

plotExampleRaster=1;
if plotExampleRaster==1
    figure();
    hax=axes();
    spikes.isSynchSpike=spikes.isInBurst;
    u=unique(spikes.assigns);
    raster_sub(filtspikes(spikes,0,'assigns',u(1)),hax,0,0,14.5,[]);
%     trialz=unique(spikes.trials);
%     if length(trialz)~=length(noThetaTrials)
%         disp('problem');
%     end
%     ss=filtspikes(spikes,0,'assigns',u(19));
%     ss.led=floor(ss.led);
%     ss.sweeps.led=floor(ss.sweeps.led);
%     raster_sub(filtspikes(ss,0,'trials',trialz(noThetaTrials),'led',0),hax,0,0,14.5,[]);
%     figure();
%     hax=axes();
%     raster_sub(filtspikes(ss,0,'trials',trialz(noThetaTrials),'led',5),hax,0,0,14.5,[]);
%     figure();
%     hax=axes();
%     raster_sub(filtspikes(ss,0,'trials',trialz(~noThetaTrials),'led',0),hax,0,0,14.5,[]);
%     figure();
%     hax=axes();
%     raster_sub(filtspikes(ss,0,'trials',trialz(~noThetaTrials),'led',5),hax,0,0,14.5,[]);
end
% figure();
burstSpikes=filtspikes(spikes,0,'isInBurst',1);
[x,psths_t,unitTrials,unitStimcond,unitLED]=getTrialByTrialUnitPSTH_sub(burstSpikes,useAssigns,useLED,bin,trialDuration,useStimcond);
x=[x max(x)+(x(2)-x(1))];
bpsth.t=x;
bpsth.psths=psths_t;
bpsth.unitTrials=unitTrials;
bpsth.unitStimcond=unitStimcond;
bpsth.unitLED=unitLED;
burstpsth=bpsth;
[x,psths_t,unitTrials,unitStimcond,unitLED]=getTrialByTrialUnitPSTH_sub(spikes,useAssigns,useLED,bin,trialDuration,useStimcond);
x=[x max(x)+(x(2)-x(1))];
bpsth.t=x;
bpsth.psths=psths_t;
bpsth.unitTrials=unitTrials;
bpsth.unitStimcond=unitStimcond;
bpsth.unitLED=unitLED;
allspikespsth=bpsth;

fractionpsth=getBurstFraction(allspikespsth,burstpsth);

end

function [varargout] = raster_sub(spikes,hAxes,bAppend,showBursts,duration,shaden)
% function [varargout] = raster(spikes,hAxes,bAppend,showBursts)
% INPUTS
%   spikes: The spikes struct. This function requires that spikes contain
%   only the following fields:
%       - .spiketimes: Vector of spike times
%       - .trials: Vector indicating trial in which each spike occurred
%   hAxes: Handles to axes;
%   bAppend: 1, append to trials already displayed; 0, use absolute trials in spikes struct 
%   showBursts: 1, plot "burst spiks" (ISIs < 4 ms)
%   duration: duration of trial
%
% OUTPUTS
%   varargout(1) = hRaster, handle to raster
%   varargout(2) = hAxes, handle to axes
%

% Created: 5/11/10 - SRO
% Modified: 11/3/11 - KR passes in trial duration

if nargin < 2
    hAxes = axes;
    bAppend = 0;
    showBursts = 0;
end

if bAppend == 0
    useTrials = 1;
else
    useTrials = 0;
end

% Set trial duration
% KR passes in

% Set spiketimes
spiketimes = spikes.spiketimes;

% Set trial numbers
if isfield(spikes,'trialsInFilter') && ~useTrials
    trials = spikes.trialsInFilter;
else
    trials = spikes.trials;
end

% Get previous trial information
t = get(hAxes,'UserData');
if ~isempty(trials)
    
    if isempty(t)
        if isfield(spikes.sweeps,'trialsInFilter')
            t.min = min(spikes.sweeps.trialsInFilter);
            t.max = max(spikes.sweeps.trialsInFilter);
        else
            t.min = min(spikes.sweeps.trials);
            t.max = max(spikes.sweeps.trials);
        end
    else
        if bAppend == 1
            % Offset by max trial number already on plot
            trials = trials + max(t.max);
        end
        t.min(end+1) = min(trials);
        t.max(end+1) = max(trials);
    end
else
    if ~isempty(t)
        if bAppend == 1
            % Offset by max trial number already on plot
            trials = 0;
            trials = trials + max(t.max);
        end
        t.min(end+1) = min(trials);
        t.max(end+1) = max(trials);
        spiketimes = NaN;
        trials = NaN;
    else
            t.min = 0;
            t.max = 1;
            spiketimes = NaN;
            trials = NaN;
    end
end
t.min=0;
t.max=10000;

% Store trial information as UserData in hAxes
set(hAxes,'UserData',t)

% Shade each trial
hold on;
if ~isempty(shaden)
    oddT=1;
    for i=trials(1):shaden:trials(end)
        if oddT==1
            fill([0 duration duration 0],[i-0.5 i-0.5 i+shaden-1+0.5 i+shaden-1+0.5],[0.9 0.9 0.9],'EdgeColor','none');
            oddT=0;
        else
            oddT=1;
        end
    end
end

% Make raster on hAxes
set(gcf,'CurrentAxes',hAxes)
hRaster = linecustommarker(spiketimes,trials);

% Set default raster properties
numtrials = length(spikes.sweeps.trials);
offset = numtrials*0.03;
offset = numtrials*0;
ymin = (min(t.min)-offset);
ymax = (max(t.max)+offset);
set(hAxes,'TickDir','out','YDir','reverse','FontSize',9, ...
    'YLim',[ymin ymax],'XLim',[0 duration])
xlabel(hAxes,'seconds')
ylabel(hAxes,'trials')

% Color synch spikes
synchtimes=spiketimes(spikes.isSynchSpike==1);
synchtrials=spikes.trials(spikes.isSynchSpike==1);
hSynch=linecustommarker(synchtimes,synchtrials);
set(hSynch,'Color',[1 0 0]);
% hSynch=line('XData',synchtimes,'YData',synchtrials,'LineStyle','none','Color',[0 1 0],'Marker','o','MarkerSize',2);

% Find bursts and plot
if showBursts
    [bTimes bTrials] = findBursts(spiketimes,trials);
    hBurst = linecustommarker(bTimes,bTrials);
    set(hBurst,'Color',[0 1 0]);
%     hBurst = line('XData',bTimes,'YData',bTrials,'LineStyle','none',...
%         'Color',[0 1 0],'Marker','o','MarkerSize',2);
else
    hBurst = [];
end
     
% Outputs
varargout{1} = hRaster;
varargout{2} = hAxes;
varargout{3} = hBurst;
end

function [firstSpikesInBurst,subsequentSpikesInBurst,preISI,postISI,isThirdSpike]=sub_findBursts(spiketimes,trialDuration)

% spikeISI=0.004;
% % spikeISI=0.006;
% preBurstISI=0.1;

% spikeISI=0.01;
% % spikeISI=0.006;
% preBurstISI=0;

% spikeISI=0.006;
% % spikeISI=0.006;
% preBurstISI=0.05;

% spikeISI=0.01;
% % spikeISI=0.006;
% preBurstISI=0.08;

% FINAL
spikeISI=0.006;
preBurstISI=0.08;

isi=[spiketimes(1) diff(spiketimes) trialDuration];
isi(isi<0)=spiketimes(isi<0);
preISI=isi(1:end-1);
postISI=isi(2:end);
% burstSpikes=isi<=spikeISI;
% precededByGapSpikes=isi>=preBurstISI;
firstSpikesInBurst=preISI>=preBurstISI & postISI<=spikeISI;


% isBurst=[precededByGapSpikes(1:end-1) & burstSpikes(2:end) logical(0)];
% firstSpikesInBurst=isBurst;

isSubsequentSpikes=zeros(1,length(spiketimes));
isThirdSpike=zeros(1,length(spiketimes));
isSubsequentSpikes([logical(0) (preISI(2:end)<=spikeISI  & firstSpikesInBurst(1:end-1))])=1;
isSubsequentSpikes([logical(0) logical(0) (preISI(3:end)<=spikeISI & preISI(2:end-1)<=spikeISI  & firstSpikesInBurst(1:end-2))])=1;
isThirdSpike([logical(0) logical(0) (preISI(3:end)<=spikeISI & preISI(2:end-1)<=spikeISI  & firstSpikesInBurst(1:end-2))])=1;
% isSubsequentSpikes([logical(0) isBurst])=1;
% stillBurst=precededByGapSpikes(1:end-2) & burstSpikes(3:end);
% isSubsequentSpikes([logical(0) logical(0) stillBurst])=1;
subsequentSpikesInBurst=isSubsequentSpikes;

end
    

function [x,psths_t,unitByUnitTrials,unitByUnitStimcond,unitByUnitLED]=getTrialByTrialUnitPSTH_sub(spikes,allAssigns,useLED,bin,trialDuration,useStimcond)

% Get trials for each unit
unitByUnitTrials=cell(1,length(allAssigns));
unitByUnitStimcond=cell(1,length(allAssigns));
unitByUnitLED=cell(1,length(allAssigns));
tt=unique(spikes.sweeps.trials);
if tt(1)>1 && length(spikes.sweeps.trials)==length(spikes.sweeps.stimcond)
    spikes.sweeps.trials=spikes.sweeps.trials-tt(1)+1;
    spikes.trials=spikes.trials-tt(1)+1;
end
% if any(isnan(spikes.sweeps.led))
%     spikes.sweeps.led(isnan(spikes.sweeps.led))=-10;
%     spikes.led(isnan(spikes.led))=-10;
% end
% if any(isnan(spikes.sweeps.stimcond))
%     spikes.sweeps.stimcond(isnan(spikes.sweeps.stimcond))=-10;
%     spikes.stimcond(isnan(spikes.stimcond))=-10;
% end

% Changed to deal with nans 7/22/15
for i=1:length(allAssigns)
    strials=spikes.sweeps.trials(~isnan(spikes.sweeps.stimcond) & ~isnan(spikes.sweeps.led));
    sstim=spikes.sweeps.stimcond(~isnan(spikes.sweeps.stimcond) & ~isnan(spikes.sweeps.led));
    sled=spikes.sweeps.led(~isnan(spikes.sweeps.stimcond) & ~isnan(spikes.sweeps.led));
    [unitByUnitTrials{i},indsin]=unique(strials);
    unitByUnitStimcond{i}=sstim(indsin);
    unitByUnitLED{i}=sled(indsin); 
end
 
% for i=1:length(allAssigns)
%     unitByUnitTrials{i}=unique(spikes.sweeps.trials);
%     unitByUnitStimcond{i}=spikes.sweeps.stimcond(unique(spikes.sweeps.trials));
%     unitByUnitLED{i}=spikes.sweeps.led(unique(spikes.sweeps.trials));
% end

psths_t=cell(length(allAssigns),length(useStimcond));
x=[];
for i=1:length(allAssigns)
    for j=1:length(useStimcond)
        useSpikes=filtspikes(spikes,0,'assigns',allAssigns(i),'stimcond',useStimcond{j});
        currTrialCon=unitByUnitTrials{i};
        currStimCon=unitByUnitStimcond{i};
        currLEDCon=unitByUnitLED{i};
        if ~isempty(useLED)
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{j}) & ismember(currLEDCon,useLED));
        else
            unitByUnitConsensus=currTrialCon(ismember(currStimCon,useStimcond{j}));
        end
        if isempty(unitByUnitConsensus)
            disp('problem with unitByUnitConsensus');
        end
        if length(unitByUnitConsensus)~=length(unitByUnitTrials{i})
            disp('stop here');
        end
        [~,~,~,x,~,psths_t{i,j}]=psth_wStd_trialByTrial(useSpikes,bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
        if size(psths_t{i,j},1)~=length(unitByUnitTrials{i})
            disp('stop');
        end
    end
end

end
    