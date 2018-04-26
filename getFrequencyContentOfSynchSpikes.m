function [outim,av_synch,synchAcrossTrials]=getFrequencyContentOfSynchSpikes(spikes,useTrials,randFraction,F1freq,usestimcond,usel,trialDuration)

alphaBand=[4 20];

backup_spikes=spikes;

if isempty(useTrials)
    tri=unique(spikes.trials(ismember(spikes.led,usel) & ismember(spikes.stimcond,usestimcond)));
    takeNTrials=floor(randFraction.*length(tri));
    takeTrialInds=randi(length(tri),[1,takeNTrials]);
    takeTrials=tri(takeTrialInds);
    spikes=filtspikes(spikes,0,'trials',takeTrials);
else
    spikes=filtspikes(spikes,0,'trials',useTrials);
end
 
usestimcond=unique(spikes.stimcond);
usel=unique(spikes.led);
[av_synch,synchAcrossTrials,synchSpikes]=sub_plotAlphaRasters(spikes,unique(spikes.assigns),[],0,F1freq,usestimcond,usel,trialDuration);

figure(); 
plot(linspace(0,trialDuration,length(av_synch)),av_synch,'Color','k');
hold on;
plot(linspace(0,trialDuration,length(av_synch)),av_synch+std(synchAcrossTrials,[],1)./sqrt(size(synchAcrossTrials,1)),'Color','k');
plot(linspace(0,trialDuration,length(av_synch)),av_synch-std(synchAcrossTrials,[],1)./sqrt(size(synchAcrossTrials,1)),'Color','k');
title('Number of Synch Spikes');

params.tapers=[3 10];
params.Fs=1./(trialDuration/size(synchAcrossTrials,2));
params.fpass=[1 50];
params.trialave=0;
movingwin=[1 0.05];

[S,t,f]=mtspecgrampb(synchAcrossTrials',movingwin,params);

figure();
outim.t=t;
outim.f=f;
outim.im=reshape(nanmean(S,3),size(S,1),size(S,2))';
imagesc(t,f,reshape(nanmean(S,3),size(S,1),size(S,2))');
title('Frequency Content of Synch Spikes');

meanS=reshape(nanmean(S,3),size(S,1),size(S,2));
figure(); 
plot(linspace(0,trialDuration,size(meanS,1)),nansum(meanS(:,f>=alphaBand(1) & f<=alphaBand(2)),2)./nansum(meanS,2));
title('Fraction of Total Power That Is In Alpha Band');


end

function [av_synch,synchAcrossTrials,spikes]=sub_plotAlphaRasters(backup_spikes,useAssigns,orderCells,combineCells,refFreq,usestimcond,usel,trialDuration)

spikes=backup_spikes;

% getn=3;
getn=10000;
% usel=[0];
% synchWindow=62.5; % in ms 
% synchWindow=15; % in ms
% synchWindow=166; % in ms 
synchWindow=10; % in ms 
% synchWindow=0; % in ms 
% % trialDuration=3.5; % in s
% trialDuration=14.5; % in s
% trialDuration=8; % in s
% trialDuration=2; % in s
unitsTrialsNear=1;
plotFigs=0;

% usestimcond=[1];

if (usel>0 & usel<1) | (usel>5 & usel<6)
    spikes.led=ceil(spikes.led);
    spikes.sweeps.led=ceil(spikes.sweeps.led);
    usel=ceil(usel);
end

if length(useAssigns)>getn
    [~,sortInd]=sort(orderCells(:,usestimcond));
    useUnits=useAssigns(sortInd(end-getn+1:end));
else
    useUnits=useAssigns;
end

if combineCells==1
    spikes.assigns(ismember(spikes.assigns,useUnits))=1000;
    useUnits=1000;
end

% Check for synchronous spiking on a trial-by-trial basis
synchWindow=synchWindow./1000;
% Pick out trials for this stimcond and LED cond
% allTrials=unique(spikes.sweeps.trials(ismember(spikes.sweeps.led,usel)));
allTrials=unique(spikes.sweeps.trials(ismember(spikes.sweeps.led,usel) & ismember(spikes.sweeps.stimcond,usestimcond)));
spikes.isSynchSpike=zeros(size(spikes.spiketimes));
edges=0:synchWindow:trialDuration;
synchAcrossTrials=zeros(length(allTrials),length(edges)-1); 
for i=1:length(allTrials)
    currTrial=allTrials(i);
    usespiketimes=spikes.spiketimes(ismember(spikes.assigns,useUnits) & ismember(spikes.trials,currTrial));
    spikeassigns=spikes.assigns(ismember(spikes.assigns,useUnits) & ismember(spikes.trials,currTrial));
    [usespiketimes,sortInd]=sort(usespiketimes);
    spikeassigns=spikeassigns(sortInd);
    first_isis=diff(usespiketimes);
    isis=first_isis;
    sameunits=diff(spikeassigns);
    isis(sameunits==0)=10000;
    nearbyspikes=zeros(2,length(isis)+1);
    nearbyspikes(1,1:end-1)=isis<=synchWindow;
    nearbyspikes(2,2:end)=isis<=synchWindow;
    donenearby=sum(nearbyspikes,1)>0;
    clear tempissynch
    tempissynch(sortInd)=donenearby;
    spikes.isSynchSpike(ismember(spikes.assigns,useUnits) & ismember(spikes.trials,currTrial))=tempissynch;
    for j=1:length(edges)-1
        synchAcrossTrials(i,j)=sum(spikes.isSynchSpike==1 & ismember(spikes.trials,currTrial) & spikes.spiketimes>=edges(j) & spikes.spiketimes<=edges(j+1));
    end
end

% synchAcrossTrials=zeros(length(allTrials),length(0:synchWindow:trialDuration)-1);    
% windowOffset=0:synchWindow./3:synchWindow;    
% for i=1:length(allTrials)
%     currTrial=allTrials(i);
%     usespiketimes=spikes.spiketimes(ismember(spikes.assigns,useUnits) & ismember(spikes.trials,currTrial));
%     countSynch=zeros(length(windowOffset),length(0:synchWindow:trialDuration)-1);
%     for j=1:length(windowOffset)
%         currw=windowOffset(j);
%         edges=0+currw:synchWindow:trialDuration;
%         for k=1:length(edges)-1
%             currStepWindow=[edges(k) edges(k+1)];
%             countSpikesinWindow=sum(usespiketimes>=currStepWindow(1) & usespiketimes<=currStepWindow(2));
%             countSynch(j,k)=countSpikesinWindow;
%             if countSpikesinWindow>1
%                 spikes.isSynchSpike(ismember(spikes.assigns,useUnits) & ismember(spikes.trials,currTrial) & (spikes.spiketimes>=currStepWindow(1) & spikes.spiketimes<=currStepWindow(2)))=1;
%             end
%         end
%     end
%     synchAcrossTrials(i,:)=nanmean(countSynch,1);
% end
 
av_synch=nanmean(synchAcrossTrials,1);
if plotFigs==1
    figure(); 
    plot(edges(1:end-1),nanmean(synchAcrossTrials,1),'Color','r');
    title('Synch Spikes Count');
end

% spikes.trialsInFilter=spikes.sweeps.trials(ismember(spikes.sweeps.led,usel) & ismember(spikes.sweeps.stimcond,usestimcond));
% %spikes.trials(ismember(spikes.led,usel) & ismember(spikes.stimcond,usestimcond));
% spikes.sweeps.trialsInFilter=spikes.sweeps.trials(ismember(spikes.sweeps.led,usel) & ismember(spikes.sweeps.stimcond,usestimcond));
if plotFigs==1
    ss=filtspikes(spikes,0,'led',usel,'stimcond',usestimcond);
    [~,~,ss.trials]=unique(ss.trials);
    [~,~,ss.sweeps.trials]=unique(ss.sweeps.trials);
    rmfield(ss,'trialsInFilter');
    if ~isempty(refFreq)
        figure();
        refTimepoints=0:1/refFreq:trialDuration;
        minTrial=min(ss.trials);
        maxTrial=max(ss.trials);
        if unitsTrialsNear==1
            hax(1)=subplot(2,1,1);
            pos=get(hax(1),'Position');
            set(hax(1),'Position',[pos(1) pos(2) pos(3) pos(4)*2]);
            for i=1:length(useUnits)
                currUnit=useUnits(i);
                currSpikes=filtspikes(ss,0,'assigns',currUnit);
                currSpikes.trials=currSpikes.trials.*length(useUnits)-length(useUnits)+i;
                currSpikes.sweeps.trials=currSpikes.sweeps.trials.*length(useUnits)-length(useUnits)+i;
                if i==1
                    raster_sub(currSpikes,hax(1),0,0,trialDuration,getn);
                else
                    raster_sub(currSpikes,hax(1),0,0,trialDuration,[]);
                end
            end
            minTrial=1;
            maxTrial=max(currSpikes.sweeps.trials);
            ylim([minTrial-0.5 maxTrial+0.5]);
            for j=1:length(refTimepoints)
                line([refTimepoints(j) refTimepoints(j)],[minTrial-0.5 maxTrial+0.5],'Color',[0.5 0.5 0.5]);
            end
            hax2=subplot(2,1,2);
            pos=get(hax2,'Position');
            set(hax2,'Position',[pos(1) pos(2)-0.05 pos(3) pos(4)/3]);
            pos=get(hax(1),'Position');
            set(hax(1),'Position',[pos(1) pos(2)-0.3 pos(3) pos(4)]);
            ss=filtspikes(spikes,0,'led',usel,'stimcond',usestimcond);
            ss=filtspikes(ss,0,'assigns',useUnits);
            [~,~,~,xpoints,ypoints]=psth_wStdev_valuesOnly_sub(ss,((1/refFreq)/4)*1000,0,trialDuration);
            plot(xpoints,ypoints);
            ylim([min(ypoints) max(ypoints)]);
        else
            for i=1:length(useUnits)
                hax(i)=subplot(length(useUnits)+1,1,i);
                if length(useUnits)==1
                    pos=get(hax(i),'Position');
                    set(hax(i),'Position',[pos(1) pos(2) pos(3) pos(4)*2]);
                end
                currUnit=useUnits(i);
                currSpikes=filtspikes(ss,0,'assigns',currUnit);
                raster_sub(currSpikes,hax(i),0,0,trialDuration,[]);
                ylim([minTrial-0.5 maxTrial+0.5]);
                for j=1:length(refTimepoints)
                    line([refTimepoints(j) refTimepoints(j)],[minTrial-0.5 maxTrial+0.5],'Color',[0.5 0.5 0.5]);
                end
            end
            hax2=subplot(length(useUnits)+1,1,length(useUnits)+1);
            pos=get(hax2,'Position');
            set(hax2,'Position',[pos(1) pos(2)-0.05 pos(3) pos(4)/3]);
            %     if length(useUnits)==1
            for i=1:length(hax)
                pos=get(hax(i),'Position');
                if length(useUnits)==1
                    set(hax(i),'Position',[pos(1) pos(2)-0.3 pos(3) pos(4)]);
                else
                    
                    set(hax(i),'Position',[pos(1) pos(2)-0.05 pos(3) pos(4)*1.3]);
                end
            end
            %     end
            ss=filtspikes(spikes,0,'led',usel,'stimcond',usestimcond);
            ss=filtspikes(ss,0,'assigns',useUnits);
            %         [~,~,~,xpoints,ypoints]=psth_wStdev_valuesOnly(ss,((1/refFreq)/4)*1000);
            [~,~,~,xpoints,ypoints]=psth_wStdev_valuesOnly_sub(ss,((1/refFreq)/4)*1000,0,trialDuration);
            plot(xpoints,ypoints);
            ylim([min(ypoints) max(ypoints)]);
        end
    else
        figure();
        minTrial=min(spikes.trials);
        maxTrial=minTrial+length(unique(spikes.trialsInFilter));
        for i=1:length(useUnits)
            subplot(length(useUnits)+1,1,i);
            hax=gca;
            currUnit=useUnits(i);
            currSpikes=filtspikes(spikes,0,'assigns',currUnit);
            raster_sub(currSpikes,hax,0,0,3.5,[]);
            ylim([minTrial maxTrial]);
        end
        subplot(length(useUnits)+1,1,length(useUnits)+1);
        ss=filtspikes(spikes,0,'assigns',useUnits,'led',usel);
        ss=filtspikes(spikes,0,'stimcond',usestimcond);
        [~,~,~,xpoints,ypoints]=psth_wStdev_valuesOnly(ss,20);
        plot(xpoints,ypoints);
    end
end
end

function [varargout] = psth_wStdev_valuesOnly_sub(spikes,binsize,bsmooth,duration)
% function [varargout] = psth(spikes,binsize,hAxes,bsmooth,duration)
%
% INPUTS
%   spiketimes:
%   binsize:
%   hAxes:
%   bsmooth:
%   duration: trial duration
%
% OUTPUTS
%   varargout{1} = hPsth;
%   varargout{2} = hAxes;
%   varargout{3} = n;
%   varargout{4} = centers;
%   varargout{5} = edges;
%   varargout{6} = xpoints;
%   varargout{7} = ypoints;

% Created:  3/14/10 - SRO
% Modified: 5/14/10 - SRO
%           6/8/10 - SRO
%           11/3/11 - KR passes in duration

if nargin < 2
    binsize = 50; 
end

% Use current axes if hAxes not supplied
% if nargin < 3
%     hAxes = gca;   
% end
% 
if nargin < 3
    bsmooth = 0;
end

% Set duration and number of trials
% duration = 3.5;     % Maybe add checking for equal sweep durations?
if isfield(spikes,'sweeps')
    a=unique(spikes.trials);
    if length(spikes.sweeps.trials)==4*length(a)
%     if length(spikes.sweeps.trials)>=4*length(a)-5 && length(spikes.sweeps.trials)<=4*length(a)+5
        numtrials=length(unique(spikes.trials));
    else
        numtrials = length(spikes.sweeps.trials);
    end
else
    numtrials=length(unique(spikes.trials));
end
% numtrials=length(unique(spikes.trials));
% numtrials = length(spikes.sweeps.trials);
% numtrials=119;

% Set spiketimes
spiketimes = spikes.spiketimes;

% Convert binsize from ms to s
binsize = binsize/1000;

% Get counts
edges = 0:binsize:duration;
n = histc(spiketimes,edges);
n = n/numtrials/binsize;

nsForStdev=zeros(length(unique(spikes.trials)),size(n,2));
allTrials=unique(spikes.trials);
for i=1:length(allTrials)
    cspikes=filtspikes(spikes,0,'trials',allTrials(i));
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
% hPsth = line('XData',centers(1:end-1),'YData',smooth(n(1:end-1),3),...
%     'Parent',hAxes,'LineWidth',1.5);
    xpoints=centers(1:end-1);
    ypoints=smooth(n(1:end-1),3);
else
%     hPsth = line('XData',centers(1:end-1),'YData',n(1:end-1),...
%     'Parent',hAxes,'LineWidth',1.5);
    xpoints=centers(1:end-1);
    ypoints=n(1:end-1);
end

% Set default axes properties
% if sum(n) == 0
%     maxN = 0.1;
% else
%     maxN = max(n);
% end
% axis([0 duration 0 maxN])
% set(hAxes,'TickDir','out','FontSize',8)
% xlabel(hAxes,'seconds')
% ylabel(hAxes,'spikes/s')

% Outputs
% varargout{1} = hPsth;
% varargout{2} = hAxes;
varargout{1} = n;
varargout{2} = centers;
varargout{3} = edges;
varargout{4} = xpoints;
varargout{5} = ypoints;
varargout{6} = 1*std(nsForStdev(:,1:end-1),1);

% figure();
% plot(centers(1:end-1),n(1:end-1),'Color','black');
% hold on; 
% plot(centers(1:end-1),n(1:end-1)-1*std(nsForStdev(1:end-1),1),'Color','green');
% plot(centers(1:end-1),n(1:end-1)+1*std(nsForStdev(1:end-1),1),'Color','green');
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
if any(isnan(spikes.sweeps.led))
    spikes.sweeps.led(isnan(spikes.sweeps.led))=-10;
    spikes.led(isnan(spikes.led))=-10;
end
if any(isnan(spikes.sweeps.stimcond))
    spikes.sweeps.stimcond(isnan(spikes.sweeps.stimcond))=-10;
    spikes.stimcond(isnan(spikes.stimcond))=-10;
end

for i=1:length(allAssigns)
    unitByUnitTrials{i}=unique(spikes.sweeps.trials);
    unitByUnitStimcond{i}=spikes.sweeps.stimcond(unique(spikes.sweeps.trials));
    unitByUnitLED{i}=spikes.sweeps.led(unique(spikes.sweeps.trials));
end

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
        [~,~,~,x,~,psths_t{i,j}]=psth_wStd_trialByTrial(useSpikes,bin,0,trialDuration,length(unitByUnitConsensus),unitByUnitConsensus);
    end
end

end