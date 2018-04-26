function onsetdelays=measureOnsetSlowing_wrapper(spikes,normInh_tau,redInh_tau)

led=[2 4 6 8];
noled=[1 3 5 7];

s=unique(spikes.stimcond);
sToUse={1:4; 5:8 ;9:12 ;13:16; 17:20 ;21:24; 25:28; 29:32};
onsetdelays=zeros(length(s),3);
% figure(); 
k=1;
for i=1:length(s)
    disp(i);
    [x,y1,y2]=muaPSTH(filtspikes(spikes,0,'stimcond',s(i)),[3:70],led,noled,1,[],[]);
%     plot(downSampAv(x,k),downSampAv(y1,k));
%     hold on;
%     plot(downSampAv(x,k),downSampAv(y2,k),'Color','r');
    [B,A]=butter(5, 0.008, 'low');
    y1=filtfilt(B,A,y1);
    y2=filtfilt(B,A,y2);
    [~,~,~,onsetdelays(i,:)]=predictCorticalResponseGivenTau(downSampAv(normInh_tau,k),downSampAv(redInh_tau,k),downSampAv(x,k),downSampAv(y1,k),downSampAv(y2,k));
end
end

function [xpoints1,ypoints1,ypoints2]=muaPSTH(spikes,useTheseFileInds,ledValue,noLedValues,binsize,useBlackTrials,useRedTrials)
% binsize in ms

nunits=1;

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

if ~isempty(useBlackTrials)
    allUnits_noLED=filtspikes(allUnits_noLED,0,'trials',sort(useBlackTrials));
end
if ~isempty(useRedTrials)
    allUnits_withLED=filtspikes(allUnits_withLED,0,'trials',sort(useRedTrials));
end
    
[n1,centers1,edges1,xpoints1,ypoints1,stds1]=psthwStdev(allUnits_noLED,binsize);

[n2,centers2,edges2,xpoints2,ypoints2,stds2]=psthwStdev(allUnits_withLED,binsize);

ypoints1=ypoints1/nunits;
ypoints2=ypoints2/nunits;

% figure();
% plot(xpoints1,ypoints1,'Color','black');
% hold on;
% plot(xpoints2,ypoints2,'Color','red');
% hold off;
end

function [varargout] = psthwStdev(spikes,binsize,bsmooth)
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
duration = 5;     % Maybe add checking for equal sweep durations?
if isfield(spikes,'sweeps')
    a=unique(spikes.trials);
    if length(spikes.sweeps.trials)==4*length(a)
        numtrials=length(unique(spikes.trials));
    else
        numtrials = length(spikes.sweeps.trials);
    end
else
    numtrials=length(unique(spikes.trials));
end
% numtrials=length(unique(spikes.trials));
% numtrials = length(spikes.sweeps.trials);

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
    xpoints=centers(1:end-1);
    ypoints=smooth(n(1:end-1),3);
else
    xpoints=centers(1:end-1);
    ypoints=n(1:end-1);
end

% Outputs
varargout{1} = n;
varargout{2} = centers;
varargout{3} = edges;
varargout{4} = xpoints;
varargout{5} = ypoints;
varargout{6} = 1*std(nsForStdev(:,1:end-1),1);
end