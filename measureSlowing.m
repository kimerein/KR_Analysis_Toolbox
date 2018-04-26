function [f,timeDelays,allPxy]=measureSlowing(spikes)

ledNorm=[1 3 5 7];
ledRed=[2 4 6 8];
ledOn=[0.51 2.99]+0.04;
stimcond=1:32;
spontToSubtract=[0.6 1]+0.04;
peak=[1.2 1.29]+0.04;
smoothBin=100;
nDerivBins=50;

xs=linspace(0,5,5000);
derivs=zeros(length(stimcond),length(xs(xs>=ledOn(1) & xs<=ledOn(2))));
heightDiffs=zeros(length(stimcond),length(xs(xs>=ledOn(1) & xs<=ledOn(2))));
timeDelays=zeros(513,length(stimcond));
allPxy=zeros(513,length(stimcond));
dontuse=zeros(1,length(stimcond));

for i=1:length(stimcond)
    [x,yNorm,yRed]=muaPSTH(filtspikes(spikes,0,'stimcond',stimcond(i)),[],ledRed,ledNorm,1,[],[]);
    if isempty(yNorm) || isempty(yRed)
        dontuse(i)=1;
        continue
    end
    yNorm=yNorm-mean(yNorm(x>=spontToSubtract(1) & x<=spontToSubtract(2)));
    m=mean(yNorm(x>=peak(1) & x<=peak(2)));
    yNorm=yNorm./m;
    yRed=yRed-mean(yRed(x>=spontToSubtract(1) & x<=spontToSubtract(2)));
    m=mean(yRed(x>=peak(1) & x<=peak(2)));
    yRed=yRed./m;
%     temp=diff(smooth(yNorm(x>=0 & x<=5),smoothBin));
    [B,A]=butter(5, 0.025, 'low');
    xFiltered=filtfilt(B,A,yNorm);
    xGradient=gradient(xFiltered);
    temp=xGradient;
%     temp=diff(yNorm(x>=0 & x<=5));
    temp=temp';

    xFiltered2=filtfilt(B,A,yRed);
    xGradient2=gradient(xFiltered2);
    temp2=xGradient2;
    temp2=temp2';
    derivs(i,:)=mean([temp(x(1:4999)>=ledOn(1) & x(1:4999)<=ledOn(2)) temp2(x(1:4999)>=ledOn(1) & x(1:4999)<=ledOn(2))],2)';
  
%     temp=yRed(x>=0 & x<=5)-yNorm(x>=0 & x<=5);
    temp=xFiltered2-xFiltered;
    heightDiffs(i,:)=temp(x(1:4999)>=ledOn(1) & x(1:4999)<=ledOn(2));
    
    [Pxy,f]=cpsd(yNorm(x>=ledOn(1) & x<=ledOn(2)),yRed(x>=ledOn(1) & x<=ledOn(2)),[],[],[],1000);
    phase=angle(Pxy);
    timeDelays(:,i)=phase./(2*pi.*f);
    allPxy(:,i)=Pxy;
end
mi=min(derivs(1:end));
ma=max(derivs(1:end));
derivBins=linspace(mi,ma+0.00001,nDerivBins);
heightsForDerivs=zeros(1,length(derivBins)-1);
stdHeights=zeros(1,length(derivBins)-1);
for i=1:length(derivBins)-1
    c=heightDiffs(derivs>=derivBins(i) & derivs<derivBins(i+1));
    heightsForDerivs(i)=mean(c);
    stdHeights(i)=std(c,[],1);
end

figure(); 
hax=axes();
hl=plot(mean([derivBins(1:end-1); derivBins(2:end)],1),heightsForDerivs);
addErrBar(mean([derivBins(1:end-1); derivBins(2:end)],1),heightsForDerivs,stdHeights,'y',hax,hl);

figure();
plot(f,mean(timeDelays(:,dontuse~=1),2));

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