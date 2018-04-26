function [mean_noled,mean_led,std_noled,std_led]=bootstrapByTrialsBETTER(spikes,useTheseFileInds)

ledValue=[2.02 4.04 6.06 8.08]; % red
noLedValue=[1.01 3.03 5.05 7.07]; % black
fractionOfTrials=0.8;
normalize=0;
% normwindow1=[0.9 1.05];
% normwindow2=[1.26 1.298];
% normwindow1=[0.95 1.1];
% normwindow1=[1.02 1.07];
% normwindow2=[1.3 1.338];
normwindow1=[1.02 1.07];
normwindow2=[1.2 1.3];
nRuns=2000;
binsize=1;
duration=5;

temp=[];
temp1=[];
for i=1:length(ledValue)
    subLed=makeTempField(spikes,'led',ledValue(i));
    temp(i,:)=subLed.temp;
    temp1(i,:)=subLed.sweeps.temp;
    subLed.temp=[];
    subLed.sweeps.temp=[];
end
subLed.temp=sum(temp,1)>=1;
subLed.sweeps.temp=sum(temp1,1)>=1;
subLed=filtspikes(subLed,0,'temp',1);
% subLed=filtspikes(spikes,0,'led',ledValue);

temp=[];
temp1=[];
for i=1:length(noLedValue)
    subNoLed=makeTempField(spikes,'led',noLedValue(i));
    temp(i,:)=subNoLed.temp;
    temp1(i,:)=subNoLed.sweeps.temp;
    subNoLed.temp=[];
    subNoLed.sweeps.temp=[];
end
subNoLed.temp=sum(temp,1)>=1;
subNoLed.sweeps.temp=sum(temp1,1)>=1;
subNoLed=filtspikes(subNoLed,0,'temp',1);
% subNoLed=filtspikes(spikes,0,'led',noLedValue);
nTrialsLed=floor(fractionOfTrials*length(unique(subLed.trials)));
nTrialsNoLed=floor(fractionOfTrials*length(unique(subNoLed.trials)));

spiketimes=spikes.spiketimes;
andTrials=spikes.trials;
[xpoints,ypoints1_all]=psth_onlySpiketimes(spiketimes,binsize,0,duration,length(unique(spikes.trials)));
noLedPSTH=zeros(nRuns,length(ypoints1_all));
ledPSTH=zeros(nRuns,length(ypoints1_all));
noLedTrials=unique(subNoLed.trials);
ledTrials=unique(subLed.trials);
for i=1:nRuns
    if mod(i,100)==1
        disp(i);
    end
    a=randsample(length(noLedTrials),nTrialsNoLed,1);
    useNoLedTrials=noLedTrials(a);
    sumPSTHs=zeros(1,length(ypoints1_all));
    for j=1:length(useNoLedTrials)
        currspiketimes=spiketimes(andTrials==useNoLedTrials(j));
        [x,y]=psth_onlySpiketimes(currspiketimes,binsize,0,duration,1);
        sumPSTHs=sumPSTHs+y;
    end
    noLedPSTH(i,:)=sumPSTHs./length(useNoLedTrials);
    
    b=randsample(length(ledTrials),nTrialsLed,1);
    useLedTrials=ledTrials(a);
    sumPSTHs=zeros(1,length(ypoints1_all));
    for j=1:length(useLedTrials)
        currspiketimes=spiketimes(andTrials==useLedTrials(j));
        [x,y]=psth_onlySpiketimes(currspiketimes,binsize,0,duration,1);
        sumPSTHs=sumPSTHs+y;
    end
    ledPSTH(i,:)=sumPSTHs./length(useLedTrials);
end

if normalize==1
    m_noled=mean(noLedPSTH,1);
    s_noled=std(noLedPSTH,[],1);
    m_led=mean(ledPSTH,1);
    s_led=std(ledPSTH,[],1);
    % Remove photoartifact
%     m_noled(1339:1343)=repmat(m_noled(1339),1,length(1339:1343));
%     m_led(1339:1343)=repmat(m_led(1339),1,length(1339:1343));
%     m_noled(1297:1304)=repmat(m_noled(1297),1,length(1297:1304));
%     m_led(1297:1304)=repmat(m_led(1297),1,length(1297:1304));
    base_noled=mean(m_noled(xpoints>normwindow1(1) & xpoints<normwindow1(2)));
    base_led=mean(m_led(xpoints>normwindow1(1) & xpoints<normwindow1(2)));
    peak_noled=mean(m_noled(xpoints>normwindow2(1) & xpoints<normwindow2(2)));
    peak_led=mean(m_led(xpoints>normwindow2(1) & xpoints<normwindow2(2)));
    m_noled=m_noled-base_noled;
    m_led=m_led-base_led;
    m_noled=m_noled./(peak_noled-base_noled);
    m_led=m_led./(peak_led-base_led);
    s_noled=s_noled./(peak_noled-base_noled);
    s_led=s_led./(peak_led-base_led);
    figure(); 
    hax=axes();
    hl=plot(xpoints,m_noled,'Color','k');
%     addErrBar(xpoints,m_noled,s_noled,'y',hax,hl);
    hold on;
    hl=plot(xpoints,m_led,'Color','r');
%     addErrBar(xpoints,m_led,s_led,'y',hax,hl);
else
    m_noled=mean(noLedPSTH,1);
    s_noled=std(noLedPSTH,[],1);
    m_led=mean(ledPSTH,1);
    s_led=std(ledPSTH,[],1);
    % Remove photoartifact
%     m_noled(1339:1343)=repmat(m_noled(1339),1,length(1339:1343));
%     m_led(1339:1343)=repmat(m_led(1339),1,length(1339:1343));
    figure(); 
    hax=axes();
    hl=plot(xpoints,mean(noLedPSTH,1),'Color','k');
%     addErrBar(xpoints,mean(noLedPSTH,1),std(noLedPSTH,[],1),'y',hax,hl);
    hold on;
    hl=plot(xpoints,mean(ledPSTH,1),'Color','r');
%     addErrBar(xpoints,mean(ledPSTH,1),std(ledPSTH,[],1),'y',hax,hl);
end

mean_noled=m_noled;
mean_led=m_led;
std_noled=s_noled;
std_led=s_led;
end

function [xpoints,ypoints]=psth_onlySpiketimes(spiketimes,binsize,bsmooth,duration,numtrials)

if nargin < 2
    binsize = 50; 
end

if nargin < 3
    bsmooth = 0;
end

% Convert binsize from ms to s
binsize = binsize/1000;

% Get counts
edges = 0:binsize:duration;
n = histc(spiketimes,edges);
n = n/numtrials/binsize;

if all(isnan(n))
    n = 0;
end

% Compute center of bins
centers = edges + diff(edges(1:2))/2;
if bsmooth
    xpoints=centers(1:end-1);
    ypoints=smooth(n(1:end-1),3);
else
    xpoints=centers(1:end-1);
    ypoints=n(1:end-1);
end
end