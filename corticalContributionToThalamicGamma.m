function corticalContributionToThalamicGamma(r,ledVal,noLedVal)

% filez=38:91;
% filez=38:91;
% filez=69:89;
% filez=82:124;
% filez=41:52;
filez=20:47;

if isfield(r,'running')
    t=unique(r.trials);
    runningTrials=nan(size(t));
    for i=1:length(t)
        if any(r.running(r.trials==t(i)))
            runningTrials(i)=1;
        else
            runningTrials(i)=0;
        end
    end
    disp(sum(runningTrials));
end

% subS=filtspikes(r,0,'fileInd',filez,'trials',t(runningTrials==1));
subS=filtspikes(r,0,'fileInd',filez);
[subS_led,subS_noLed]=sortSpikesByLED(subS,ledVal,noLedVal);
[x,pmu,avSpec_noLED,f,ps]=gammaInThalamus_fromPSTH(subS_noLed,[]);
% [x,pmu,avSpec_noLED,f,ps]=gammaInThalamus_fromPSTH(subS,[]);
[x,pmu,avSpec_LED,f,ps_led]=gammaInThalamus_fromPSTH(subS_led,[]);
% [x,pmu,avSpec_LED,f,ps_led]=gammaInThalamus_fromPSTH(subS,[]);
% avSpec_LED=avSpec_noLED; ps_led=ps;
figure(); 
plot(f,avSpec_noLED,'Color','k');
hold on; 
plot(f,avSpec_LED,'Color','r');
figure(); 
hax=axes();
plotLineAndErr(f,ps,'k',hax);
% plotLineAndErr(f,ps.*repmat(f,size(ps,1),1),'k',hax);
hold on; 
% plotLineAndErr(f,ps_led.*repmat(f,size(ps,1),1),'r',hax);
plotLineAndErr(f,ps_led,'r',hax);

% highRange=[30 60];
% lowRange=[1 30];

lowRange=[1 50];
highRange=[50 70];
% lowRange=[1 50];
% lowRange=[80 100];
% highRange=[60 65];
lowPow=nanmean(ps(:,f>=lowRange(1) & f<=lowRange(2)),2);
highPow=nanmean(ps(:,f>=highRange(1) & f<=highRange(2)),2);
lowPow_led=nanmean(ps_led(:,f>=lowRange(1) & f<=lowRange(2)),2);
highPow_led=nanmean(ps_led(:,f>=highRange(1) & f<=highRange(2)),2);
figure();
scatter([ones(1,length(lowPow)) 2.*ones(1,length(lowPow_led))],[highPow./lowPow; highPow_led./lowPow_led]);
xlim([0 3]);
figure();
hax=axes();
plotLineAndErr(1,highPow./lowPow,'k',hax);
hold on;
plotLineAndErr(2,highPow_led./lowPow_led,'r',hax);
xlim([0 3]);

% ps=ps./repmat(max(ps,[],2),1,size(ps,2));
% ps_led=ps_led./repmat(max(ps_led,[],2),1,size(ps_led,2));
ps=ps-repmat(min(ps(:,f>=45 & f<=85),[],2),1,size(ps,2));
ps_led=ps_led-repmat(min(ps_led(:,f>=45 & f<=85),[],2),1,size(ps_led,2));
ps=ps./repmat(max(ps(:,f>=45 & f<=85),[],2),1,size(ps,2));
ps_led=ps_led./repmat(max(ps_led(:,f>=45 & f<=85),[],2),1,size(ps_led,2));
figure(); 
hax=axes();
plotLineAndErr(f,ps.*repmat(f,size(ps,1),1),'k',hax);
hold on; 
plotLineAndErr(f,ps_led.*repmat(f,size(ps_led,1),1),'r',hax);


% subS=filtspikes(r,0,'fileInd',filez,'trials',t(runningTrials==0));
% [subS_led,subS_noLed]=sortSpikesByLED(subS,ledVal,noLedVal);
% [x,pmu,avSpec_noLED,f,ps]=gammaInThalamus_fromPSTH(subS_noLed,[]);
% [x,pmu,avSpec_LED,f,ps_led]=gammaInThalamus_fromPSTH(subS_led,[]);
% figure(); 
% plot(f,avSpec_noLED,'Color','k');
% hold on; 
% plot(f,avSpec_LED,'Color','r');
% figure(); 
% hax=axes();
% plotLineAndErr(f,ps,'k',hax);
% hold on; 
% plotLineAndErr(f,ps_led,'r',hax);

end

function hl=plotLineAndErr(x,data,c,hax)

hl=plot(x,nanmean(data,1),'Color',c);
addErrBar(x,nanmean(data,1),nanstd(data,[],1)./sqrt(size(data,1)),'y',hax,hl);

end

function [allUnits_withLED,allUnits_noLED]=sortSpikesByLED(spikes,ledValue,noLedValues)
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