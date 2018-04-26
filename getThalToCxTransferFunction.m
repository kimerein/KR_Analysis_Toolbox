function [allPoly,freqTransfer,allTimepoints]=getThalToCxTransferFunction(thalSpikes,cxSpikes)

nBins=20;
freqs=[1 2 4 6 8 10 12 14 16 18 20 30 40 50 60];
ledForFreqs=freqs+0.05;
freqTransfer=cell(1,15);
allThal=zeros(15,5000);
allCx=zeros(15,5000);
allPoly=zeros(15,2);
for i=1:2:length(freqs)
    disp(i);
    % Freq 1 
    if i==length(freqs)
        ledVals.noLedValues=ledForFreqs(i);
        ledVals.ledValues=ledForFreqs(i);
    else
        ledVals.noLedValues=ledForFreqs(i);
        ledVals.ledValues=ledForFreqs(i+1);
    end
    [thalx,thaly1,thaly2]=compareMUA(thalSpikes,[],[],[],ledVals);
    [cxx,cxy1,cxy2]=compareMUA(cxSpikes,[],[],[],ledVals);
    p1=polyfit(thaly1,cxy1,1);
    p2=polyfit(thaly2,cxy2,1);
    mi=min(thaly1);
    ma=max(thaly1)+0.0001;
    useBins=linspace(mi,ma,nBins+1);
    cxSteps=zeros(1,nBins);
    thSteps=zeros(1,nBins);
    for j=1:length(useBins)-1
        if isempty(thaly1(thaly1>=useBins(j) & thaly1<useBins(j+1)))
            if j==1
                continue
            else
                thSteps(j)=mean([useBins(j) useBins(j+1)]);
                cxSteps(j)=cxSteps(j-1);
            end
        else
            thSteps(j)=mean(thaly1(thaly1>=useBins(j) & thaly1<useBins(j+1)));
        end
        if isempty(cxy1(thaly1>=useBins(j) & thaly1<useBins(j+1)))
            if j==1
                continue
            else
                cxSteps(j)=cxSteps(j-1);
            end
        else
            cxSteps(j)=mean(cxy1(thaly1>=useBins(j) & thaly1<useBins(j+1)));
        end
    end
    freqTransfer{i}=[thSteps; cxSteps];
    if i~=length(freqs)
        mi=min(thaly2);
        ma=max(thaly2)+0.0001;
        useBins=linspace(mi,ma,nBins+1);
        cxSteps=zeros(1,nBins);
        thSteps=zeros(1,nBins);
        for j=1:length(useBins)-1
            if isempty(thaly2(thaly2>=useBins(j) & thaly2<useBins(j+1)))
                if j==1
                    continue
                else
                    thSteps(j)=mean([useBins(j) useBins(j+1)]);
                    cxSteps(j)=cxSteps(j-1);
                end
            else
                thSteps(j)=mean(thaly2(thaly2>=useBins(j) & thaly2<useBins(j+1)));
            end
            if isempty(cxy2(thaly2>=useBins(j) & thaly2<useBins(j+1)))
                if j==1
                    continue
                else
                    cxSteps(j)=cxSteps(j-1);
                end
            else
                cxSteps(j)=mean(cxy2(thaly2>=useBins(j) & thaly2<useBins(j+1)));
            end
        end
        freqTransfer{i+1}=[thSteps; cxSteps];
    end
    if i==length(freqs)
        allThal(i,:)=thaly1;
        allCx(i,:)=cxy1;
        allPoly(i,:)=p1;
    else
        allThal(i,:)=thaly1;
        allThal(i+1,:)=thaly2;
        allCx(i,:)=cxy1;
        allCx(i+1,:)=cxy2;
        allPoly(i,:)=p1;
        allPoly(i+1,:)=p2;
    end
end
allTimepoints.allThal=allThal;
allTimepoints.allCx=allCx;
allTimepoints.x=cxx;

x=0:500; figure();
for i=1:size(allPoly,1)
y=allPoly(i,1).*x+allPoly(i,2);
plot(x,y);
hold on;
end
for i=1:length(freqTransfer)
curr=freqTransfer{i};
plot(curr(1,:),curr(2,:),'Color','r');
end

end

function [xpoints1,ypoints1,ypoints2]=compareMUA(spikes,useTheseFileInds,useBlackTrials,useRedTrials,ledVals)

ledValue=ledVals.ledValues; % red
noLedValues=ledVals.noLedValues; % black
binsize=1; % ms

stimconds=1:10000;
nunits=1;

spikes=filtspikes(spikes,0,'stimcond',stimconds); 

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
    
[n1,centers1,edges1,xpoints1,ypoints1,stds1]=psth_wStdev_valuesOnly(allUnits_noLED,binsize);

[n2,centers2,edges2,xpoints2,ypoints2,stds2]=psth_wStdev_valuesOnly(allUnits_withLED,binsize);

ypoints1=ypoints1/nunits;
ypoints2=ypoints2/nunits;

% figure();
% plot(xpoints1,ypoints1,'Color','black');
% hold on;
% plot(xpoints2,ypoints2,'Color','red');
% hold off;
end